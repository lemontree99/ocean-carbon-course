#!/usr/bin/env python
"""计算1982–2014年CMIP6与Landschützer海气CO2通量空间分布。"""

from __future__ import annotations

import argparse
import calendar
import json
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr


MODELS = ["GFDL-ESM4", "MPI-ESM1-2-LR", "CanESM5"]
MOLAR_MASS_C = 12.0107  # g mol-1
R_EARTH = 6_371_000.0

parser = argparse.ArgumentParser()
parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[2])
parser.add_argument("--output-dir", type=Path, default=None)
args = parser.parse_args()

root = args.project_root.resolve()
raw = root / "data" / "raw"
result_dir = (args.output_dir or root / "results" / "02_spatial_distribution").resolve()
tables, intermediate = result_dir / "tables", result_dir / "intermediate"
tables.mkdir(parents=True, exist_ok=True)
intermediate.mkdir(parents=True, exist_ok=True)
time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)


def record_seconds(ds: xr.Dataset) -> np.ndarray:
    """取得每个记录代表的秒数。"""
    name = ds.time.attrs.get("bounds")
    if name in ds:
        seconds = []
        for start, end in ds[name].values:
            delta = end - start
            seconds.append(
                float(delta / np.timedelta64(1, "s"))
                if isinstance(delta, np.timedelta64)
                else float(delta.total_seconds())
            )
        return np.asarray(seconds)
    years = np.array([time.year for time in ds.time.values])
    cal = str(ds.time.encoding.get("calendar", "standard")).lower()
    days = np.full(len(years), 365) if cal in {"noleap", "365_day"} else np.array(
        [366 if calendar.isleap(int(year)) else 365 for year in years]
    )
    return days * 86400.0


def model_mean_field(path: Path, start_year: int | None, end_year: int | None):
    """把kg C m-2 s-1积分成多年平均mol C m-2 yr-1。"""
    with xr.open_dataset(path, decode_times=time_coder) as ds:
        years = np.array([time.year for time in ds.time.values])
        keep = np.ones(len(years), dtype=bool)
        if start_year is not None:
            keep &= years >= start_year
        if end_year is not None:
            keep &= years <= end_year
        indices = np.where(keep)[0]
        seconds = record_seconds(ds)[indices]
        selected_years = years[indices]

        shape = tuple(ds.fgco2.sizes[dim] for dim in ds.fgco2.dims if dim != "time")
        total_kg_m2 = np.zeros(shape)
        valid = np.zeros(shape, dtype=bool)
        for start in range(0, len(indices), 120):
            part = indices[start : start + 120]
            values = ds.fgco2.isel(time=part).values.astype("float64")
            total_kg_m2 += np.nansum(values * seconds[start : start + len(part), None, None], axis=0)
            valid |= np.any(np.isfinite(values), axis=0)

        field = total_kg_m2 * 1000.0 / MOLAR_MASS_C / len(np.unique(selected_years))
        field[~valid] = np.nan
        if "lat" in ds.coords:
            lon, lat = np.meshgrid(ds.lon.values, ds.lat.values)
        else:
            lon, lat = ds.longitude.values, ds.latitude.values
    return field, lon % 360.0, lat


def conservative_bin(field, lon, lat, area, step=2.0):
    """按源格点面积守恒汇总到规则2°格点。"""
    nlat, nlon = int(180 / step), int(360 / step)
    valid = np.isfinite(field) & np.isfinite(lon) & np.isfinite(lat) & np.isfinite(area) & (area > 0)
    iy = np.clip(((lat[valid] + 90) // step).astype(int), 0, nlat - 1)
    ix = np.clip(((lon[valid] % 360) // step).astype(int), 0, nlon - 1)
    index = iy * nlon + ix
    area_sum = np.bincount(index, weights=area[valid], minlength=nlat * nlon)
    flux_sum = np.bincount(index, weights=field[valid] * area[valid], minlength=nlat * nlon)
    result = np.full(nlat * nlon, np.nan)
    result[area_sum > 0] = flux_sum[area_sum > 0] / area_sum[area_sum > 0]
    return result.reshape(nlat, nlon), area_sum.reshape(nlat, nlon)


def spatial_metrics(field_set, obs_area, regions):
    """计算共同格点上的全球与分区统计。"""
    global_rows, regional_rows = [], []
    for model in MODELS:
        model_field, obs_field = field_set[model], field_set["Landschuetzer"]
        common = np.isfinite(model_field) & np.isfinite(obs_field) & (obs_area > 0)
        weight = obs_area[common]
        model_values, obs_values = model_field[common], obs_field[common]
        model_mean = np.average(model_values, weights=weight)
        obs_mean = np.average(obs_values, weights=weight)
        covariance = np.average((model_values - model_mean) * (obs_values - obs_mean), weights=weight)
        correlation = covariance / np.sqrt(
            np.average((model_values - model_mean) ** 2, weights=weight)
            * np.average((obs_values - obs_mean) ** 2, weights=weight)
        )
        global_rows.append(
            {
                "model": model,
                "period": "1982-2014",
                "model_area_mean_mol_m2_yr": model_mean,
                "obs_area_mean_mol_m2_yr": obs_mean,
                "mean_bias_mol_m2_yr": np.average(model_values - obs_values, weights=weight),
                "RMSE_mol_m2_yr": np.sqrt(np.average((model_values - obs_values) ** 2, weights=weight)),
                "pattern_correlation": correlation,
                "model_integral_PgC_yr": np.sum(model_values * weight) * MOLAR_MASS_C / 1e15,
                "obs_integral_PgC_yr": np.sum(obs_values * weight) * MOLAR_MASS_C / 1e15,
            }
        )

        for region, region_mask in regions.items():
            mask = common & region_mask
            weight = obs_area[mask]
            m, o = model_field[mask], obs_field[mask]
            m_mean, o_mean = np.average(m, weights=weight), np.average(o, weights=weight)
            covariance = np.average((m - m_mean) * (o - o_mean), weights=weight)
            region_correlation = covariance / np.sqrt(
                np.average((m - m_mean) ** 2, weights=weight)
                * np.average((o - o_mean) ** 2, weights=weight)
            )
            regional_rows.append(
                {
                    "model": model,
                    "region": region,
                    "model_mean_mol_m2_yr": m_mean,
                    "obs_mean_mol_m2_yr": o_mean,
                    "bias_mol_m2_yr": np.average(m - o, weights=weight),
                    "RMSE_mol_m2_yr": np.sqrt(np.average((m - o) ** 2, weights=weight)),
                    "pattern_correlation": region_correlation,
                    "bias_integral_PgC_yr": np.sum((m - o) * weight) * MOLAR_MASS_C / 1e15,
                }
            )
    return global_rows, regional_rows


# 1. Landschützer：原始月场、向上为正，转为向下为正；不做任何漂移订正或时间平滑。
obs_file = root / "data" / "obs" / "landschutzer" / "spco2_1982-2015_MPI_SOM-FFN_v2016.nc"
with xr.open_dataset(obs_file, decode_times=time_coder, mask_and_scale=False) as ds:
    years = np.array([time.year for time in ds.time.values])
    indices = np.where((years >= 1982) & (years <= 2014))[0]
    seconds = record_seconds(ds)[indices]
    numerator = np.zeros((ds.sizes["lat"], ds.sizes["lon"]))
    denominator = np.zeros_like(numerator)
    for start in range(0, len(indices), 120):
        part = indices[start : start + 120]
        values = ds.fgco2_raw.isel(time=part).values.astype("float64")
        good = np.isfinite(values) & (np.abs(values) < 100)
        weight = seconds[start : start + len(part), None, None]
        numerator += np.sum(np.where(good, values, 0.0) * weight, axis=0)
        denominator += np.sum(good * weight, axis=0)
    obs_native = -np.divide(numerator, denominator, out=np.full_like(numerator, np.nan), where=denominator > 0)
    obs_lon, obs_lat = np.meshgrid(ds.lon.values % 360.0, ds.lat.values)

lat_lower = np.deg2rad(np.clip(obs_lat - 0.5, -90, 90))
lat_upper = np.deg2rad(np.clip(obs_lat + 0.5, -90, 90))
obs_area_native = R_EARTH**2 * np.deg2rad(1.0) * (np.sin(lat_upper) - np.sin(lat_lower))
obs_2deg, obs_area_2deg = conservative_bin(obs_native, obs_lon, obs_lat, obs_area_native)

# 2. 计算模式原始总净通量；漂移检验和局地订正在脚本05中完成。
fields = {"Landschuetzer": obs_2deg}
mapped_areas = {"Landschuetzer": obs_area_2deg}

for model in MODELS:
    with xr.open_dataset(raw / model / "historical" / "areacello.nc") as ds:
        area = ds.areacello.values.astype("float64")

    historical, lon, lat = model_mean_field(
        raw / model / "historical" / "fgco2.nc", 1982, 2014
    )
    fields[model], mapped_areas[model] = conservative_bin(historical, lon, lat, area)

# 3. 共同2°格点上的全球和分区指标。
lat_common = np.arange(-89, 90, 2.0)
lon_common = np.arange(1, 360, 2.0)
lon2, lat2 = np.meshgrid(lon_common, lat_common)
regions = {
    "Southern Ocean (60S-40S)": (lat2 >= -60) & (lat2 < -40),
    "Equatorial Pacific": (lat2 >= -10) & (lat2 <= 10) & (lon2 >= 140) & (lon2 < 280),
    "North Atlantic": (lat2 >= 30) & (lat2 < 60) & ((lon2 >= 280) | (lon2 < 20)),
    "North Pacific": (lat2 >= 30) & (lat2 < 60) & (lon2 >= 120) & (lon2 < 240),
    "Indian Ocean": (lat2 >= -30) & (lat2 < 20) & (lon2 >= 20) & (lon2 < 120),
}

global_rows, regional_rows = spatial_metrics(fields, obs_area_2deg, regions)

# 4. 保存供绘图和复核使用的场与统计表。
np.savez_compressed(
    intermediate / "spatial_fields_2deg.npz",
    lat=lat_common,
    lon=lon_common,
    obs_area=obs_area_2deg,
    **fields,
)
pd.DataFrame(global_rows).to_csv(tables / "global_spatial_metrics.csv", index=False)
pd.DataFrame(regional_rows).to_csv(tables / "regional_spatial_metrics.csv", index=False)
with (tables / "calculation_metadata.json").open("w", encoding="utf-8") as file:
    json.dump(
        {
            "period": "1982-2014",
            "sign": "positive downward into ocean",
            "observation": "Landschutzer v2016 fgco2_raw, sign multiplied by -1",
            "model_piControl_correction": (
                "not applied here; script 05 tests native-grid piControl trends and "
                "applies local correction only where FDR q<0.05 and |drift|>=0.10"
            ),
            "regridding": "source-area-weighted conservative binning to 2 degree cells",
            "smoothing": "none",
        },
        file,
        ensure_ascii=False,
        indent=2,
    )

print(pd.DataFrame(global_rows).to_string(index=False))
print(f"\n结果目录：{result_dir}")
