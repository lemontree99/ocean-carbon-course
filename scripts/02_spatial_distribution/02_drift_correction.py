#!/usr/bin/env python
"""检验三个CMIP6模式分支对齐piControl段的局地fgco2时间漂移。"""

from __future__ import annotations

import argparse
import calendar
from pathlib import Path

import cftime
import numpy as np
import pandas as pd
from scipy.stats import norm
import xarray as xr


MODELS = ["GFDL-ESM4", "MPI-ESM1-2-LR", "CanESM5"]
MOLAR_MASS_C = 12.0107
PRACTICAL_LIMIT = 0.10  # mol C m-2 yr-1，1982–2014平均场累计漂移

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


def record_seconds(ds):
    """取得每个年记录代表的秒数。"""
    name = ds.time.attrs.get("bounds")
    if name in ds:
        return np.asarray([(end - start).total_seconds() for start, end in ds[name].values])
    years = np.array([time.year for time in ds.time.values])
    cal = str(ds.time.encoding.get("calendar", "standard")).lower()
    days = np.full(len(years), 365) if cal in {"noleap", "365_day"} else np.array(
        [366 if calendar.isleap(int(year)) else 365 for year in years]
    )
    return days * 86400.0


def branch_year(path):
    """把historical分支时刻转换为piControl日历年。"""
    with xr.open_dataset(path, decode_times=False) as ds:
        date = cftime.num2date(
            ds.attrs["branch_time_in_parent"],
            ds.attrs["parent_time_units"],
            calendar=ds.time.attrs.get("calendar", "standard"),
        )
    return int(date.year)


def fdr_mask(p, q=0.05):
    """Benjamini–Hochberg FDR检验。"""
    mask = np.zeros(p.shape, dtype=bool)
    valid = np.isfinite(p)
    ordered = np.sort(p[valid])
    passed = ordered <= q * np.arange(1, len(ordered) + 1) / len(ordered)
    cutoff = ordered[np.where(passed)[0][-1]] if np.any(passed) else np.nan
    if np.isfinite(cutoff):
        mask[valid] = p[valid] <= cutoff
    return mask, cutoff


def conservative_bin(field, lon, lat, area, step=2.0):
    """按源格点面积汇总到规则2°格点。"""
    nlat, nlon = int(180 / step), int(360 / step)
    valid = np.isfinite(field) & np.isfinite(lon) & np.isfinite(lat) & np.isfinite(area) & (area > 0)
    iy = np.clip(((lat[valid] + 90) // step).astype(int), 0, nlat - 1)
    ix = np.clip(((lon[valid] % 360) // step).astype(int), 0, nlon - 1)
    index = iy * nlon + ix
    area_sum = np.bincount(index, weights=area[valid], minlength=nlat * nlon)
    value_sum = np.bincount(index, weights=field[valid] * area[valid], minlength=nlat * nlon)
    result = np.full(nlat * nlon, np.nan)
    result[area_sum > 0] = value_sum[area_sum > 0] / area_sum[area_sum > 0]
    return result.reshape(nlat, nlon)


def spatial_metrics(field_set, obs_area, regions):
    """重新计算订正后的全球和分区空间指标。"""
    global_rows, regional_rows = [], []
    for model in MODELS:
        model_field, obs_field = field_set[model], field_set["Landschuetzer"]
        common = np.isfinite(model_field) & np.isfinite(obs_field) & (obs_area > 0)
        weight = obs_area[common]
        m, o = model_field[common], obs_field[common]
        m_mean, o_mean = np.average(m, weights=weight), np.average(o, weights=weight)
        correlation = np.average((m - m_mean) * (o - o_mean), weights=weight) / np.sqrt(
            np.average((m - m_mean) ** 2, weights=weight)
            * np.average((o - o_mean) ** 2, weights=weight)
        )
        global_rows.append(
            {
                "model": model,
                "period": "1982-2014",
                "model_area_mean_mol_m2_yr": m_mean,
                "obs_area_mean_mol_m2_yr": o_mean,
                "mean_bias_mol_m2_yr": np.average(m - o, weights=weight),
                "RMSE_mol_m2_yr": np.sqrt(np.average((m - o) ** 2, weights=weight)),
                "pattern_correlation": correlation,
                "model_integral_PgC_yr": np.sum(m * weight) * MOLAR_MASS_C / 1e15,
                "obs_integral_PgC_yr": np.sum(o * weight) * MOLAR_MASS_C / 1e15,
            }
        )
        for region, region_mask in regions.items():
            mask = common & region_mask
            weight = obs_area[mask]
            m, o = model_field[mask], obs_field[mask]
            m_mean, o_mean = np.average(m, weights=weight), np.average(o, weights=weight)
            correlation = np.average((m - m_mean) * (o - o_mean), weights=weight) / np.sqrt(
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
                    "pattern_correlation": correlation,
                    "bias_integral_PgC_yr": np.sum((m - o) * weight) * MOLAR_MASS_C / 1e15,
                }
            )
    return global_rows, regional_rows


lat_common = np.arange(-89, 90, 2.0)
lon_common = np.arange(1, 360, 2.0)
mapped = {"lat": lat_common, "lon": lon_common}
summary_rows, regional_rows, grid_tables = [], [], []
comparison_rmse = pd.read_csv(tables / "global_spatial_metrics.csv").set_index("model")

for model in MODELS:
    branch = branch_year(raw / model / "historical" / "fgco2.nc")
    with xr.open_dataset(raw / model / "historical" / "areacello.nc") as ds:
        area = ds.areacello.values.astype("float64")

    with xr.open_dataset(raw / model / "piControl" / "fgco2.nc", decode_times=time_coder) as ds:
        years = np.array([time.year for time in ds.time.values])
        indices = np.where((years >= branch) & (years <= branch + 250))[0]
        if len(indices) != 251:
            raise ValueError(f"{model}: 分支对齐piControl时段不足251年")
        seconds = record_seconds(ds)[indices]
        annual = ds.fgco2.isel(time=indices).values.astype("float64")
        annual *= seconds[:, None, None] * 1000.0 / MOLAR_MASS_C
        if "lat" in ds.coords:
            lon, lat = np.meshgrid(ds.lon.values % 360.0, ds.lat.values)
        else:
            lon, lat = ds.longitude.values % 360.0, ds.latitude.values

    valid = np.isfinite(area) & (area > 0) & np.all(np.isfinite(annual), axis=0)
    y = annual[:, valid]
    t = np.arange(251, dtype="float64")
    tc = t - t.mean()
    sxx = np.sum(tc**2)
    mean = y.mean(axis=0)
    slope = np.sum(tc[:, None] * y, axis=0) / sxx
    residual = y - mean - tc[:, None] * slope

    # Newey–West HAC标准误，Bartlett核最大滞后30年。
    score = residual * tc[:, None]
    long_run_variance = np.sum(score**2, axis=0)
    for lag in range(1, 31):
        weight = 1.0 - lag / 31.0
        long_run_variance += 2.0 * weight * np.sum(score[lag:] * score[:-lag], axis=0)
    standard_error = np.sqrt(np.maximum(long_run_variance, 0.0)) / sxx
    p_values = np.full(slope.shape, np.nan)
    good_se = standard_error > 0
    p_values[good_se] = 2.0 * norm.sf(np.abs(slope[good_se] / standard_error[good_se]))

    slope_native = np.full(area.shape, np.nan)
    p_native = np.full(area.shape, np.nan)
    slope_native[valid], p_native[valid] = slope, p_values
    raw_p05 = np.isfinite(p_native) & (p_native < 0.05)
    significant, cutoff = fdr_mask(p_native)
    drift_native = slope_native * (1998 - 1850)
    required = significant & (np.abs(drift_native) >= PRACTICAL_LIMIT)

    rows, columns = np.where(raw_p05)
    grid_tables.append(
        pd.DataFrame(
            {
                "model": model,
                "native_row": rows,
                "native_column": columns,
                "longitude": lon[raw_p05],
                "latitude": lat[raw_p05],
                "p_value": p_native[raw_p05],
                "trend_mol_m2_yr_per_decade": 10 * slope_native[raw_p05],
                "drift_1982_2014_mol_m2_yr": drift_native[raw_p05],
                "FDR_retained": significant[raw_p05],
                "correction_required": required[raw_p05],
                "area_m2": area[raw_p05],
            }
        )
    )

    significant_float = np.where(valid, significant.astype(float), np.nan)
    raw_p05_float = np.where(valid, raw_p05.astype(float), np.nan)
    required_float = np.where(valid, required.astype(float), np.nan)
    correction_native = np.where(valid, np.where(required, -drift_native, 0.0), np.nan)
    mapped[f"{model}_drift"] = conservative_bin(drift_native, lon, lat, area)
    mapped[f"{model}_significant_fraction"] = conservative_bin(
        significant_float, lon, lat, area
    )
    mapped[f"{model}_raw_p05_fraction"] = conservative_bin(raw_p05_float, lon, lat, area)
    mapped[f"{model}_required_fraction"] = conservative_bin(required_float, lon, lat, area)
    mapped[f"{model}_required_correction"] = conservative_bin(
        correction_native, lon, lat, area
    )

    ocean_area = area[valid].sum()
    drift_rms = np.sqrt(np.sum(drift_native[valid] ** 2 * area[valid]) / ocean_area)
    bias_rmse = comparison_rmse.loc[model, "RMSE_mol_m2_yr"]
    summary_rows.append(
        {
            "model": model,
            "branch_year_in_piControl": branch,
            "control_start_year": branch,
            "control_end_year": branch + 250,
            "n_years": len(indices),
            "HAC_max_lag_years": 30,
            "FDR_q": 0.05,
            "FDR_p_cutoff": cutoff,
            "practical_limit_mol_m2_yr": PRACTICAL_LIMIT,
            "area_mean_drift_mol_m2_yr": np.sum(drift_native[valid] * area[valid]) / ocean_area,
            "area_RMS_drift_mol_m2_yr": drift_rms,
            "drift_RMS_over_model_obs_RMSE": drift_rms / bias_rmse,
            "FDR_significant_area_percent": 100 * area[significant].sum() / ocean_area,
            "raw_p05_area_percent": 100 * area[raw_p05].sum() / ocean_area,
            "correction_required_area_percent": 100 * area[required].sum() / ocean_area,
            "absolute_drift_95th_percentile_mol_m2_yr": np.nanpercentile(
                np.abs(drift_native[valid]), 95
            ),
        }
    )

    regions = {
        "Southern Ocean (60S-40S)": (lat >= -60) & (lat < -40),
        "Equatorial Pacific": (lat >= -10) & (lat <= 10) & (lon >= 140) & (lon < 280),
        "North Atlantic": (lat >= 30) & (lat < 60) & ((lon >= 280) | (lon < 20)),
        "North Pacific": (lat >= 30) & (lat < 60) & (lon >= 120) & (lon < 240),
        "Indian Ocean": (lat >= -30) & (lat < 20) & (lon >= 20) & (lon < 120),
    }
    for region, region_mask in regions.items():
        mask = valid & region_mask
        region_area = area[mask].sum()
        regional_rows.append(
            {
                "model": model,
                "region": region,
                "area_mean_drift_mol_m2_yr": np.sum(drift_native[mask] * area[mask]) / region_area,
                "area_RMS_drift_mol_m2_yr": np.sqrt(
                    np.sum(drift_native[mask] ** 2 * area[mask]) / region_area
                ),
                "FDR_significant_area_percent": 100 * area[mask & significant].sum() / region_area,
                "raw_p05_area_percent": 100 * area[mask & raw_p05].sum() / region_area,
                "correction_required_area_percent": 100 * area[mask & required].sum() / region_area,
            }
        )

# 将局地订正施加到原始1982–2014空间场，并按相同方法重新计算指标。
base = np.load(intermediate / "spatial_fields_2deg.npz")
corrected = {"Landschuetzer": base["Landschuetzer"]}
for model in MODELS:
    corrected[model] = base[model] + mapped[f"{model}_required_correction"]

lon2, lat2 = np.meshgrid(lon_common, lat_common)
common_regions = {
    "Southern Ocean (60S-40S)": (lat2 >= -60) & (lat2 < -40),
    "Equatorial Pacific": (lat2 >= -10) & (lat2 <= 10) & (lon2 >= 140) & (lon2 < 280),
    "North Atlantic": (lat2 >= 30) & (lat2 < 60) & ((lon2 >= 280) | (lon2 < 20)),
    "North Pacific": (lat2 >= 30) & (lat2 < 60) & (lon2 >= 120) & (lon2 < 240),
    "Indian Ocean": (lat2 >= -30) & (lat2 < 20) & (lon2 >= 20) & (lon2 < 120),
}
corrected_global, corrected_regional = spatial_metrics(
    corrected, base["obs_area"], common_regions
)

np.savez_compressed(intermediate / "spatial_drift_fields_2deg.npz", **mapped)
np.savez_compressed(
    intermediate / "spatial_fields_2deg_drift_corrected.npz",
    lat=lat_common,
    lon=lon_common,
    obs_area=base["obs_area"],
    **corrected,
)
pd.DataFrame(summary_rows).to_csv(tables / "spatial_drift_summary.csv", index=False)
pd.DataFrame(regional_rows).to_csv(tables / "spatial_drift_regional_summary.csv", index=False)
pd.DataFrame(corrected_global).to_csv(
    tables / "global_spatial_metrics_drift_corrected.csv", index=False
)
pd.DataFrame(corrected_regional).to_csv(
    tables / "regional_spatial_metrics_drift_corrected.csv", index=False
)
pd.concat(grid_tables, ignore_index=True).to_csv(
    tables / "spatial_raw_p05_native_gridcells.csv", index=False
)
print(pd.DataFrame(summary_rows).to_string(index=False))
print(f"\n结果目录：{result_dir}")
