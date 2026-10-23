#!/usr/bin/env python
"""计算全球海洋碳汇，并用分支对齐的piControl趋势订正时间漂移。"""

from __future__ import annotations

import argparse
import calendar
import json
from pathlib import Path

import cftime
import numpy as np
import pandas as pd
import xarray as xr


MODELS = ["GFDL-ESM4", "MPI-ESM1-2-LR", "CanESM5"]
KG_PER_PG = 1e12

parser = argparse.ArgumentParser()
parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[2])
parser.add_argument("--output-dir", type=Path, default=None)
args = parser.parse_args()

root = args.project_root.resolve()
raw = root / "data" / "raw"
result_dir = (args.output_dir or root / "results" / "01_global_ocean_sink").resolve()
out = result_dir / "tables"
out.mkdir(parents=True, exist_ok=True)
time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)


def record_seconds(ds: xr.Dataset) -> np.ndarray:
    """读取CF时间边界；若年频文件没有边界，则按日历推算一年秒数。"""
    bounds_name = ds.time.attrs.get("bounds")
    if bounds_name in ds:
        result = []
        for start, end in ds[bounds_name].values:
            delta = end - start
            result.append(
                float(delta / np.timedelta64(1, "s"))
                if isinstance(delta, np.timedelta64)
                else float(delta.total_seconds())
            )
        return np.asarray(result)

    years = np.array([time.year for time in ds.time.values])
    cal = str(ds.time.encoding.get("calendar", "standard")).lower()
    if cal in {"noleap", "365_day"}:
        days = np.full(len(years), 365)
    elif cal == "360_day":
        days = np.full(len(years), 360)
    else:
        days = np.array([366 if calendar.isleap(int(year)) else 365 for year in years])
    return days * 86400.0


def annual_sink(flux_file: Path, area: xr.DataArray) -> pd.DataFrame:
    """在原生网格积分fgco2，并按年份汇总为Pg C/yr。"""
    with xr.open_dataset(flux_file, decode_times=time_coder) as ds:
        flux = ds.fgco2
        seconds = record_seconds(ds)
        years = np.array([time.year for time in ds.time.values])

        # 分块读取，避免一次载入CanESM5的完整文件。
        rate = np.empty(len(years))
        for start in range(0, len(years), 120):
            stop = min(start + 120, len(years))
            block = flux.isel(time=slice(start, stop)).astype("float64")
            rate[start:stop] = (block * area).sum(area.dims, skipna=True).values

    mass = rate * seconds / KG_PER_PG
    return (
        pd.DataFrame({"year": years, "sink_raw_PgC_yr": mass})
        .groupby("year", as_index=False)["sink_raw_PgC_yr"]
        .sum()
    )


def branch_year(flux_file: Path) -> int:
    """把historical的分支时刻转换为piControl中的日历年。"""
    with xr.open_dataset(flux_file, decode_times=False) as ds:
        date = cftime.num2date(
            ds.attrs["branch_time_in_parent"],
            ds.attrs["parent_time_units"],
            calendar=ds.time.attrs.get("calendar", "standard"),
        )
    return int(date.year)


def block_slope_test(values, block=30, n_boot=5000, seed=0):
    """用循环移动块bootstrap检验线性趋势，同时保留年际自相关。"""
    y = np.asarray(values, dtype="float64")
    t = np.arange(len(y), dtype="float64")
    tc = t - t.mean()
    denominator = np.sum(tc**2)
    slope = np.sum(tc * y) / denominator
    residual = y - (y.mean() + slope * tc)
    rng = np.random.default_rng(seed)
    slopes = np.empty(n_boot)
    n_blocks = int(np.ceil(len(y) / block))
    offsets = np.arange(block)
    for i in range(n_boot):
        starts = rng.integers(0, len(y), n_blocks)
        sampled = residual[((starts[:, None] + offsets) % len(y)).ravel()[: len(y)]]
        slopes[i] = slope + np.sum(tc * sampled) / denominator
    lower, upper = np.percentile(slopes, [2.5, 97.5])
    null_slopes = slopes - slope
    p_value = (np.count_nonzero(np.abs(null_slopes) >= abs(slope)) + 1) / (n_boot + 1)
    return slope, lower, upper, p_value


# 1. 面积积分historical、ssp585和piControl。
experiment_parts, control_parts, area_rows, branch_years = [], [], [], {}

for model in MODELS:
    with xr.open_dataset(raw / model / "historical" / "areacello.nc") as ds:
        area = ds.areacello.astype("float64").load()

    ocean_area = float(area.where(area > 0).sum(skipna=True))
    assert 3.0e14 < ocean_area < 4.0e14
    area_rows.append({"model": model, "ocean_area_m2": ocean_area})
    branch_years[model] = branch_year(raw / model / "historical" / "fgco2.nc")

    for experiment in ["historical", "ssp585", "piControl"]:
        table = annual_sink(raw / model / experiment / "fgco2.nc", area)
        table.insert(0, "experiment", experiment)
        table.insert(0, "model", model)
        if experiment == "ssp585":
            table = table[table.year <= 2100].copy()
        (control_parts if experiment == "piControl" else experiment_parts).append(table)

annual = pd.concat(experiment_parts, ignore_index=True)
control = pd.concat(control_parts, ignore_index=True)

# 2. 在与1850–2100平行的piControl段拟合趋势，只减去相对分支时刻的漂移。
control_rows, stability_rows, subperiod_rows = [], [], []
annual["piControl_drift_PgC_yr"] = np.nan
for model_index, model in enumerate(MODELS):
    full = control[control.model == model]
    branch = branch_years[model]
    fit = full[full.year.between(branch, branch + 250)].copy()
    if len(fit) != 251:
        raise ValueError(f"{model}: 分支对齐的piControl时段不足251年")

    slope, intercept = np.polyfit(
        fit.year - branch, fit.sink_raw_PgC_yr, 1
    )
    slope, ci_low, ci_high, p_value = block_slope_test(
        fit.sink_raw_PgC_yr, seed=20261022 + model_index
    )
    mask = annual.model == model
    annual.loc[mask, "piControl_drift_PgC_yr"] = slope * (annual.loc[mask, "year"] - 1850)
    control_rows.append(
        {
            "model": model,
            "n_years_total": len(full),
            "branch_year_in_piControl": branch,
            "fit_start_year_in_piControl": branch,
            "fit_end_year_in_piControl": branch + 250,
            "n_fit_years": len(fit),
            "piControl_full_mean_PgC_yr": full.sink_raw_PgC_yr.mean(),
            "piControl_fit_std_PgC_yr": fit.sink_raw_PgC_yr.std(),
            "piControl_fit_intercept_PgC_yr": intercept,
            "piControl_trend_PgC_yr_per_decade": slope * 10,
            "drift_at_2014_PgC_yr": slope * (2014 - 1850),
            "drift_at_2100_PgC_yr": slope * (2100 - 1850),
        }
    )

    thirds = np.array_split(np.arange(len(fit)), 3)
    third_slopes = []
    for number, indices in enumerate(thirds, 1):
        part = fit.iloc[indices]
        part_slope = np.polyfit(
            part.year - branch, part.sink_raw_PgC_yr, 1
        )[0]
        third_slopes.append(part_slope)
        subperiod_rows.append(
            {
                "model": model,
                "subperiod": f"third_{number}",
                "experiment_equivalent_start_year": 1850 + int(indices[0]),
                "experiment_equivalent_end_year": 1850 + int(indices[-1]),
                "control_start_year": int(part.year.iloc[0]),
                "control_end_year": int(part.year.iloc[-1]),
                "n_years": len(part),
                "trend_PgC_yr_per_decade": part_slope * 10,
            }
        )

    standard_deviation = fit.sink_raw_PgC_yr.std()
    cumulative_drift = slope * 250
    stability_rows.append(
        {
            "model": model,
            "n_years": len(fit),
            "bootstrap_block_years": 30,
            "bootstrap_replicates": 5000,
            "trend_PgC_yr_per_decade": slope * 10,
            "CI95_low_PgC_yr_per_decade": ci_low * 10,
            "CI95_high_PgC_yr_per_decade": ci_high * 10,
            "p_value_two_sided": p_value,
            "trend_significant_95pct": p_value < 0.05,
            "first_50yr_mean_PgC_yr": fit.sink_raw_PgC_yr.iloc[:50].mean(),
            "last_50yr_mean_PgC_yr": fit.sink_raw_PgC_yr.iloc[-50:].mean(),
            "last_minus_first_50yr_mean_PgC_yr": (
                fit.sink_raw_PgC_yr.iloc[-50:].mean()
                - fit.sink_raw_PgC_yr.iloc[:50].mean()
            ),
            "control_std_PgC_yr": standard_deviation,
            "cumulative_linear_drift_1850_2100_PgC_yr": cumulative_drift,
            "absolute_drift_over_control_std": abs(cumulative_drift) / standard_deviation,
            "three_subperiod_trends_same_sign": (
                np.all(np.asarray(third_slopes) > 0)
                or np.all(np.asarray(third_slopes) < 0)
            ),
            "stability_conclusion": (
                "no significant linear drift"
                if p_value >= 0.05
                else "significant linear drift"
            ),
        }
    )

control_summary = pd.DataFrame(control_rows)
stability = pd.DataFrame(stability_rows)
subperiod_trends = pd.DataFrame(subperiod_rows)
annual["sink_corrected_PgC_yr"] = (
    annual.sink_raw_PgC_yr - annual.piControl_drift_PgC_yr
)
annual = annual.sort_values(["model", "year"]).reset_index(drop=True)

# 3. 关键时期统计。
period_rows = []
periods = {
    "historical_full": (1850, 2014),
    "historical_recent": (2005, 2014),
    "ssp585_near_term": (2015, 2024),
    "ssp585_full": (2015, 2100),
    "ssp585_late_century": (2091, 2100),
}
for model in MODELS:
    for label, (start, end) in periods.items():
        x = annual[(annual.model == model) & annual.year.between(start, end)]
        period_rows.append(
            {
                "model": model,
                "period": label,
                "start_year": start,
                "end_year": end,
                "n_years": len(x),
                "raw_mean_PgC_yr": x.sink_raw_PgC_yr.mean(),
                "corrected_mean_PgC_yr": x.sink_corrected_PgC_yr.mean(),
                "corrected_std_PgC_yr": x.sink_corrected_PgC_yr.std(),
                "corrected_trend_PgC_yr_per_decade": np.polyfit(
                    x.year, x.sink_corrected_PgC_yr, 1
                )[0]
                * 10,
            }
        )
period_summary = pd.DataFrame(period_rows)

# 4. 与GCB逐年配对。
gcb = pd.read_csv(root / "data" / "obs" / "gcb" / "gcb2025_ocean_sink.csv")
wide = annual.pivot(index="year", columns="model", values="sink_corrected_PgC_yr")
wide["Three-model mean"] = wide[MODELS].mean(axis=1)
comparison_rows = []
comparison_periods = {
    "1959-2014": (1959, 2014),
    "1982-2014": (1982, 2014),
    "2005-2014": (2005, 2014),
    "2015-2024": (2015, 2024),
    "1959-2024": (1959, 2024),
}
for label, (start, end) in comparison_periods.items():
    for model in MODELS + ["Three-model mean"]:
        x = wide[[model]].rename(columns={model: "model_sink"}).reset_index()
        x = x.merge(gcb[["year", "ocean_sink_GCB_PgC_yr"]], on="year")
        x = x[x.year.between(start, end)].dropna()
        error = x.model_sink - x.ocean_sink_GCB_PgC_yr
        comparison_rows.append(
            {
                "model": model,
                "period": label,
                "n_years": len(x),
                "model_mean_PgC_yr": x.model_sink.mean(),
                "GCB_mean_PgC_yr": x.ocean_sink_GCB_PgC_yr.mean(),
                "bias_PgC_yr": error.mean(),
                "bias_percent": 100 * error.mean() / x.ocean_sink_GCB_PgC_yr.mean(),
                "RMSE_PgC_yr": np.sqrt(np.mean(error**2)),
                "correlation": x.model_sink.corr(x.ocean_sink_GCB_PgC_yr),
                "model_trend_PgC_yr_per_decade": np.polyfit(x.year, x.model_sink, 1)[0] * 10,
                "GCB_trend_PgC_yr_per_decade": np.polyfit(
                    x.year, x.ocean_sink_GCB_PgC_yr, 1
                )[0]
                * 10,
            }
        )
comparison = pd.DataFrame(comparison_rows)

# 5. 保存结果。
annual.to_csv(out / "global_ocean_sink_annual.csv", index=False)
control.to_csv(out / "picontrol_global_sink_annual.csv", index=False)
control_summary.to_csv(out / "picontrol_correction_summary.csv", index=False)
stability.to_csv(out / "picontrol_stability_test.csv", index=False)
subperiod_trends.to_csv(out / "picontrol_subperiod_trends.csv", index=False)
period_summary.to_csv(out / "period_summary.csv", index=False)
comparison.to_csv(out / "gcb_comparison_stats.csv", index=False)
with (out / "calculation_metadata.json").open("w", encoding="utf-8") as file:
    json.dump(
        {
            "formula": "sum(fgco2 * areacello * record_seconds) / 1e12",
            "units": "Pg C yr-1",
            "positive": "downward into ocean",
            "correction": (
                "branch-aligned linear piControl drift fitted over the 251-year segment "
                "parallel to 1850-2100; correction is zero at 1850"
            ),
            "branch_year_in_piControl": branch_years,
            "stability_test": (
                "30-year circular moving-block bootstrap with 5000 replicates; "
                "three non-overlapping subperiod trend check"
            ),
            "ocean_area": area_rows,
        },
        file,
        ensure_ascii=False,
        indent=2,
    )

print("分支对齐的piControl时间漂移")
print(
    control_summary[
        [
            "model",
            "branch_year_in_piControl",
            "piControl_trend_PgC_yr_per_decade",
            "drift_at_2014_PgC_yr",
            "drift_at_2100_PgC_yr",
        ]
    ].to_string(index=False)
)
print(f"\n结果目录：{result_dir}")
