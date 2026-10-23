#!/usr/bin/env python
"""计算年度、固定窗口和长期累计海洋碳汇效率。"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


MODELS = ["GFDL-ESM4", "MPI-ESM1-2-LR", "CanESM5"]
WINDOWS = [20, 30, 40]

parser = argparse.ArgumentParser()
parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[2])
parser.add_argument("--output-dir", type=Path, default=None)
args = parser.parse_args()

root = args.project_root.resolve()
result_dir = (args.output_dir or root / "results" / "03_sink_efficiency").resolve()
tables = result_dir / "tables"
tables.mkdir(parents=True, exist_ok=True)

sink = pd.read_csv(root / "results" / "01_global_ocean_sink" / "tables" / "global_ocean_sink_annual.csv")
co2 = pd.read_csv(root / "data" / "obs" / "co2" / "co2_concentration_1850_2100.csv")
control = pd.read_csv(
    root / "results" / "01_global_ocean_sink" / "tables" / "picontrol_correction_summary.csv"
).set_index("model")

expected_years = np.arange(1850, 2101)
assert np.array_equal(co2.year.to_numpy(), expected_years)
for model in MODELS:
    years = sink.loc[sink.model == model, "year"].to_numpy()
    assert np.array_equal(years, expected_years), f"{model}: 年份不连续"

data = sink.merge(co2[["year", "co2_ppm"]], on="year", validate="many_to_one")
data["annual_co2_change_ppm"] = data.groupby("model").co2_ppm.diff()
data["annual_efficiency_PgC_per_ppm"] = (
    data.sink_corrected_PgC_yr / data.annual_co2_change_ppm.replace(0, np.nan)
)

parts = []
for model in MODELS:
    x = data[data.model == model].sort_values("year").copy()
    x["cumulative_sink_since_1850_PgC"] = x.sink_corrected_PgC_yr.where(x.year > 1850, 0).cumsum()
    x["cumulative_co2_increase_since_1850_ppm"] = x.co2_ppm - x.co2_ppm.iloc[0]
    x["cumulative_efficiency_PgC_per_ppm"] = (
        x.cumulative_sink_since_1850_PgC
        / x.cumulative_co2_increase_since_1850_ppm.replace(0, np.nan)
    )
    for window in WINDOWS:
        x[f"efficiency_{window}yr_PgC_per_ppm"] = (
            x.sink_corrected_PgC_yr.rolling(window).sum()
            / (x.co2_ppm - x.co2_ppm.shift(window))
        )
    baseline = control.loc[model, "piControl_fit_intercept_PgC_yr"]
    x["efficiency_30yr_control_baseline_removed_PgC_per_ppm"] = (
        (x.sink_corrected_PgC_yr - baseline).rolling(30).sum()
        / (x.co2_ppm - x.co2_ppm.shift(30))
    )
    parts.append(x)
annual = pd.concat(parts, ignore_index=True)


def interval_efficiency(x: pd.DataFrame, first_flux_year: int, last_flux_year: int, column: str):
    """年度通量求和；CO2基准为首个通量年份的前一年。"""
    start = first_flux_year - 1
    numerator = x.loc[x.year.between(first_flux_year, last_flux_year), column].sum()
    denominator = x.loc[x.year == last_flux_year, "co2_ppm"].iloc[0] - x.loc[x.year == start, "co2_ppm"].iloc[0]
    return numerator / denominator, numerator, denominator


periods = {
    "historical_1851_2014": (1851, 2014),
    "ssp585_2015_2100": (2015, 2100),
    "recent_historical_1985_2014": (1985, 2014),
    "late_century_2071_2100": (2071, 2100),
}
summary_rows, sensitivity_rows = [], []
for model in MODELS:
    x = annual[annual.model == model].copy()
    baseline = control.loc[model, "piControl_fit_intercept_PgC_yr"]
    x["sink_control_baseline_removed_PgC_yr"] = x.sink_corrected_PgC_yr - baseline

    for period, (first, last) in periods.items():
        efficiency, uptake, increase = interval_efficiency(x, first, last, "sink_corrected_PgC_yr")
        summary_rows.append(
            {
                "model": model,
                "period": period,
                "first_flux_year": first,
                "last_flux_year": last,
                "cumulative_ocean_uptake_PgC": uptake,
                "co2_increase_ppm": increase,
                "efficiency_PgC_per_ppm": efficiency,
            }
        )

    series = {
        "raw": "sink_raw_PgC_yr",
        "drift_corrected": "sink_corrected_PgC_yr",
        "drift_and_control_baseline_corrected": "sink_control_baseline_removed_PgC_yr",
    }
    for label, column in series.items():
        for window in WINDOWS:
            early, _, _ = interval_efficiency(x, 2015 - window, 2014, column)
            late, _, _ = interval_efficiency(x, 2101 - window, 2100, column)
            sensitivity_rows.append(
                {
                    "model": model,
                    "series": label,
                    "window_years": window,
                    "efficiency_ending_2014_PgC_per_ppm": early,
                    "efficiency_ending_2100_PgC_per_ppm": late,
                    "absolute_change_PgC_per_ppm": late - early,
                    "percent_change": 100 * (late / early - 1),
                }
            )

summary = pd.DataFrame(summary_rows)
sensitivity = pd.DataFrame(sensitivity_rows)
annual.to_csv(tables / "sink_efficiency_annual.csv", index=False)
summary.to_csv(tables / "sink_efficiency_summary.csv", index=False)
sensitivity.to_csv(tables / "efficiency_window_sensitivity.csv", index=False)

nonpositive = co2.loc[co2.co2_ppm.diff() <= 0, "year"].tolist()
with (tables / "calculation_metadata.json").open("w", encoding="utf-8") as file:
    json.dump(
        {
            "formula": "sum annual ocean sink / (CO2_end - CO2_start)",
            "units": "Pg C per ppm",
            "sink": "analysis 1 branch-aligned piControl drift-corrected annual fgco2 integral",
            "co2": (
                "official input4MIPs global annual mean: UoM-CMIP-1-2-0 historical "
                "joined to UoM-REMIND-MAGPIE-ssp585-1-2-1"
            ),
            "co2_data_dois": [
                "10.22033/ESGF/input4MIPs.1118",
                "10.22033/ESGF/input4MIPs.9868",
            ],
            "primary_window_years": 30,
            "sensitivity_windows_years": WINDOWS,
            "annual_efficiency": "sink(year) / [CO2(year) - CO2(year-1)]; no smoothing",
            "annual_nonpositive_co2_increment_years": nonpositive,
            "plot_start_year": 1850,
            "plot_y_scale": (
                "symmetric logarithmic because official input4MIPs annual CO2 growth is "
                "very small in 1944 and 1946 and slightly negative in 1945"
            ),
        },
        file,
        ensure_ascii=False,
        indent=2,
    )

print(
    sensitivity[(sensitivity.series == "drift_corrected") & (sensitivity.window_years == 30)]
    .drop(columns="series")
    .to_string(index=False)
)
print(f"\n结果目录：{result_dir}")
