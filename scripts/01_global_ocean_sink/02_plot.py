#!/usr/bin/env python
"""绘制全球海洋碳汇时间序列；--series raw 使用未订正值。"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


MODELS = ["GFDL-ESM4", "MPI-ESM1-2-LR", "CanESM5"]
COLORS = {"GFDL-ESM4": "#0072B2", "MPI-ESM1-2-LR": "#D55E00", "CanESM5": "#009E73"}

parser = argparse.ArgumentParser()
parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[2])
parser.add_argument("--results-dir", type=Path, default=None)
parser.add_argument("--series", choices=["corrected", "raw"], default="corrected")
args = parser.parse_args()

root = args.project_root.resolve()
results = (args.results_dir or root / "results" / "01_global_ocean_sink").resolve()
tables, figures = results / "tables", results / "figures"
figures.mkdir(parents=True, exist_ok=True)
annual = pd.read_csv(tables / "global_ocean_sink_annual.csv")
gcb = pd.read_csv(root / "data" / "obs" / "gcb" / "gcb2025_ocean_sink.csv")

column = "sink_corrected_PgC_yr" if args.series == "corrected" else "sink_raw_PgC_yr"
label = "drift-corrected" if args.series == "corrected" else "uncorrected"
suffix = "" if args.series == "corrected" else "_raw"

fig, ax = plt.subplots(figsize=(11.5, 5.8), constrained_layout=True)

means = {}
for model in MODELS:
    x = annual[annual.model == model].sort_values("year")
    means[model] = x[x.year.between(1959, 2014)][column].mean()
    ax.plot(
        x.year,
        x[column],
        color=COLORS[model],
        linewidth=1.0,
        label=f"{model} ({means[model]:.2f})",
    )

gcb = gcb.dropna(subset=["ocean_sink_GCB_PgC_yr"])
means["GCB 2025"] = gcb[gcb.year.between(1959, 2014)].ocean_sink_GCB_PgC_yr.mean()
ax.fill_between(
    gcb.year,
    gcb.ocean_sink_GCB_PgC_yr - 0.4,
    gcb.ocean_sink_GCB_PgC_yr + 0.4,
    color="0.45",
    alpha=0.16,
    linewidth=0,
    label="GCB uncertainty (±0.4)",
)
ax.plot(
    gcb.year,
    gcb.ocean_sink_GCB_PgC_yr,
    color="black",
    linewidth=1.7,
    marker="o",
    markersize=2.2,
    label=f"GCB 2025 ({means['GCB 2025']:.2f})",
    zorder=5,
)

ax.axvline(2015, color="0.35", linestyle="--", linewidth=1.0)
ax.axhline(0, color="0.55", linewidth=0.8)
ax.text(2017, 0.25, "SSP5-8.5", color="0.3", fontsize=9)
ax.set(xlim=(1850, 2100), xlabel="Year", ylabel="Ocean CO$_2$ sink (Pg C yr$^{-1}$)")
ax.set_title(f"Global ocean carbon sink: {label} annual CMIP6 and GCB 2025")
ax.grid(axis="y", color="0.90", linewidth=0.7)
ax.legend(
    loc="upper left",
    ncol=2,
    frameon=False,
    fontsize=9,
    title="1959–2014 mean in parentheses (Pg C yr$^{-1}$)",
    title_fontsize=9,
)

fig.savefig(figures / f"global_ocean_sink_timeseries{suffix}.png", dpi=300, bbox_inches="tight")
fig.savefig(figures / f"global_ocean_sink_timeseries{suffix}.pdf", bbox_inches="tight")
plt.close(fig)

# 分支对齐piControl段稳定性：原始年度值和线性趋势，不做时间平滑。
control = pd.read_csv(tables / "picontrol_global_sink_annual.csv")
summary = pd.read_csv(tables / "picontrol_correction_summary.csv").set_index("model")
stability = pd.read_csv(tables / "picontrol_stability_test.csv").set_index("model")
fig, axes = plt.subplots(1, 3, figsize=(13, 3.8), sharex=True, sharey=True, constrained_layout=True)
for ax, model in zip(axes, MODELS):
    row, test = summary.loc[model], stability.loc[model]
    branch = int(row.branch_year_in_piControl)
    x = control[(control.model == model) & control.year.between(branch, branch + 250)].copy()
    elapsed = x.year - branch
    slope = test.trend_PgC_yr_per_decade / 10
    fitted = row.piControl_fit_intercept_PgC_yr + slope * elapsed
    ax.plot(elapsed, x.sink_raw_PgC_yr, color="0.45", linewidth=0.65)
    ax.plot(elapsed, fitted, color=COLORS[model], linewidth=2.0)
    ax.axhline(0, color="0.75", linewidth=0.7)
    ax.set_title(
        f"{model}\ntrend = {test.trend_PgC_yr_per_decade:.4f}, "
        f"p = {test.p_value_two_sided:.3f}"
    )
    ax.set_xlabel("Years since historical branch")
    ax.grid(axis="y", color="0.90", linewidth=0.6)
axes[0].set_ylabel("piControl ocean CO$_2$ flux (Pg C yr$^{-1}$)")
fig.suptitle("Stability of branch-aligned piControl segments (annual values, no smoothing)")
fig.savefig(figures / "picontrol_stability_check.png", dpi=300, bbox_inches="tight")
fig.savefig(figures / "picontrol_stability_check.pdf", bbox_inches="tight")
plt.close(fig)
print(f"图已保存到：{results}")
