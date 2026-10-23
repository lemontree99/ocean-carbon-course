#!/usr/bin/env python
"""绘制未经平滑的年度效率及30年累计窗口效率。"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D


MODELS = ["GFDL-ESM4", "MPI-ESM1-2-LR", "CanESM5"]
COLORS = {"GFDL-ESM4": "#0072B2", "MPI-ESM1-2-LR": "#D55E00", "CanESM5": "#009E73"}

parser = argparse.ArgumentParser()
parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[2])
parser.add_argument("--results-dir", type=Path, default=None)
args = parser.parse_args()

root = args.project_root.resolve()
results = (args.results_dir or root / "results" / "03_sink_efficiency").resolve()
tables, figures = results / "tables", results / "figures"
figures.mkdir(parents=True, exist_ok=True)

annual = pd.read_csv(tables / "sink_efficiency_annual.csv")
sensitivity = pd.read_csv(tables / "efficiency_window_sensitivity.csv")
change = sensitivity[(sensitivity.series == "drift_corrected") & (sensitivity.window_years == 30)].set_index("model")

fig, ax = plt.subplots(figsize=(11.5, 5.8), constrained_layout=True)
for model in MODELS:
    x = annual[annual.model == model]
    row = change.loc[model]
    label = (
        f"{model}: {row.efficiency_ending_2014_PgC_per_ppm:.2f} → "
        f"{row.efficiency_ending_2100_PgC_per_ppm:.2f} "
        f"({row.percent_change:.0f}%)"
    )
    ax.plot(
        x.year,
        x.efficiency_30yr_PgC_per_ppm,
        color=COLORS[model],
        linewidth=4.0,
        alpha=0.28,
        zorder=1,
    )
    ax.plot(
        x.year,
        x.annual_efficiency_PgC_per_ppm,
        color=COLORS[model],
        linewidth=0.9,
        label=label,
        zorder=2,
    )

ax.axvline(2015, color="0.35", linestyle="--", linewidth=1.0)
ax.axhline(0, color="0.65", linewidth=0.8)
ax.text(2017, 0.94, "SSP5-8.5", color="0.3", fontsize=9, transform=ax.get_xaxis_transform())
ax.set_yscale("symlog", linthresh=0.5, linscale=1.0, base=10)
ax.set_yticks([-100, -30, -10, -3, -1, -0.5, 0, 0.5, 1, 3, 10, 30, 100])
ax.set_yticklabels(["−100", "−30", "−10", "−3", "−1", "−0.5", "0", "0.5", "1", "3", "10", "30", "100"])
ax.set(
    xlim=(1850, 2100),
    ylim=(-180, 180),
    xlabel="Year",
    ylabel="Ocean sink efficiency (Pg C ppm$^{-1}$; symmetric-log scale)",
    title="Ocean carbon-sink efficiency: annual values and 30-year cumulative windows",
)
ax.grid(axis="y", color="0.90", linewidth=0.7)
handles, labels = ax.get_legend_handles_labels()
handles.append(Line2D([0], [0], color="0.35", linewidth=4.0, alpha=0.28))
labels.append("Pale broad lines: 30-year cumulative windows")
ax.legend(
    handles,
    labels,
    loc="upper right",
    frameon=False,
    fontsize=9,
    title="30-year efficiency ending in 2014 → 2100",
    title_fontsize=9,
)

fig.savefig(figures / "sink_efficiency_timeseries.png", dpi=300, bbox_inches="tight")
fig.savefig(figures / "sink_efficiency_timeseries.pdf", bbox_inches="tight")
plt.close(fig)
print(f"图已保存到：{figures}")
