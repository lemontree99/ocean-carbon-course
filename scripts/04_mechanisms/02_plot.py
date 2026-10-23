#!/usr/bin/env python
"""绘制机制量级和两个30年窗口的海洋状态对比。"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


root = Path(__file__).resolve().parents[2]
folder = root / "results/04_mechanisms"
table = pd.read_csv(folder / "tables/mechanism_magnitude_summary.csv")
state = pd.read_csv(folder / "tables/state_variable_comparison.csv")
hist = state[state.period == "1985-2014"].set_index("model").loc[table.model]
late = state[state.period == "2071-2100"].set_index("model").loc[table.model]

plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
x = np.arange(len(table))
fig = plt.figure(figsize=(11, 7.2), constrained_layout=True)
gs = fig.add_gridspec(2, 3, height_ratios=[1.35, 1])

ax = fig.add_subplot(gs[0, :])
w = 0.23
bars = [
    ax.bar(x - w, table.actual_total_loss_PgC_per_ppm, w, color="#263238", label="Actual efficiency loss"),
    ax.bar(x, table.buffer_estimated_loss_PgC_per_ppm, w, color="#168C8C", label="Buffer weakening estimate"),
    ax.bar(x + w, table.warming_joint_estimated_loss_PgC_per_ppm, w, color="#E8892E", label="Warming-related physical estimate"),
]
for group in bars:
    ax.bar_label(group, fmt="%.2f", padding=3)
ax.set_xticks(x, table.model)
ax.set_ylabel(r"Efficiency loss (Pg C ppm$^{-1}$)")
ax.set_ylim(0, table.actual_total_loss_PgC_per_ppm.max() * 1.28)
ax.grid(axis="y", color="#DCE3E6", linewidth=0.8)
ax.set_axisbelow(True)
ax.legend(frameon=False, ncol=3, loc="upper right")

items = [
    ("surface_temperature_degC", "Surface temperature (°C)", ".2f"),
    ("surface_salinity", "Surface salinity", ".2f"),
    ("surface_TA_mol_m3", r"Surface TA (mol m$^{-3}$)", ".3f"),
]
for i, (column, ylabel, fmt) in enumerate(items):
    ax = fig.add_subplot(gs[1, i])
    h, f = hist[column].to_numpy(), late[column].to_numpy()
    ax.vlines(x, h, f, color="#AAB7BF", linewidth=2)
    ax.scatter(x, h, s=42, color="#3569A8", label="1985–2014", zorder=3)
    ax.scatter(x, f, s=42, color="#D97722", label="2071–2100", zorder=3)
    for j in x:
        h_offset, f_offset = ((0, -13), (0, 7)) if h[j] <= f[j] else ((0, 7), (0, -13))
        ax.annotate(format(h[j], fmt), (j, h[j]), xytext=h_offset, textcoords="offset points", ha="center", fontsize=7)
        ax.annotate(format(f[j], fmt), (j, f[j]), xytext=f_offset, textcoords="offset points", ha="center", fontsize=7)
    ax.set_xticks(x, table.model, rotation=12)
    ax.set_xlim(-0.32, len(x) - 0.68)
    ax.set_ylabel(ylabel)
    ax.grid(axis="y", color="#E3E8EA", linewidth=0.7)
    ax.set_axisbelow(True)
    ax.margins(y=0.14)
    if i == 0:
        ax.legend(frameon=False, loc="best")

out = folder / "figures"
out.mkdir(parents=True, exist_ok=True)
for suffix in ["png", "pdf"]:
    fig.savefig(out / f"mechanism_magnitude_comparison.{suffix}", dpi=400, bbox_inches="tight")
plt.close(fig)
