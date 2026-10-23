#!/usr/bin/env python
"""绘制DIC库存变化与fgco2累计通量的未订正自洽性比较。"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


parser = argparse.ArgumentParser()
parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[2])
args = parser.parse_args()
root = args.project_root.resolve()
table = pd.read_csv(root / "results" / "05_self_consistency" / "tables" / "raw_closure.csv")
out = root / "results" / "05_self_consistency" / "figures"
out.mkdir(parents=True, exist_ok=True)

labels = [f"{m}\n{'Historical' if e == 'historical' else 'SSP5-8.5'}" for m, e in zip(table.model, table.experiment)]
x = np.arange(len(table)); width = 0.34
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
fig, ax = plt.subplots(figsize=(10.2, 5.5))
fig.subplots_adjust(left=0.09, right=0.98, top=0.82, bottom=0.23)
b1 = ax.bar(x - width / 2, table.fgco2_integral_PgC, width, color="#263238", label=r"Integrated $fgco_2$")
b2 = ax.bar(x + width / 2, table.DIC_inventory_change_PgC, width, color="#E8892E", label="DIC inventory change")
ax.bar_label(b1, fmt="%.1f", padding=3, fontsize=8.5)
ax.bar_label(b2, fmt="%.1f", padding=3, fontsize=8.5)
for i, row in table.iterrows():
    top = max(row.fgco2_integral_PgC, row.DIC_inventory_change_PgC)
    ax.text(i, top + 44, f"residual {row.residual_percent_of_flux:+.1f}%", ha="center", fontsize=8.5, color="#53666F", bbox=dict(facecolor="white", edgecolor="none", pad=1.5))
fig.text(0.09, 0.95, "Uncorrected ocean-carbon budget closure", ha="left", va="top", fontsize=15, fontweight="bold")
fig.text(0.09, 0.895, r"No piControl correction; first and last monthly flux weighted by 0.5", ha="left", va="top", color="#53666F")
ax.set_ylabel("Cumulative carbon (Pg C)")
ax.set_xticks(x, labels)
ax.set_ylim(0, max(table.DIC_inventory_change_PgC.max(), table.fgco2_integral_PgC.max()) * 1.22)
ax.grid(axis="y", color="#DCE3E6", linewidth=0.8); ax.set_axisbelow(True)
ax.axvline(1.5, color="#BBC6CA", lw=1); ax.axvline(3.5, color="#BBC6CA", lw=1)
ax.legend(frameon=False, ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.15))
for suffix in ["png", "pdf"]:
    fig.savefig(out / f"dic_fgco2_closure.{suffix}", dpi=400, bbox_inches="tight")
plt.close(fig)