#!/usr/bin/env python
"""绘制piControl局地漂移及需要订正区域。"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


MODELS = ["GFDL-ESM4", "MPI-ESM1-2-LR", "CanESM5"]

parser = argparse.ArgumentParser()
parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[2])
parser.add_argument("--results-dir", type=Path, default=None)
args = parser.parse_args()

root = args.project_root.resolve()
results = (args.results_dir or root / "results" / "02_spatial_distribution").resolve()
tables, intermediate, figures = results / "tables", results / "intermediate", results / "figures"
figures.mkdir(parents=True, exist_ok=True)
data = np.load(intermediate / "spatial_drift_fields_2deg.npz")
summary = pd.read_csv(tables / "spatial_drift_summary.csv").set_index("model")
lon_edges, lat_edges = np.arange(0, 362, 2), np.arange(-90, 92, 2)

with (root / "data" / "obs" / "shapefiles" / "ne_110m_coastline.geojson").open(
    encoding="utf-8"
) as file:
    coast_json = json.load(file)
coastlines = []
for feature in coast_json["features"]:
    geometry = feature["geometry"]
    parts = [geometry["coordinates"]] if geometry["type"] == "LineString" else geometry["coordinates"]
    for part in parts:
        points = np.asarray(part)
        x, y = points[:, 0] % 360, points[:, 1]
        for segment in np.split(np.arange(len(x)), np.where(np.abs(np.diff(x)) > 180)[0] + 1):
            if len(segment) > 1:
                coastlines.append((x[segment], y[segment]))


def draw_map(ax, field, title, limit=0.3):
    mesh = ax.pcolormesh(
        lon_edges, lat_edges, field, cmap="RdBu", vmin=-limit, vmax=limit,
        shading="flat", rasterized=True
    )
    for x, y in coastlines:
        ax.plot(x, y, color="0.15", linewidth=0.35)
    ax.set(xlim=(0, 360), ylim=(-90, 90), title=title)
    ax.set_xticks([0, 60, 120, 180, 240, 300, 360])
    ax.set_xticklabels(["0°", "60°E", "120°E", "180°", "120°W", "60°W", "0°"])
    ax.set_yticks([-60, -30, 0, 30, 60])
    ax.grid(color="0.6", linewidth=0.3, alpha=0.35)
    return mesh


def draw_fraction(ax, field, title):
    mesh = ax.pcolormesh(
        lon_edges, lat_edges, 100 * field, cmap="Reds", vmin=0, vmax=100,
        shading="flat", rasterized=True
    )
    for x, y in coastlines:
        ax.plot(x, y, color="0.15", linewidth=0.35)
    ax.set(xlim=(0, 360), ylim=(-90, 90), title=title)
    ax.set_xticks([0, 60, 120, 180, 240, 300, 360])
    ax.set_xticklabels(["0°", "60°E", "120°E", "180°", "120°W", "60°W", "0°"])
    ax.set_yticks([-60, -30, 0, 30, 60])
    ax.grid(color="0.6", linewidth=0.3, alpha=0.35)
    return mesh


fig, axes = plt.subplots(2, 3, figsize=(14, 7.0), constrained_layout=True)
for column, model in enumerate(MODELS):
    row = summary.loc[model]
    drift = data[f"{model}_drift"]
    fraction = data[f"{model}_required_fraction"]
    top = draw_map(
        axes[0, column], drift,
        f"{model}\nRMS drift / spatial RMSE = {row.drift_RMS_over_model_obs_RMSE:.1%}"
    )
    correction = np.where(fraction >= 0.5, -drift, np.nan)
    bottom = draw_map(
        axes[1, column], correction,
        f"Required ocean area = {row.correction_required_area_percent:.2f}%"
    )
axes[0, 0].set_ylabel("All estimated drift\nLatitude")
axes[1, 0].set_ylabel("Correction-required regions\nLatitude")
fig.suptitle(
    "Local piControl drift affecting the 1982–2014 climatology\n"
    "Bottom: FDR q < 0.05 and |drift| ≥ 0.10 mol C m$^{-2}$ yr$^{-1}$"
)
fig.colorbar(
    top, ax=axes[0], orientation="horizontal", shrink=0.5, pad=0.04,
    fraction=0.025, aspect=50
).set_label("Estimated drift contribution (mol C m$^{-2}$ yr$^{-1}$); blue = positive")
fig.colorbar(
    bottom, ax=axes[1], orientation="horizontal", shrink=0.5, pad=0.08,
    fraction=0.025, aspect=50
).set_label("Correction to apply (mol C m$^{-2}$ yr$^{-1}$); blue = increase uptake")
fig.savefig(figures / "spatial_drift_required_regions.png", dpi=300, bbox_inches="tight")
fig.savefig(figures / "spatial_drift_required_regions.pdf", bbox_inches="tight")
plt.close(fig)

# 未订正p<0.05与FDR保留区域对比。
fig, axes = plt.subplots(2, 3, figsize=(14, 6.8), constrained_layout=True)
for column, model in enumerate(MODELS):
    row = summary.loc[model]
    mesh = draw_fraction(
        axes[0, column], data[f"{model}_raw_p05_fraction"],
        f"{model}\nraw p < 0.05 area = {row.raw_p05_area_percent:.1f}%"
    )
    mesh = draw_fraction(
        axes[1, column], data[f"{model}_significant_fraction"],
        f"FDR-retained area = {row.FDR_significant_area_percent:.2f}%"
    )
axes[0, 0].set_ylabel("Raw p < 0.05\nLatitude")
axes[1, 0].set_ylabel("After FDR q = 0.05\nLatitude")
fig.suptitle("Grid-cell trend significance before and after multiple-testing correction")
fig.colorbar(
    mesh, ax=axes, orientation="horizontal", shrink=0.5, pad=0.06,
    fraction=0.025, aspect=50
).set_label("Fraction of source-ocean area significant within each 2° cell (%)")
fig.savefig(figures / "spatial_raw_p05_vs_fdr.png", dpi=300, bbox_inches="tight")
fig.savefig(figures / "spatial_raw_p05_vs_fdr.pdf", bbox_inches="tight")
plt.close(fig)

# 根据漂移检验标准进行局地订正后的空间对比。
corrected = np.load(intermediate / "spatial_fields_2deg_drift_corrected.npz")
metrics = pd.read_csv(tables / "global_spatial_metrics_drift_corrected.csv").set_index("model")
fig, axes = plt.subplots(2, 2, figsize=(13, 6.8), constrained_layout=True)
for ax, name, key in zip(
    axes.flat,
    ["Landschützer v2016", *MODELS],
    ["Landschuetzer", *MODELS],
):
    mesh = draw_map(ax, corrected[key], name, limit=3.0)
for ax in axes[:, 0]:
    ax.set_ylabel("Latitude")
fig.suptitle(
    "Air–sea CO$_2$ flux climatology, 1982–2014\n"
    "local drift correction only where FDR q < 0.05 and |drift| ≥ 0.10"
)
fig.colorbar(
    mesh, ax=axes, orientation="horizontal", shrink=0.55, pad=0.06,
    fraction=0.025, aspect=50
).set_label("mol C m$^{-2}$ yr$^{-1}$")
fig.savefig(figures / "spatial_flux_climatology_drift_corrected.png", dpi=300, bbox_inches="tight")
fig.savefig(figures / "spatial_flux_climatology_drift_corrected.pdf", bbox_inches="tight")
plt.close(fig)

fig, axes = plt.subplots(1, 3, figsize=(14, 3.6), constrained_layout=True)
for ax, model in zip(axes, MODELS):
    row = metrics.loc[model]
    mesh = draw_map(
        ax,
        corrected[model] - corrected["Landschuetzer"],
        f"{model}\nr = {row.pattern_correlation:.2f}, RMSE = {row.RMSE_mol_m2_yr:.2f}",
        limit=2.0,
    )
axes[0].set_ylabel("Latitude")
fig.suptitle("Drift-corrected CMIP6 minus Landschützer, 1982–2014")
fig.colorbar(
    mesh, ax=axes, orientation="horizontal", shrink=0.5, pad=0.12,
    fraction=0.035, aspect=45
).set_label("Bias (mol C m$^{-2}$ yr$^{-1}$); blue = excessive uptake")
fig.savefig(figures / "spatial_bias_maps_drift_corrected.png", dpi=300, bbox_inches="tight")
fig.savefig(figures / "spatial_bias_maps_drift_corrected.pdf", bbox_inches="tight")
plt.close(fig)
print(f"图已保存到：{results}")
