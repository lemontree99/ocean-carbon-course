#!/usr/bin/env python
"""用dissic×volcello独立核算全球DIC库存，并与fgco2积分比较。

库存变化与累计海气通量来自同一试验，不扣除piControl背景收支，
与最终报告表1的口径一致。
"""

from __future__ import annotations

import argparse
import calendar
import json
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr


MODELS = ["GFDL-ESM4", "MPI-ESM1-2-LR", "CanESM5"]
PERIODS = {"historical": (1850, 2014), "ssp585": (2015, 2100)}
CARBON_G_PER_MOL = 12.0

parser = argparse.ArgumentParser()
parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[2])
args = parser.parse_args()
root = args.project_root.resolve()
raw = root / "data" / "raw"
out = root / "results" / "05_self_consistency" / "tables"
out.mkdir(parents=True, exist_ok=True)
time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)


def record_seconds(ds):
    """每个fgco2记录代表的秒数。"""
    bounds = ds.time.attrs.get("bounds")
    if bounds and bounds in ds:
        return np.array([(b - a).total_seconds() for a, b in ds[bounds].values])
    years = np.array([t.year for t in ds.time.values])
    cal = str(ds.time.encoding.get("calendar", "standard")).lower()
    if cal in {"noleap", "365_day"}:
        days = np.full(len(years), 365)
    elif cal == "360_day":
        days = np.full(len(years), 360)
    else:
        days = np.array([366 if calendar.isleap(int(y)) else 365 for y in years])
    return days * 86400.0


def inventory(model, experiment):
    """返回端点日期和Σ(dissic×volcello)库存，单位Pg C。"""
    dpath = raw / model / experiment / "dissic_slices.nc"
    vpath = raw / model / experiment / "volcello.nc"

    with xr.open_dataset(dpath, decode_times=time_coder) as dds, xr.open_dataset(vpath, decode_times=time_coder) as vds:
        dic = dds.dissic.astype("float64")
        volume = vds.volcello.astype("float64")
        values = (dic * volume).sum([d for d in dic.dims if d != "time"], skipna=True).values * CARBON_G_PER_MOL / 1e15
        dates = dds.time.values.copy()
    return dates, values


def monthly_flux(model, experiment, start, end):
    """月平均库存中心之间的海气通量：首尾月各取一半。"""
    with xr.open_dataset(raw / model / "historical" / "areacello.nc") as ads:
        area = ads.areacello.astype("float64").load()
    with xr.open_dataset(raw / model / experiment / "fgco2.nc", decode_times=time_coder) as ds:
        indices = np.where([(t >= start) and (t <= end) for t in ds.time.values])[0]
        seconds = record_seconds(ds)[indices]
        rate = np.empty(len(indices))
        for first in range(0, len(indices), 120):
            ix = indices[first:first + 120]
            rate[first:first + len(ix)] = (ds.fgco2.isel(time=ix).astype("float64") * area).sum(area.dims, skipna=True).values
    weights = np.ones(len(indices))
    weights[[0, -1]] = 0.5
    return float(np.sum(rate * seconds * weights) / 1e12)


inventory_rows, primary_rows = [], []
for model in MODELS:
    for experiment in PERIODS:
        dates, stock = inventory(model, experiment)
        raw_delta = float(stock[-1] - stock[0])
        raw_flux = monthly_flux(model, experiment, dates[0], dates[-1])
        residual = raw_delta - raw_flux
        cal = str(type(dates[0]).__name__).lower()
        year_days = 365 if "noleap" in cal else 365.2425
        elapsed = (dates[-1] - dates[0]).total_seconds() / (year_days * 86400)

        inventory_rows.extend([
            {"model": model, "experiment": experiment, "endpoint": "start", "date": str(dates[0]), "DIC_inventory_PgC": stock[0]},
            {"model": model, "experiment": experiment, "endpoint": "end", "date": str(dates[-1]), "DIC_inventory_PgC": stock[-1]},
        ])
        primary_rows.append({
            "model": model,
            "experiment": experiment,
            "elapsed_years": elapsed,
            "DIC_inventory_change_PgC": raw_delta,
            "fgco2_integral_PgC": raw_flux,
            "residual_inventory_minus_flux_PgC": residual,
            "residual_percent_of_flux": 100 * residual / abs(raw_flux),
            "closure_ratio_DIC_to_fgco2": raw_delta / raw_flux,
        })

primary = pd.DataFrame(primary_rows)
primary.to_csv(out / "raw_closure.csv", index=False)
pd.DataFrame(inventory_rows).to_csv(out / "DIC_inventory_endpoints.csv", index=False)

(out / "calculation_metadata.json").write_text(json.dumps({
    "inventory_formula": "sum(dissic [mol m-3] * official volcello [m3]) * 12 g C mol-1 / 1e15 g Pg-1",
    "flux_formula": "sum(fgco2 [kg C m-2 s-1] * areacello [m2] * seconds) / 1e12 kg Pg-1",
    "time_alignment": "monthly-mean stock centers (mid-month); first and last monthly flux each receive half weight",
    "drift_correction": "none; inventory and flux are compared consistently within each experiment",
    "residual_sign": "positive means DIC inventory gain exceeds air-sea uptake",
    "official_volume_files": "data/raw/<model>/<experiment>/volcello.nc",
    "warning": "DIC is not total ocean carbon; residual contains omitted carbon pools and external boundary terms",
}, ensure_ascii=False, indent=2), encoding="utf-8")

print(primary.to_string(index=False))
