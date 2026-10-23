#!/usr/bin/env python
"""提取机理分析所需的两个30年窗口平均场。"""

import argparse
import json
from pathlib import Path

import fsspec
import numpy as np
import pandas as pd
import xarray as xr


MODELS = ["GFDL-ESM4", "MPI-ESM1-2-LR", "CanESM5"]
WINDOWS = {"historical": (1985, 2014), "ssp585": (2071, 2100)}
DEPTHS = [0.0, 200.0, 1000.0]

parser = argparse.ArgumentParser()
parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[2])
parser.add_argument("--force", action="store_true")
parser.add_argument("--models", nargs="+", choices=MODELS, default=MODELS)
parser.add_argument("--experiments", nargs="+", choices=list(WINDOWS), default=list(WINDOWS))
args = parser.parse_args()
root = args.project_root.resolve()
catalog = pd.read_csv(root / "data/preprocessing/pangeo-cmip6.csv")
manifest = []


def source(model, experiment, variable):
    q = catalog[
        (catalog.source_id == model)
        & (catalog.experiment_id == experiment)
        & (catalog.member_id == "r1i1p1f1")
        & (catalog.table_id == "Omon")
        & (catalog.variable_id == variable)
    ]
    q = q[q.grid_label == ("gr" if model == "GFDL-ESM4" else "gn")]
    if len(q) != 1:
        raise ValueError(f"{model} {experiment} {variable}: 找到{len(q)}条记录")
    return q.iloc[0].zstore


def window_mean(remote, variable, years, depths=None):
    ds = xr.open_zarr(fsspec.get_mapper(remote, token="anon"), consolidated=True, chunks=None)
    da = ds[variable]
    if depths is not None:
        idx = [int(np.abs(ds.lev.values - z).argmin()) for z in depths]
        da = da.isel(lev=idx)

    total, ndays = None, 0.0
    for year in range(years[0], years[1] + 1):
        one = da.where(da.time.dt.year == year, drop=True)
        if one.sizes["time"] != 12:
            raise ValueError(f"{variable} {year}: 月数={one.sizes['time']}")
        days = one.time.dt.days_in_month
        part = (one * days).sum("time").load()
        total = part if total is None else total + part
        ndays += float(days.sum())
        if year == years[0] or year % 5 == 0 or year == years[1]:
            print(f"    {variable} {year}: 完成", flush=True)

    mean = total / ndays
    mean.attrs.update(da.attrs)
    mean.attrs["processing"] = f"calendar-day-weighted mean, {years[0]}-{years[1]}"
    ds.close()
    return mean


for model in args.models:
    for experiment in args.experiments:
        years = WINDOWS[experiment]
        folder = root / "data/raw" / model / experiment
        carbonate_file = folder / "mechanism_carbonate_surface_30yr_mean.nc"
        strat_file = folder / "mechanism_stratification_30yr_mean.nc"
        carbonate_sources = {v: source(model, experiment, v) for v in ["dissic", "talk"]}
        strat_sources = {v: source(model, experiment, v) for v in ["thetao", "so"]}
        print(f"\n{model} {experiment} {years[0]}-{years[1]}", flush=True)

        if args.force or not carbonate_file.exists():
            data = {v: window_mean(u, v, years, [0.0]) for v, u in carbonate_sources.items()}
            ds = xr.Dataset(data, attrs={"model": model, "experiment": experiment,
                "member_id": "r1i1p1f1", "window": f"{years[0]}-{years[1]}"})
            tmp = carbonate_file.with_suffix(".tmp.nc")
            ds.to_netcdf(tmp, encoding={v: {"zlib": True, "complevel": 4} for v in data})
            tmp.replace(carbonate_file)
            print(f"  写出 {carbonate_file.name}", flush=True)
        else:
            print(f"  {carbonate_file.name} 已存在，跳过", flush=True)

        if args.force or not strat_file.exists():
            data = {v: window_mean(u, v, years, DEPTHS) for v, u in strat_sources.items()}
            ds = xr.Dataset(data, attrs={"model": model, "experiment": experiment,
                "member_id": "r1i1p1f1", "window": f"{years[0]}-{years[1]}",
                "requested_depths_m": "0, 200, 1000"})
            tmp = strat_file.with_suffix(".tmp.nc")
            ds.to_netcdf(tmp, encoding={v: {"zlib": True, "complevel": 4} for v in data})
            tmp.replace(strat_file)
            print(f"  写出 {strat_file.name}", flush=True)
        else:
            print(f"  {strat_file.name} 已存在，跳过", flush=True)

        manifest.append({"model": model, "experiment": experiment,
            "window": [years[0], years[1]],
            "carbonate_file": str(carbonate_file.relative_to(root)),
            "stratification_file": str(strat_file.relative_to(root)),
            "carbonate_sources": carbonate_sources, "stratification_sources": strat_sources})

(root / "data/preprocessing/mechanism_30yr_manifest.json").write_text(
    json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
)
print("\n30年窗口数据全部完成", flush=True)
