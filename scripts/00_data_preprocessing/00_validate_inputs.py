#!/usr/bin/env python
"""检查最终报告计算所需的数据文件和关键变量。"""

from pathlib import Path

import xarray as xr


root = Path(__file__).resolve().parents[2]
raw = root / "data/raw"
models = ["GFDL-ESM4", "MPI-ESM1-2-LR", "CanESM5"]
required = {
    "historical": {
        "areacello.nc": "areacello", "fgco2.nc": "fgco2",
        "dissic_slices.nc": "dissic", "volcello.nc": "volcello",
        "mechanism_carbonate_surface_30yr_mean.nc": "dissic",
        "mechanism_stratification_30yr_mean.nc": "thetao",
    },
    "ssp585": {
        "fgco2.nc": "fgco2", "dissic_slices.nc": "dissic",
        "volcello.nc": "volcello",
        "mechanism_carbonate_surface_30yr_mean.nc": "dissic",
        "mechanism_stratification_30yr_mean.nc": "thetao",
    },
    "piControl": {"fgco2.nc": "fgco2"},
}

missing = []
for model in models:
    for experiment, files in required.items():
        for name, variable in files.items():
            path = raw / model / experiment / name
            if not path.exists():
                missing.append(str(path.relative_to(root)))
                continue
            with xr.open_dataset(path, decode_times=False) as ds:
                if variable not in ds:
                    missing.append(f"{path.relative_to(root)}: 缺少变量 {variable}")

other = [
    "data/obs/gcb/gcb2025_ocean_sink.csv",
    "data/obs/landschutzer/spco2_1982-2015_MPI_SOM-FFN_v2016.nc",
    "data/obs/shapefiles/ne_110m_coastline.geojson",
    "data/obs/co2/co2_concentration_1850_2100.csv",
    "data/preprocessing/pangeo-cmip6.csv",
]
missing.extend(path for path in other if not (root / path).exists())

if missing:
    raise SystemExit("缺少报告输入：\n- " + "\n- ".join(missing))

print(f"输入检查通过：3个模式，{sum(len(v) for v in required.values()) * 3}个模式文件及5个公共输入。")
