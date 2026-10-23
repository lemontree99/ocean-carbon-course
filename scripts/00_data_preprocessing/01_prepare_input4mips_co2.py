#!/usr/bin/env python
"""下载并提取官方 CMIP6 input4MIPs 全球年均大气 CO2 浓度。"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from urllib.request import urlretrieve

import numpy as np
import pandas as pd
import xarray as xr


root = Path(__file__).resolve().parents[2]
out = root / "data" / "obs" / "co2"
raw = out / "input4mips"
raw.mkdir(parents=True, exist_ok=True)

files = {
    "historical": {
        "name": "mole-fraction-of-carbon-dioxide-in-air_input4MIPs_GHGConcentrations_CMIP_UoM-CMIP-1-2-0_gr1-GMNHSH_0000-2014.nc",
        "url": "https://zenodo.org/api/records/6386979/files/data/raw/cmip6-prescribed-concs/mole-fraction-of-carbon-dioxide-in-air_input4MIPs_GHGConcentrations_CMIP_UoM-CMIP-1-2-0_gr1-GMNHSH_0000-2014.nc/content",
        "md5": "63bf264ef75413ebb365218130fd00c2",
        "file_start": 0,
        "file_end": 2014,
        "keep_start": 1850,
        "keep_end": 2014,
        "source": "input4MIPs UoM-CMIP-1-2-0 historical",
        "doi": "10.22033/ESGF/input4MIPs.1118",
        "reference_doi": "10.5194/gmd-10-2057-2017",
    },
    "ssp585": {
        "name": "mole-fraction-of-carbon-dioxide-in-air_input4MIPs_GHGConcentrations_ScenarioMIP_UoM-REMIND-MAGPIE-ssp585-1-2-1_gr1-GMNHSH_2015-2500.nc",
        "url": "https://zenodo.org/api/records/6386979/files/data/raw/cmip6-prescribed-concs/mole-fraction-of-carbon-dioxide-in-air_input4MIPs_GHGConcentrations_ScenarioMIP_UoM-REMIND-MAGPIE-ssp585-1-2-1_gr1-GMNHSH_2015-2500.nc/content",
        "md5": "61867afcfca7a7f0adb896a59d86bb26",
        "file_start": 2015,
        "file_end": 2500,
        "keep_start": 2015,
        "keep_end": 2100,
        "source": "input4MIPs UoM-REMIND-MAGPIE-ssp585-1-2-1",
        "doi": "10.22033/ESGF/input4MIPs.9868",
        "reference_doi": "10.5194/gmd-13-3571-2020",
    },
}

parts = []
for experiment, info in files.items():
    path = raw / info["name"]
    if not path.exists():
        print(f"下载 {experiment}: {info['url']}")
        urlretrieve(info["url"], path)

    digest = hashlib.md5(path.read_bytes()).hexdigest()
    assert digest == info["md5"], f"{path.name}: MD5不匹配"

    ds = xr.open_dataset(path, decode_times=False)
    var = ds["mole_fraction_of_carbon_dioxide_in_air"]
    assert var.attrs["units"] == "1.e-6"
    assert var.sizes["sector"] == 3 and "0: Global" in ds.sector.attrs["ids"]

    years = np.arange(info["file_start"], info["file_end"] + 1)
    assert len(years) == var.sizes["time"]
    part = pd.DataFrame(
        {
            "year": years,
            # 文件数值以 1e-6 摩尔分数为单位，数值本身即 ppm。
            "co2_ppm": var.isel(sector=0).values,
            "source": info["source"],
        }
    ).query("@info['keep_start'] <= year <= @info['keep_end']")
    parts.append(part)
    ds.close()

co2 = pd.concat(parts, ignore_index=True)
assert np.array_equal(co2.year, np.arange(1850, 2101))
assert co2.co2_ppm.notna().all()

co2.to_csv(out / "co2_concentration_1850_2100.csv", index=False, float_format="%.8f")
co2.query("year <= 2014")[["year", "co2_ppm"]].to_csv(
    out / "co2_historical_annual_1850_2014.txt",
    sep=" ",
    index=False,
    header=False,
    float_format="%.8f",
)
co2.query("year >= 2015")[["year", "co2_ppm"]].to_csv(
    out / "co2_ssp585_annual_2015_2100.txt",
    sep=" ",
    index=False,
    header=False,
    float_format="%.8f",
)

metadata = {
    "description": "CMIP6 input4MIPs global annual-mean atmospheric CO2 mole fraction",
    "variable": "mole_fraction_of_carbon_dioxide_in_air",
    "sector": "0: Global",
    "native_units": "1.e-6 mole fraction",
    "output_units": "ppm",
    "retrieval_record": "https://doi.org/10.5281/zenodo.6386979",
    "retrieval_note": "Exact input4MIPs NetCDF copies archived with IPCC AR6 processing files",
    "files": files,
}
(out / "input4mips_source_metadata.json").write_text(
    json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8"
)

check = co2[co2.year.isin([1850, 1940, 1941, 2014, 2015, 2100])].copy()
check["annual_change_ppm"] = co2.co2_ppm.diff().loc[check.index]
print(check.to_string(index=False))
print(f"\n已写出：{out / 'co2_concentration_1850_2100.csv'}")
