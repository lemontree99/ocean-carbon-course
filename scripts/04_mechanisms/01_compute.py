#!/usr/bin/env python
"""量化缓冲能力减弱与变暖相关物理效应的量级。"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import gsw
import numpy as np
import pandas as pd
import PyCO2SYS as pyco2
import xarray as xr


MODELS = ["GFDL-ESM4", "MPI-ESM1-2-LR", "CanESM5"]

parser = argparse.ArgumentParser()
parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[2])
args = parser.parse_args()
root = args.project_root.resolve()
raw = root / "data" / "raw"
out = root / "results" / "04_mechanisms" / "tables"
out.mkdir(parents=True, exist_ok=True)


def mean(x, area, mask):
    """共同有效海区的面积加权平均。"""
    w = np.where(mask, area, 0.0)
    return float(np.nansum(np.where(mask, x * w, 0.0)) / np.sum(w))


def weiss_k0(temp, sal):
    """Weiss (1974) CO2溶解度，mol kg-1 atm-1。"""
    tk = temp + 273.15
    return np.exp(
        -60.2409
        + 93.4517 * (100 / tk)
        + 23.3585 * np.log(tk / 100)
        + sal * (0.023517 - 0.023656 * (tk / 100) + 0.0047036 * (tk / 100) ** 2)
    )


def carbonate(dic, alk, temp, sal, rho, mask):
    """在有效格点用TA+DIC求Revelle因子和pCO2。"""
    result = pyco2.sys(
        par1=alk[mask] / rho[mask] * 1e6,
        par2=dic[mask] / rho[mask] * 1e6,
        par1_type=1,
        par2_type=2,
        salinity=sal[mask],
        temperature=temp[mask],
        pressure=0,
        opt_pH_scale=1,
        opt_k_carbonic=10,
        opt_k_bisulfate=1,
    )
    revelle = np.full(dic.shape, np.nan)
    pco2 = np.full(dic.shape, np.nan)
    revelle[mask] = result["revelle_factor"]
    pco2[mask] = result["pCO2"]
    return revelle, pco2


eff = pd.read_csv(root / "results" / "03_sink_efficiency" / "tables" / "efficiency_window_sensitivity.csv")
eff = eff[(eff.series == "drift_corrected") & (eff.window_years == 30)].set_index("model")
carbonate_rows, warming_rows, mechanism_rows, regional_rows, state_rows = [], [], [], [], []

for model in MODELS:
    area = xr.open_dataset(raw / model / "historical" / "areacello.nc").areacello.load().values.astype(float)
    hc = xr.open_dataset(raw / model / "historical" / "mechanism_carbonate_surface_30yr_mean.nc").load()
    fc = xr.open_dataset(raw / model / "ssp585" / "mechanism_carbonate_surface_30yr_mean.nc").load()
    hs = xr.open_dataset(raw / model / "historical" / "mechanism_stratification_30yr_mean.nc").load()
    fs = xr.open_dataset(raw / model / "ssp585" / "mechanism_stratification_30yr_mean.nc").load()

    # 部分模式的表层文件仍保留长度为 1 的 lev 维，统一压成二维海表场。
    dic_h, alk_h = hc.dissic.squeeze(drop=True).values, hc.talk.squeeze(drop=True).values
    dic_f, alk_f = fc.dissic.squeeze(drop=True).values, fc.talk.squeeze(drop=True).values
    temp_h, sal_h = hs.thetao.values, hs.so.values
    temp_f, sal_f = fs.thetao.values, fs.so.values
    th, sh, tf, sf = temp_h[0], sal_h[0], temp_f[0], sal_f[0]

    mask = (
        np.isfinite(area + dic_h + alk_h + dic_f + alk_f + th + sh + tf + sf)
        & (area > 0)
        & (dic_h > 0.5) & (dic_f > 0.5)
        & (alk_h > 0.5) & (alk_f > 0.5)
        & (sh > 2) & (sh < 45) & (sf > 2) & (sf < 45)
        & (th > -3) & (th < 40) & (tf > -3) & (tf < 40)
    )
    rho_h = gsw.rho(sh, th, 0)
    rho_f = gsw.rho(sf, tf, 0)
    rh, ph = carbonate(dic_h, alk_h, th, sh, rho_h, mask)
    rf, pf = carbonate(dic_f, alk_f, tf, sf, rho_f, mask)
    # 主缓冲反事实：固定1985-2014年平均温度，DIC、TA和盐度使用各自时期数据。
    rho_f_buffer = gsw.rho(sf, th, 0)
    rf_buffer, pf_buffer = carbonate(dic_f, alk_f, th, sf, rho_f_buffer, mask)
    # 固定历史窗口平均温盐的结果仅保留为敏感性对照。
    rf_chem, pf_chem = carbonate(dic_f, alk_f, th, sh, rho_h, mask)
    cmask = (
        mask
        & np.isfinite(rh + ph + rf + pf + rf_buffer + pf_buffer + rf_chem + pf_chem)
        & (rh > 0) & (rf_buffer > 0) & (rf_chem > 0)
    )

    r_h = mean(rh, area, cmask)
    r_f = mean(rf, area, cmask)
    r_f_buffer = mean(rf_buffer, area, cmask)
    r_f_chem = mean(rf_chem, area, cmask)
    p_h = mean(ph, area, cmask)
    p_f = mean(pf, area, cmask)
    p_f_buffer = mean(pf_buffer, area, cmask)
    p_f_chem = mean(pf_chem, area, cmask)
    buffer_retention = mean(1 / rf_buffer, area, cmask) / mean(1 / rh, area, cmask)
    buffer_retention_fixed_ts = mean(1 / rf_chem, area, cmask) / mean(1 / rh, area, cmask)
    beta_h = dic_h / (rh * ph)
    beta_f = dic_f / (rf_buffer * pf_buffer)
    beta_f_fixed_ts = dic_f / (rf_chem * pf_chem)

    for period, temp, sal, alk in [
        ("1985-2014", th, sh, alk_h),
        ("2071-2100", tf, sf, alk_f),
    ]:
        state_rows.append({
            "model": model,
            "period": period,
            "surface_temperature_degC": mean(temp, area, cmask),
            "surface_salinity": mean(sal, area, cmask),
            "surface_TA_mol_m3": mean(alk, area, cmask),
        })

    carbonate_rows.append({
        "model": model,
        "valid_ocean_area_percent": 100 * np.sum(area[cmask]) / np.sum(area[area > 0]),
        "surface_DIC_historical_mol_m3": mean(dic_h, area, cmask),
        "surface_DIC_late_mol_m3": mean(dic_f, area, cmask),
        "surface_TA_historical_mol_m3": mean(alk_h, area, cmask),
        "surface_TA_late_mol_m3": mean(alk_f, area, cmask),
        "Revelle_historical": r_h,
        "Revelle_late_actual_TS": r_f,
        "Revelle_late_fixed_historical_T_actual_S": r_f_buffer,
        "Revelle_late_fixed_historical_TS": r_f_chem,
        "Revelle_change_percent_main_fixed_historical_T": 100 * (r_f_buffer / r_h - 1),
        "Revelle_change_percent_actual_TS": 100 * (r_f / r_h - 1),
        "Revelle_change_percent_fixed_TS": 100 * (r_f_chem / r_h - 1),
        "carbonate_pCO2_historical_uatm": p_h,
        "carbonate_pCO2_late_actual_TS_uatm": p_f,
        "carbonate_pCO2_late_fixed_historical_T_actual_S_uatm": p_f_buffer,
        "carbonate_pCO2_late_fixed_historical_TS_uatm": p_f_chem,
        "inverse_Revelle_retention": buffer_retention,
        "buffer_loss_fraction": 1 - buffer_retention,
        "inverse_Revelle_retention_fixed_historical_TS": buffer_retention_fixed_ts,
        "buffer_loss_fraction_fixed_historical_TS": 1 - buffer_retention_fixed_ts,
        "full_carbonate_sensitivity_historical_mol_m3_uatm": mean(beta_h, area, cmask),
        "full_carbonate_sensitivity_late_fixed_historical_T_actual_S_mol_m3_uatm": mean(beta_f, area, cmask),
        "full_carbonate_sensitivity_late_fixed_historical_TS_mol_m3_uatm": mean(beta_f_fixed_ts, area, cmask),
        "full_carbonate_sensitivity_change_percent": 100 * (mean(beta_f, area, cmask) / mean(beta_h, area, cmask) - 1),
    })

    # 直接热力学效应：只改变温度，盐度保持历史窗口平均水平。
    k_h = weiss_k0(th, sh)
    k_f_temp = weiss_k0(tf, sh)
    k_f_actual = weiss_k0(tf, sf)
    k_ratio = mean(k_f_temp, area, cmask) / mean(k_h, area, cmask)

    strat = {}
    for deep_index, label in [(1, "200m"), (2, "1000m")]:
        dz = float(hs.lev.values[deep_index] - hs.lev.values[0])
        sigma_h = gsw.sigma0(sal_h, temp_h)
        sigma_f = gsw.sigma0(sal_f, temp_f)
        sigma_f_temp = gsw.sigma0(sal_h, temp_f)
        n2_h = 9.81 / 1025 * (sigma_h[deep_index] - sigma_h[0]) / dz
        n2_f = 9.81 / 1025 * (sigma_f[deep_index] - sigma_f[0]) / dz
        n2_f_temp = 9.81 / 1025 * (sigma_f_temp[deep_index] - sigma_f_temp[0]) / dz
        nmask = cmask & np.isfinite(n2_h + n2_f + n2_f_temp) & (n2_h > 0) & (n2_f > 0) & (n2_f_temp > 0)
        strat[label] = {
            "h": mean(n2_h, area, nmask),
            "f": mean(n2_f, area, nmask),
            "f_temp": mean(n2_f_temp, area, nmask),
            "area": 100 * np.sum(area[nmask]) / np.sum(area[area > 0]),
            "n2_h": n2_h,
            "n2_f_temp": n2_f_temp,
            "mask": nmask,
        }

    # 主估算采用0–200 m。假设垂向输送能力与N²成反比，仅作为一阶量级代理。
    transport_retention = strat["200m"]["h"] / strat["200m"]["f_temp"]
    warming_retention = k_ratio * transport_retention
    transport_retention_1000 = strat["1000m"]["h"] / strat["1000m"]["f_temp"]
    warming_retention_1000 = k_ratio * transport_retention_1000
    e_h = float(eff.loc[model, "efficiency_ending_2014_PgC_per_ppm"])
    e_f = float(eff.loc[model, "efficiency_ending_2100_PgC_per_ppm"])
    buffer_loss = e_h * (1 - buffer_retention)
    buffer_loss_fixed_ts = e_h * (1 - buffer_retention_fixed_ts)
    solubility_loss = e_h * (1 - k_ratio)
    stratification_loss = e_h * (1 - transport_retention)
    warming_loss = e_h * (1 - warming_retention)

    latitude = (
        np.broadcast_to(hc.lat.values[:, None], area.shape)
        if "lat" in hc.coords
        else hc.latitude.values
    )
    for region, lo, hi in [
        ("Southern high latitudes", -90, -40),
        ("Low and mid latitudes", -40, 40),
        ("Northern high latitudes", 40, 90.1),
    ]:
        rmask = cmask & (latitude >= lo) & (latitude < hi)
        nmask = rmask & strat["200m"]["mask"]
        bkeep = mean(1 / rf_buffer, area, rmask) / mean(1 / rh, area, rmask)
        bkeep_fixed_ts = mean(1 / rf_chem, area, rmask) / mean(1 / rh, area, rmask)
        kkeep = mean(k_f_temp, area, rmask) / mean(k_h, area, rmask)
        vkeep = mean(strat["200m"]["n2_h"], area, nmask) / mean(strat["200m"]["n2_f_temp"], area, nmask)
        wkeep = kkeep * vkeep
        regional_rows.append({
            "model": model,
            "region": region,
            "latitude_bounds_degrees": f"{lo:g} to {hi:g}",
            "valid_ocean_area_m2": np.sum(area[rmask]),
            "Revelle_historical": mean(rh, area, rmask),
            "Revelle_late_fixed_historical_T_actual_S": mean(rf_buffer, area, rmask),
            "Revelle_late_fixed_historical_TS": mean(rf_chem, area, rmask),
            "buffer_loss_percent": 100 * (1 - bkeep),
            "buffer_loss_percent_fixed_historical_TS": 100 * (1 - bkeep_fixed_ts),
            "surface_warming_degC": mean(tf - th, area, rmask),
            "solubility_loss_percent": 100 * (1 - kkeep),
            "N2_0_200m_temperature_only_change_percent": 100 * (1 / vkeep - 1),
            "warming_joint_loss_percent": 100 * (1 - wkeep),
            "larger_standalone_mechanism": "buffer weakening" if (1 - bkeep) > (1 - wkeep) else "warming-related physical effect",
        })

    warming_rows.append({
        "model": model,
        "surface_temperature_historical_degC": mean(th, area, cmask),
        "surface_temperature_late_degC": mean(tf, area, cmask),
        "surface_warming_degC": mean(tf - th, area, cmask),
        "surface_salinity_change": mean(sf - sh, area, cmask),
        "K0_historical_mol_kg_atm": mean(k_h, area, cmask),
        "K0_late_temperature_only_mol_kg_atm": mean(k_f_temp, area, cmask),
        "K0_late_actual_TS_mol_kg_atm": mean(k_f_actual, area, cmask),
        "solubility_retention": k_ratio,
        "solubility_loss_fraction": 1 - k_ratio,
        "N2_0_200m_historical_s_2": strat["200m"]["h"],
        "N2_0_200m_late_actual_TS_s_2": strat["200m"]["f"],
        "N2_0_200m_late_temperature_only_s_2": strat["200m"]["f_temp"],
        "N2_0_200m_temperature_only_change_percent": 100 * (strat["200m"]["f_temp"] / strat["200m"]["h"] - 1),
        "N2_0_1000m_historical_s_2": strat["1000m"]["h"],
        "N2_0_1000m_late_actual_TS_s_2": strat["1000m"]["f"],
        "N2_0_1000m_late_temperature_only_s_2": strat["1000m"]["f_temp"],
        "N2_0_1000m_temperature_only_change_percent": 100 * (strat["1000m"]["f_temp"] / strat["1000m"]["h"] - 1),
        "valid_stratification_area_0_200m_percent": strat["200m"]["area"],
        "transport_retention_inverse_N2_proxy": transport_retention,
        "warming_joint_retention": warming_retention,
        "warming_joint_loss_fraction": 1 - warming_retention,
        "warming_joint_retention_0_1000m_sensitivity": warming_retention_1000,
        "warming_joint_loss_fraction_0_1000m_sensitivity": 1 - warming_retention_1000,
    })
    mechanism_rows.append({
        "model": model,
        "recent_efficiency_PgC_per_ppm": e_h,
        "late_efficiency_PgC_per_ppm": e_f,
        "actual_total_loss_PgC_per_ppm": e_h - e_f,
        "actual_total_loss_percent": 100 * (1 - e_f / e_h),
        "buffer_estimated_loss_PgC_per_ppm": buffer_loss,
        "buffer_estimated_loss_percent_of_baseline": 100 * buffer_loss / e_h,
        "buffer_fixed_historical_TS_sensitivity_PgC_per_ppm": buffer_loss_fixed_ts,
        "solubility_internal_diagnostic_PgC_per_ppm": solubility_loss,
        "stratification_internal_proxy_PgC_per_ppm": stratification_loss,
        "warming_joint_estimated_loss_PgC_per_ppm": warming_loss,
        "warming_joint_estimated_loss_percent_of_baseline": 100 * warming_loss / e_h,
        "warming_joint_1000m_sensitivity_PgC_per_ppm": e_h * (1 - warming_retention_1000),
        "counterfactual_separation": "buffer fixes historical-window T and uses period-specific S; warming fixes historical-window S and uses period-specific T",
        "larger_standalone_mechanism": "buffer weakening" if buffer_loss > warming_loss else "warming-related physical effect",
    })

pd.DataFrame(carbonate_rows).to_csv(out / "carbonate_buffering_summary.csv", index=False)
pd.DataFrame(warming_rows).to_csv(out / "warming_physical_summary.csv", index=False)
pd.DataFrame(mechanism_rows).to_csv(out / "mechanism_magnitude_summary.csv", index=False)
pd.DataFrame(regional_rows).to_csv(out / "regional_mechanism_summary.csv", index=False)
pd.DataFrame(state_rows).to_csv(out / "state_variable_comparison.csv", index=False)

metadata = {
    "state_windows": {"historical": "1985-2014", "ssp585": "2071-2100"},
    "efficiency_windows": {"historical": "1985-2014", "ssp585": "2071-2100"},
    "carbonate_solver": "PyCO2SYS 1.8.3.4; TA+DIC; Lueker et al. carbonic constants (option 10)",
    "DIC_TA_conversion": "mol m-3 / gsw.rho(S,T,0) * 1e6 -> umol kg-1",
    "buffer_primary": "period-specific DIC, TA and salinity, with temperature fixed at the 1985-2014 mean",
    "buffer_sensitivity": "late-period DIC and TA with both temperature and salinity fixed at historical-window means",
    "solubility": "Weiss (1974); late-period temperature with historical-window salinity",
    "stratification": "N2=(g/rho0)*(sigma0_deep-sigma0_surface)/dz; late-period temperature with historical-window salinity",
    "warming_joint_scaling": "K0_late/K0_historical multiplied by N2_historical/N2_late_temperature_only",
    "regional_bands": ["90S-40S", "40S-40N", "40N-90N"],
    "warning": "The counterfactuals isolate direct temperature and salinity effects, but DIC and TA remain coupled model responses; inverse-N2 transport scaling is an order-of-magnitude proxy, and mechanism estimates are not strictly additive",
}
(out / "calculation_metadata.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
print(pd.DataFrame(mechanism_rows).to_string(index=False))
