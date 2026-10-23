# 02 空间分布（报告 3.2.3、4.2、附图 2、4、5）

## 脚本与顺序

| 步骤 | 脚本 | 作用 |
|---|---|---|
| 1 | `01_prepare_fields.py` | 计算 1982–2014 年模式与 Landschützer 的平均通量场，统一到 2° 共同网格；输出原始空间指标 |
| 2 | `02_drift_correction.py` | 用分支对齐 piControl 计算格点漂移，Newey–West 标准误 + FDR 0.05；只订正 FDR 显著且累计幅度 ≥ 0.10 mol C m⁻² yr⁻¹ 的格点；输出订正场和订正后指标 |
| 3 | `03_plot.py` | 绘制漂移订正区域、raw p=0.05 与 FDR 对比、气候态和偏差图 |

```powershell
python scripts/02_spatial_distribution/01_prepare_fields.py
python scripts/02_spatial_distribution/02_drift_correction.py
python scripts/02_spatial_distribution/03_plot.py
```

## 输入

- `data/raw/<model>/historical/fgco2.nc`、`areacello.nc`
- `data/raw/<model>/piControl/fgco2.nc`
- `data/obs/` 下的 Landschützer v2016 月产品

## 输出

- `results/02_spatial_distribution/intermediate/spatial_fields_2deg.npz`
- `results/02_spatial_distribution/intermediate/spatial_drift_fields_2deg.npz`
- `results/02_spatial_distribution/intermediate/spatial_fields_2deg_drift_corrected.npz`
- `results/02_spatial_distribution/tables/global_spatial_metrics.csv`
- `results/02_spatial_distribution/tables/regional_spatial_metrics.csv`
- `results/02_spatial_distribution/tables/spatial_drift_summary.csv`
- `results/02_spatial_distribution/tables/spatial_drift_regional_summary.csv`
- `results/02_spatial_distribution/tables/global_spatial_metrics_drift_corrected.csv`
- `results/02_spatial_distribution/tables/regional_spatial_metrics_drift_corrected.csv`
- `results/02_spatial_distribution/tables/spatial_raw_p05_native_gridcells.csv`
- `results/02_spatial_distribution/figures/spatial_drift_required_regions.png/pdf`
- `results/02_spatial_distribution/figures/spatial_raw_p05_vs_fdr.png/pdf`
- `results/02_spatial_distribution/figures/spatial_flux_climatology_drift_corrected.png/pdf`
- `results/02_spatial_distribution/figures/spatial_bias_maps_drift_corrected.png/pdf`

## 口径

- 报告 4.2 的型相关、RMSE、区域偏差使用 **订正后** 的
  `global_spatial_metrics_drift_corrected.csv` / `regional_spatial_metrics_drift_corrected.csv`。
- 订正只针对局地显著漂移格点，不减去 piControl 完整空间气候态，避免删除赤道自然放气和高纬自然吸收。
