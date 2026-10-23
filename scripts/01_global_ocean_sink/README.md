# 01 全球海洋碳汇（报告 3.2.1、4.1、附图 1、附图 3）

## 脚本

| 脚本 | 作用 |
|---|---|
| `01_compute.py` | 读 `fgco2 × areacello`，按模式原生网格积分，汇总为 Pg C yr⁻¹；用分支对齐的 piControl 线性趋势做时间漂移订正；输出年度序列、piControl 稳定性、时期均值和 GCB 对比 |
| `02_plot.py` | 绘制 1850–2100 年全球碳汇与 GCB 对比，以及三模式 piControl 稳定性图；`--series corrected\|raw` 控制用订正或未订正序列 |

## 运行

```powershell
python scripts/01_global_ocean_sink/01_compute.py
python scripts/01_global_ocean_sink/02_plot.py --series raw
python scripts/01_global_ocean_sink/02_plot.py
```

## 输入

- `data/raw/<model>/historical/fgco2.nc`、`areacello.nc`
- `data/raw/<model>/ssp585/fgco2.nc`
- `data/raw/<model>/piControl/fgco2.nc`
- `data/obs/gcb/gcb2025_ocean_sink.csv`

## 输出

- `results/01_global_ocean_sink/tables/global_ocean_sink_annual.csv`
- `results/01_global_ocean_sink/tables/picontrol_global_sink_annual.csv`
- `results/01_global_ocean_sink/tables/picontrol_correction_summary.csv`
- `results/01_global_ocean_sink/tables/picontrol_stability_test.csv`
- `results/01_global_ocean_sink/tables/picontrol_subperiod_trends.csv`
- `results/01_global_ocean_sink/tables/period_summary.csv`
- `results/01_global_ocean_sink/tables/gcb_comparison_stats.csv`
- `results/01_global_ocean_sink/figures/global_ocean_sink_timeseries[_raw].png/pdf`
- `results/01_global_ocean_sink/figures/picontrol_stability_check.png/pdf`

## 口径

- 报告 4.1 和附图 3 使用 **未订正** 序列（`sink_raw_PgC_yr`，`--series raw`）。
- 4.3 的效率计算使用 **订正后** 的全球碳汇序列（`sink_corrected_PgC_yr`）；4.2 使用 `02_spatial_distribution` 的空间订正场。
- 漂移订正只扣分支以后目标量的线性趋势，不扣控制试验平均气候态和逐年波动。
- 计算脚本读取约 2 GB 月通量，运行时间较长；绘图只读 CSV，几秒完成。
