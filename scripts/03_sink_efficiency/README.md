# 03 碳汇效率（报告 3.3、4.3、附图 6、附图 7）

## 脚本

| 脚本 | 作用 |
|---|---|
| `01_compute.py` | 读全球碳汇年度序列与大气 CO₂，计算年度效率、累计效率和 20/30/40 年窗口效率；输出汇总表和窗口敏感性 |
| `02_plot.py` | 绘制年度效率与 30 年窗口效率 |
| `03_plot_co2.py` | 绘制 1850–2100 年大气 CO₂ 浓度 |

```powershell
python scripts/03_sink_efficiency/01_compute.py
python scripts/03_sink_efficiency/02_plot.py
python scripts/03_sink_efficiency/03_plot_co2.py
```

## 输入

- `results/01_global_ocean_sink/tables/global_ocean_sink_annual.csv`
- `data/obs/co2/co2_concentration_1850_2100.csv`

## 输出

- `results/03_sink_efficiency/tables/sink_efficiency_annual.csv`
- `results/03_sink_efficiency/tables/sink_efficiency_summary.csv`
- `results/03_sink_efficiency/tables/efficiency_window_sensitivity.csv`
- `results/03_sink_efficiency/tables/calculation_metadata.json`
- `results/03_sink_efficiency/figures/sink_efficiency_timeseries.png/pdf`
- `results/03_sink_efficiency/figures/atmospheric_co2_concentration.png/pdf`

## 口径

- 效率定义为累计碳汇 / 同期大气 CO₂ 增量（Pg C ppm⁻¹）；年度效率用逐年增量。
- 主结果使用 piControl 趋势订正后的全球碳汇；`efficiency_window_sensitivity.csv` 同时给出 raw、
  drift-corrected 和 control-baseline-removed 三套序列。
- 报告 4.3 的 60%–62% 对应 1985–2014 与 2071–2100 两个 30 年窗口，不是 1950–2100 的全程降幅。
