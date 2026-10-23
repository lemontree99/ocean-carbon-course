# 最终报告代码对应关系

| 报告内容 | 计算代码 | 绘图代码 | 正式输出 |
|---|---|---|---|
| 输入检查 | `00_data_preprocessing/00_validate_inputs.py` | — | 终端检查结果 |
| 全球海洋碳汇与piControl稳定性 | `01_global_ocean_sink/01_compute.py` | `01_global_ocean_sink/02_plot.py` | `results/01_global_ocean_sink/` |
| 空间分布、FDR漂移检验与偏差 | `02_spatial_distribution/01_prepare_fields.py`、`02_drift_correction.py` | `02_spatial_distribution/03_plot.py` | `results/02_spatial_distribution/` |
| 碳汇效率与大气CO₂ | `03_sink_efficiency/01_compute.py` | `02_plot.py`、`03_plot_co2.py` | `results/03_sink_efficiency/` |
| 缓冲能力与变暖机制量级 | `04_mechanisms/01_compute.py` | `04_mechanisms/02_plot.py` | `results/04_mechanisms/` |
| DIC库存—海气通量自洽性 | `05_self_consistency/01_compute.py` | `05_self_consistency/02_plot.py` | `results/05_self_consistency/` |
| DOCX/PDF排版 | `report_build/` | — | `docs/02_final_report/` |

## 完整运行

```powershell
python scripts/run_report_pipeline.py --list
python scripts/run_report_pipeline.py --dry-run
python scripts/run_report_pipeline.py
```

可以用`--only 11,12`只运行指定步骤，或用`--from-step 11`从机制分析开始。

下载脚本不进入日常重算流水线：`01_prepare_input4mips_co2.py`准备官方CO₂浓度，`02_download_mechanism_window_means.py`准备两个30年窗口的碳酸盐与温盐场。已有文件会自动跳过。
