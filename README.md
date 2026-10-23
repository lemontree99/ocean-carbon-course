# CMIP6 海洋碳汇与吸收效率课程项目

使用 GFDL-ESM4、MPI-ESM1-2-LR 和 CanESM5 分析全球碳汇、空间分布、吸收效率、机制量级与DIC库存收支。本仓库来自课程最终提交包 `submission.zip`。

## 报告内容与代码

| 内容 | 代码目录 |
|---|---|
| 输入检查和必要下载 | `scripts/00_data_preprocessing/` |
| 全球海洋碳汇与piControl漂移 | `scripts/01_global_ocean_sink/` |
| 空间分布与局地漂移检验 | `scripts/02_spatial_distribution/` |
| 碳汇效率 | `scripts/03_sink_efficiency/` |
| 缓冲与变暖机制量级 | `scripts/04_mechanisms/` |
| DIC库存自洽性检查 | `scripts/05_self_consistency/` |

各目录包含独立的计算和绘图脚本。正式图、统计表和分节分析报告保存在`results/`，详细对应关系见`scripts/CODE_MAP.md`。

## 运行环境

```powershell
conda env create -f environment.yml
conda activate ocean-carbon-miniproject
python scripts/run_report_pipeline.py --dry-run
python scripts/run_report_pipeline.py
```

## 数据

原始CMIP6和观测文件未放入提交包。输入文件结构、变量和单位见`data_manifest/数据目录与格式说明.md`；观测清单和30年机理数据下载来源也保存在`data_manifest/`。将数据放回说明中的相对路径后，运行`00_validate_inputs.py`即可检查完整性。

仓库提供已计算的表格、图件及部分中间数据，**不包含完整原始数据，克隆后不能直接完成全流程重算**。原始数据须从原提供方获取并遵守其使用条款。最终报告及预估表见 `docs/`。

## 验证与使用边界

发布时检查 Python 语法及流水线 dry-run，未重跑全部计算。机制缩放为一阶代理，不是严格因果归因；结果作为课程研究快照使用。

## 许可

原创代码采用 MIT 许可（见 LICENSE）。第三方数据、底图及课程资料遵循各自原始条款，不因本代码许可而重新授权。
