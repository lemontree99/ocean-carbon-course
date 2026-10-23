# 数据准备

| 文件 | 用途 |
|---|---|
| `00_validate_inputs.py` | 检查最终报告所需模式文件、变量和公共资料；不修改数据 |
| `01_prepare_input4mips_co2.py` | 下载并拼接官方CMIP6 input4MIPs historical与SSP5-8.5 CO₂浓度 |
| `02_download_mechanism_window_means.py` | 从Pangeo CMIP6提取1985–2014和2071–2100年日数加权平均的表层DIC/TA及0/200/1000 m温盐 |

正式分析直接读取`data/raw/`和`data/obs/`。下载来源记录在`data/preprocessing/mechanism_30yr_manifest.json`和`data/obs/co2/input4mips_source_metadata.json`。
