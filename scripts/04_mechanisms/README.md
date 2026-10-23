# 分析4：效率下降机制

在项目根目录运行：

```powershell
python scripts/04_mechanisms/01_compute.py
python scripts/04_mechanisms/02_plot.py
```

`01_compute.py`使用1985–2014和2071–2100年平均场，逐格点调用PyCO2SYS计算Revelle因子，用Weiss（1974）公式计算CO₂溶解度，并用GSW状态方程计算0–200 m和0–1000 m的密度层化。缓冲反事实固定历史温度，变暖反事实固定历史盐度。`02_plot.py`读取计算表，生成机制量级及两时期海温、盐度和总碱度对比图。

主要输出位于`results/04_mechanisms/`。计算脚本不依赖当前工作目录，可用`--project-root`指定项目根目录。
