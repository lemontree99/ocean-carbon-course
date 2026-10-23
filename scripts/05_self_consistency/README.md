# 分析5：DIC库存—海气通量自洽性检查

```powershell
python scripts/05_self_consistency/01_compute.py
python scripts/05_self_consistency/02_plot.py
```

`01_compute.py`以`dissic × volcello`计算全球DIC库存变化，并将`fgco2 × areacello`在相同端点之间积分。库存端点为月中平均，因此首尾月份的通量各取一半权重。主结果不做piControl订正，与最终报告表格口径一致。

输入均位于`data/raw/<model>/<experiment>/`：`dissic_slices.nc`、`volcello.nc`、`fgco2.nc`，面积读取historical的`areacello.nc`。

输出包括`raw_closure.csv`、`DIC_inventory_endpoints.csv`、计算元数据和`dic_fgco2_closure.png/.pdf`。DIC并非海洋全部碳库，残差还包含河流输入、沉积交换以及DOC/POC/PIC等边界项。
