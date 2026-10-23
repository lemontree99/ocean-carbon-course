#!/usr/bin/env python
"""绘制官方 input4MIPs 大气 CO2 浓度时间序列。"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


root = Path(__file__).resolve().parents[2]
data = pd.read_csv(root / "data" / "obs" / "co2" / "co2_concentration_1850_2100.csv")
figures = root / "results" / "03_sink_efficiency" / "figures"
figures.mkdir(parents=True, exist_ok=True)

historical = data[data.year <= 2014]
ssp585 = data[data.year >= 2015]

fig, ax = plt.subplots(figsize=(10.5, 5.4), constrained_layout=True)
ax.plot(historical.year, historical.co2_ppm, color="#0072B2", linewidth=2.2, label="Historical")
ax.plot(ssp585.year, ssp585.co2_ppm, color="#D55E00", linewidth=2.2, label="SSP5-8.5")
ax.axvline(2015, color="0.4", linestyle="--", linewidth=1.0)
ax.set(
    xlim=(1850, 2100),
    ylim=(250, 1200),
    xlabel="Year",
    ylabel="Atmospheric CO$_2$ concentration (ppm)",
    title="Prescribed atmospheric CO$_2$ concentration in CMIP6",
)
ax.set_xticks([1850, 1900, 1950, 2000, 2015, 2050, 2100])
ax.grid(axis="y", color="0.90", linewidth=0.7)
ax.legend(frameon=False, loc="upper left")

fig.savefig(figures / "atmospheric_co2_concentration.png", dpi=300, bbox_inches="tight")
fig.savefig(figures / "atmospheric_co2_concentration.pdf", bbox_inches="tight")
plt.close(fig)
print(f"图已保存到：{figures}")
