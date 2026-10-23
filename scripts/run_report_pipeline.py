#!/usr/bin/env python
"""按最终报告顺序运行数据检查、五项分析及全部绘图。"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"

STEPS = [
    {"id": 1, "label": "validate_inputs", "script": "00_data_preprocessing/00_validate_inputs.py", "args": []},
    {"id": 2, "label": "global_compute", "script": "01_global_ocean_sink/01_compute.py", "args": []},
    {"id": 3, "label": "global_plot_raw", "script": "01_global_ocean_sink/02_plot.py", "args": ["--series", "raw"],
     "note": "报告 4.1 / 附图3 使用未订正图"},
    {"id": 4, "label": "global_plot_corrected", "script": "01_global_ocean_sink/02_plot.py", "args": []},
    {"id": 5, "label": "spatial_prepare", "script": "02_spatial_distribution/01_prepare_fields.py", "args": []},
    {"id": 6, "label": "spatial_drift", "script": "02_spatial_distribution/02_drift_correction.py", "args": []},
    {"id": 7, "label": "spatial_plot", "script": "02_spatial_distribution/03_plot.py", "args": []},
    {"id": 8, "label": "efficiency_compute", "script": "03_sink_efficiency/01_compute.py", "args": []},
    {"id": 9, "label": "efficiency_plot", "script": "03_sink_efficiency/02_plot.py", "args": []},
    {"id": 10, "label": "co2_plot", "script": "03_sink_efficiency/03_plot_co2.py", "args": []},
    {"id": 11, "label": "mechanism_compute", "script": "04_mechanisms/01_compute.py", "args": []},
    {"id": 12, "label": "mechanism_plot", "script": "04_mechanisms/02_plot.py", "args": []},
    {"id": 13, "label": "closure_compute", "script": "05_self_consistency/01_compute.py", "args": []},
    {"id": 14, "label": "closure_plot", "script": "05_self_consistency/02_plot.py", "args": []},
]


def parse_ids(text: str) -> set[int]:
    if not text:
        return set()
    return {int(item) for item in text.split(",") if item.strip()}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--python", default=sys.executable, help="Python 解释器")
    parser.add_argument("--dry-run", action="store_true", help="只打印命令，不执行")
    parser.add_argument("--list", action="store_true", help="列出步骤后退出")
    parser.add_argument("--only", default="", help="只运行指定编号，逗号分隔，例如 2,3,11")
    parser.add_argument("--from-step", type=int, default=1, help="从指定编号开始")
    args = parser.parse_args()

    if args.list:
        for step in STEPS:
            print(f"{step['id']:>2}  {step['label']:<24} {step['script']}")
        return 0

    only_ids = parse_ids(args.only)
    selected = []
    for step in STEPS:
        if step["id"] < args.from_step and not only_ids:
            continue
        if only_ids and step["id"] not in only_ids:
            continue
        selected.append(step)

    if not selected:
        print("没有选中任何步骤。用 --list 查看编号。")
        return 1

    for step in selected:
        script = SCRIPTS / step["script"]
        if not script.exists():
            print(f"[error] 缺少脚本: {script}")
            return 2
        cmd = [args.python, str(script), *step.get("args", [])]
        print(f"\n[{step['id']:>2}] {step['label']}")
        print("  " + " ".join(cmd))
        if step.get("note"):
            print(f"  note: {step['note']}")
        if args.dry_run:
            continue
        subprocess.run(cmd, cwd=ROOT, check=True)

    print(f"\n完成 {len(selected)} 个步骤。" + ("（dry-run，未执行）" if args.dry_run else ""))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
