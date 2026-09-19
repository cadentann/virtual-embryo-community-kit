from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def env() -> dict[str, str]:
    value = os.environ.copy()
    srcs = [str(ROOT / "modules/task1_walkthrough/src"), str(ROOT / "modules/data_inspector/src")]
    value["PYTHONPATH"] = os.pathsep.join(srcs + ([value["PYTHONPATH"]] if value.get("PYTHONPATH") else []))
    return value


def test_synthetic_a_to_b_path(tmp_path: Path) -> None:
    output = tmp_path / "task1_synthetic_prediction.h5ad"
    result = subprocess.run(
        [sys.executable, "-m", "vec_t1_walkthrough.cli", "demo", "--output", str(output), "--seed", "20260911"],
        cwd=tmp_path, env=env(), capture_output=True, text=True, check=False,
    )
    assert result.returncode == 0, result.stderr
    sidecar = Path(str(output) + ".selection.json")
    assert output.is_file() and sidecar.is_file()
    selection = json.loads(sidecar.read_text())
    assert selection["selection_policy"] == "uniform_without_replacement"
    report_dir = tmp_path / "inspection"
    result = subprocess.run(
        [sys.executable, "-m", "vec_data_inspector.cli", str(output), "--report-dir", str(report_dir)],
        cwd=tmp_path, env=env(), capture_output=True, text=True, check=False,
    )
    assert result.returncode == 0, result.stderr
    assert (report_dir / "inspection_report.md").is_file()
    assert (report_dir / "inspection_report.json").is_file()
