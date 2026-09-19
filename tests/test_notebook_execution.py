"""Execute notebooks as a new participant would, using synthetic inputs only."""

from __future__ import annotations

import shutil
import os
import sys
from pathlib import Path

import pytest

@pytest.mark.parametrize("relative_notebook", [
    Path("modules/task1_walkthrough/notebooks/task1_baseline_to_submission.ipynb"),
    Path("modules/data_inspector/notebooks/vec_data_inspector_demo.ipynb"),
])
def test_notebook_is_clean_and_local_execution_ready(tmp_path: Path, relative_notebook: Path):
    root = Path(__file__).resolve().parents[1]
    notebook = root / relative_notebook
    assert notebook.is_file(), f"Configured notebook missing: {notebook}"
    nbformat = pytest.importorskip("nbformat")
    nbclient = pytest.importorskip("nbclient")
    # Run from a disposable copy of the contribution root, as a new clone
    # would. This preserves the release candidate while retaining its relative
    # examples/ and notebooks/ paths.
    candidate_relative = relative_notebook.parents[1]
    candidate_root = root / candidate_relative
    copied_root = tmp_path / candidate_relative.name
    shutil.copytree(candidate_root, copied_root)
    copied_notebook = copied_root / relative_notebook.relative_to(candidate_relative)
    document = nbformat.read(copied_notebook, as_version=4)
    client = nbclient.NotebookClient(document, timeout=180, kernel_name="python3", resources={"metadata": {"path": str(copied_root)}})
    # nbclient inherits process environment through the spawned kernel. Avoid external data paths.
    with pytest.MonkeyPatch.context() as monkeypatch:
        monkeypatch.setenv("VEC_USE_SYNTHETIC_FIXTURES", "1")
        client.execute(cwd=str(copied_root))
