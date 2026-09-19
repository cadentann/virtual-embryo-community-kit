"""Public-interface test support.  Does not import either release candidate."""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pytest


HERE = Path(__file__).resolve().parent
TEST_ROOT = HERE.parent


def require_anndata():
    return pytest.importorskip("anndata", reason="AnnData is required for release acceptance tests")


def file_fingerprint(path: Path) -> dict[str, Any]:
    stat = path.stat()
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return {"sha256": digest.hexdigest(), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def load_config() -> dict[str, Any]:
    path = TEST_ROOT / "blackbox_config.json"
    assert path.is_file(), (
        f"Missing {path}. Copy blackbox_config.example.json and configure only public commands."
    )
    data = json.loads(path.read_text())
    for name in ("release_root", "inspector_command", "walkthrough_submission_command"):
        assert name in data, f"blackbox_config.json is missing {name!r}"
    assert isinstance(data["inspector_command"], list) and data["inspector_command"], "inspector_command must be argv"
    assert isinstance(data["walkthrough_submission_command"], list), "walkthrough_submission_command must be argv"
    root = (TEST_ROOT / data["release_root"]).resolve()
    assert root.is_dir(), f"Configured release_root does not exist: {root}"
    data["_release_root"] = root
    return data


def resolve_argv(argv: list[str], root: Path, **substitutions: str) -> list[str]:
    """Expand deliberate public placeholders and make relative script paths portable."""
    values = {"python": sys.executable, **substitutions}
    result = [str(piece).format(**values) for piece in argv]
    if len(result) > 1 and not result[1].startswith("-"):
        possible_path = root / result[1]
        if possible_path.exists():
            result[1] = str(possible_path)
    return result


def candidate_env(config: dict[str, Any]) -> dict[str, str]:
    """Make source-layout candidates importable without a machine-local venv path.

    The public black-box configuration invokes documented ``python -m`` module
    commands.  This small harness shim makes both ``src/`` trees available from
    a source checkout; an installed release works equally well because its
    installed package takes precedence where appropriate.
    """
    roots = sorted(path for path in config["_release_root"].glob("modules/*/src") if path.is_dir())
    inherited = os.environ.get("PYTHONPATH")
    entries = [str(path) for path in roots]
    if inherited:
        entries.append(inherited)
    return {"PYTHONPATH": os.pathsep.join(entries)}


def run(argv: list[str], *, cwd: Path, extra_env: dict[str, str] | None = None, timeout: int = 90):
    env = os.environ.copy()
    env.update(extra_env or {})
    return subprocess.run(argv, cwd=cwd, env=env, capture_output=True, text=True, timeout=timeout)


def assert_success(result: subprocess.CompletedProcess[str]) -> None:
    assert result.returncode == 0, (
        f"Command failed ({result.returncode}): {result.args}\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
    )


def ordered_panel(tmp_path: Path, names: tuple[str, ...] = ("Gata4", "Sox2", "Tbx5", "Pax6")) -> Path:
    panel = tmp_path / "ordered_panel.txt"
    panel.write_text("\n".join(names) + "\n")
    return panel


def make_adata(
    tmp_path: Path,
    *,
    name: str,
    genes: tuple[str, ...] = ("Gata4", "Sox2", "Tbx5", "Pax6"),
    n_obs: int = 6,
    sparse: bool = False,
    sparse_format: str = "csr",
    dtype: str = "float32",
    include_celltype: bool = True,
    spatial: np.ndarray | None = None,
    values: np.ndarray | None = None,
    uns: dict | None = None,
) -> Path:
    ad = require_anndata()
    matrix = values if values is not None else np.arange(n_obs * len(genes), dtype=dtype).reshape(n_obs, len(genes))
    matrix = np.asarray(matrix, dtype=dtype)
    if sparse:
        scipy_sparse = pytest.importorskip("scipy.sparse")
        if sparse_format == "csr":
            matrix = scipy_sparse.csr_matrix(matrix)
        elif sparse_format == "csc":
            matrix = scipy_sparse.csc_matrix(matrix)
        else:
            raise ValueError(f"Unsupported test sparse format: {sparse_format}")
    obs = {"celltype": ["type_a" if n % 2 == 0 else "type_b" for n in range(n_obs)]} if include_celltype else {}
    item = ad.AnnData(X=matrix, obs=obs)
    item.var_names = list(genes)
    if uns:
        item.uns.update(uns)
    if spatial is not None:
        item.obsm["spatial_3D"] = spatial
    output = tmp_path / f"{name}.h5ad"
    item.write_h5ad(output)
    return output


def make_large_sparse_sentinel(tmp_path: Path) -> Path:
    """Small on disk, but 1.6 GB if converted to a 100k x 2k float64 dense array."""
    ad = require_anndata()
    scipy_sparse = pytest.importorskip("scipy.sparse")
    matrix = scipy_sparse.csr_matrix((100_000, 2_000), dtype=np.float64)
    item = ad.AnnData(X=matrix)
    item.var_names = [f"g{index}" for index in range(item.n_vars)]
    output = tmp_path / "large_sparse_sentinel.h5ad"
    item.write_h5ad(output)
    return output
