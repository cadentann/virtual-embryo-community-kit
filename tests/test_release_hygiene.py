from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FORBIDDEN_PARTS = {"00_sources", "06_audit", "07_handoff"}
RUNTIME_GENERATED_PARTS = {
    "__pycache__",
    ".pytest_cache",
    ".venv",
    "local_outputs",
    "artifacts",
    "reports",
    "generated",
    "official_contract_cache",
}
FORBIDDEN_SUFFIXES = {".h5ad", ".h5", ".hdf5", ".zarr", ".loom"}
FORBIDDEN_TEXT = ("/Users/", "/private/", "Downloads", "ghp_", "BEGIN OPENSSH PRIVATE KEY")


def test_allowlist_tree_has_no_private_or_generated_residue() -> None:
    for path in ROOT.rglob("*"):
        rel = path.relative_to(ROOT).as_posix()
        parts = path.relative_to(ROOT).parts
        # The documented install, synthetic quickstarts, and tests can create
        # these ignored local-only directories. They are skipped here; the
        # pre-ZIP packaging gate scans the clean release tree separately and
        # rejects them before shipment.
        if any(part in RUNTIME_GENERATED_PARTS or part.endswith(".egg-info") for part in parts):
            continue
        assert not path.is_symlink(), f"symlink is not allowed: {path}"
        assert not any(part in FORBIDDEN_PARTS or Path(part).suffix.lower() in FORBIDDEN_SUFFIXES for part in parts), rel
        if path.is_file():
            if path == Path(__file__):
                continue
            if path.suffix.lower() in {".md", ".json", ".py", ".toml", ".txt", ".ipynb", ".gitignore"}:
                text = path.read_text(errors="ignore")
                assert not any(token in text for token in FORBIDDEN_TEXT), rel


def test_notebooks_are_clean() -> None:
    for path in ROOT.rglob("*.ipynb"):
        data = json.loads(path.read_text())
        assert all(not cell.get("outputs") for cell in data.get("cells", [])), path
        assert all(cell.get("execution_count") is None for cell in data.get("cells", [])), path
