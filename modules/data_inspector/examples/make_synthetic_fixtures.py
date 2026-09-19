"""Create tiny deterministic AnnData-shaped HDF5 fixtures; no challenge data is used."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np


def _strings(group: h5py.Group, name: str, values: list[str]) -> h5py.Dataset:
    return group.create_dataset(name, data=np.asarray(values, dtype=h5py.string_dtype("utf-8")))


def _axis(group: h5py.File, name: str, index: list[str], extra: dict[str, list[str]] | None = None) -> None:
    axis = group.create_group(name)
    axis.attrs["_index"] = "_index"
    _strings(axis, "_index", index)
    for key, values in (extra or {}).items():
        _strings(axis, key, values)


def dense_fixture(path: Path) -> None:
    """A T1-like synthetic dense file: no spatial field, ordered fake genes, celltype."""
    with h5py.File(path, "w") as handle:
        x = handle.create_dataset("X", data=np.array([[0, 1, 2, 0], [3, 0, 0, 4], [0, 5, 0, 6]], dtype=np.float32))
        x.attrs["encoding-type"] = "array"
        _axis(handle, "obs", ["synthetic-cell-1", "synthetic-cell-2", "synthetic-cell-3"], {"celltype": ["type_a", "type_b", "type_a"]})
        _axis(handle, "var", ["SYN_GENE_A", "SYN_GENE_B", "SYN_GENE_C", "SYN_GENE_D"])
        handle.create_group("obsm")
        handle.create_group("layers")
        handle.create_group("uns")


def sparse_spatial_fixture(path: Path) -> None:
    """A T2/T3-like synthetic CSR file with local-frame 3D coordinates."""
    with h5py.File(path, "w") as handle:
        x = handle.create_group("X")
        x.attrs["encoding-type"] = "csr_matrix"
        x.attrs["shape"] = np.array([4, 3], dtype=np.int64)
        x.create_dataset("data", data=np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float32))
        x.create_dataset("indices", data=np.array([0, 2, 1, 2], dtype=np.int32))
        x.create_dataset("indptr", data=np.array([0, 2, 2, 3, 4], dtype=np.int32))
        _axis(handle, "obs", ["synthetic-cell-1", "synthetic-cell-2", "synthetic-cell-3", "synthetic-cell-4"], {"celltype": ["mesoderm", "mesoderm", "other", "other"]})
        _axis(handle, "var", ["SYN_PANEL_1", "SYN_PANEL_2", "SYN_PANEL_3"])
        obsm = handle.create_group("obsm")
        obsm.create_dataset("spatial_3D", data=np.array([[0, 0, 0], [2, 0, 1], [2, 4, 1], [4, 4, 3]], dtype=np.float32))
        handle.create_group("layers")
        handle.create_group("uns")


def _sparse_x(handle: h5py.File, n_obs: int, n_vars: int, *, entries_per_row: int = 3) -> None:
    """Write a tiny deterministic CSR matrix without constructing a dense array."""
    x = handle.create_group("X")
    x.attrs["encoding-type"] = "csr_matrix"
    x.attrs["shape"] = np.array([n_obs, n_vars], dtype=np.int64)
    positions = np.arange(n_obs * entries_per_row, dtype=np.int64)
    x.create_dataset("data", data=(positions + 1).astype(np.float32))
    x.create_dataset("indices", data=((positions * 37) % n_vars).astype(np.int32))
    x.create_dataset("indptr", data=np.arange(0, n_obs * entries_per_row + 1, entries_per_row, dtype=np.int32))


def extended_wide_sparse_fixture(path: Path) -> None:
    """A 2,048-feature sparse file with labels and deliberately no spatial field."""
    n_obs, n_vars = 8, 2_048
    with h5py.File(path, "w") as handle:
        _sparse_x(handle, n_obs, n_vars)
        _axis(handle, "obs", [f"wide-cell-{i}" for i in range(n_obs)], {"celltype": ["kind_a", "kind_b"] * (n_obs // 2)})
        _axis(handle, "var", [f"SYN_WIDE_{i:04d}" for i in range(n_vars)])
        handle.create_group("obsm")
        handle.create_group("layers")
        handle.create_group("uns")


def extended_spatial_500_fixture(path: Path) -> None:
    """A 500-feature sparse file with labels and a local 3D coordinate field."""
    n_obs, n_vars = 10, 500
    with h5py.File(path, "w") as handle:
        _sparse_x(handle, n_obs, n_vars)
        _axis(handle, "obs", [f"panel-cell-{i}" for i in range(n_obs)], {"celltype": ["kind_a"] * 5 + ["kind_b"] * 5})
        _axis(handle, "var", [f"SYN_PANEL_{i:03d}" for i in range(n_vars)])
        obsm = handle.create_group("obsm")
        coordinates = np.column_stack((np.arange(n_obs), np.arange(n_obs) * 2, np.arange(n_obs) * 3)).astype(np.float32)
        obsm.create_dataset("spatial_3D", data=coordinates)
        handle.create_group("layers")
        handle.create_group("uns")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("generated"))
    parser.add_argument("--extended", action="store_true", help="Also generate wider generated-only sparse fixtures for shape coverage")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    dense_fixture(args.output_dir / "synthetic_dense_t1_like.h5ad")
    sparse_spatial_fixture(args.output_dir / "synthetic_sparse_spatial_t2_like.h5ad")
    (args.output_dir / "synthetic_panel_reordered.txt").write_text("SYN_GENE_D\nSYN_GENE_C\nSYN_GENE_B\nSYN_GENE_A\n", encoding="utf-8")
    if args.extended:
        extended_wide_sparse_fixture(args.output_dir / "synthetic_sparse_2048_genes_no_spatial.h5ad")
        extended_spatial_500_fixture(args.output_dir / "synthetic_sparse_500_genes_spatial.h5ad")
    print(f"Wrote deterministic synthetic fixtures to {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
