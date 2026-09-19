"""Small, read-only HDF5 inspection helpers for AnnData files."""

from .inspector import inspect_h5ad, write_reports

__all__ = ["inspect_h5ad", "write_reports"]
