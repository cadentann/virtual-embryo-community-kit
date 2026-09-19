"""Format-oriented helpers; these do not score a prediction or access held-out truth."""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import ssl
from typing import Any, Iterable
from urllib.request import Request, urlopen

import anndata as ad
import certifi
import numpy as np
import pandas as pd
from scipy import sparse


OFFICIAL_PANEL_INDEX_URL = "https://virtualembryo.ai/challenge/panels/index.json"
OFFICIAL_PANEL_BASE_URL = "https://virtualembryo.ai/challenge/panels/"
FETCH_USER_AGENT = "VEC-T1-Walkthrough/0.1 (+https://virtualembryo.ai/challenge)"
DEFAULT_SAMPLING_SEED = 0


class ContractError(ValueError):
    """Raised when board metadata or a gene panel cannot support a safe output."""


@dataclass
class ValidationReport:
    """Structural checks only; passing does not mean an upload will score well."""

    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.errors

    def require_ok(self) -> None:
        if self.errors:
            raise ContractError("Task 1 structural validation failed:\n- " + "\n- ".join(self.errors))


def _unique_or_error(names: Iterable[str], label: str) -> list[str]:
    values = [str(x) for x in names]
    duplicate = pd.Index(values)[pd.Index(values).duplicated()].unique().tolist()
    if duplicate:
        preview = ", ".join(map(str, duplicate[:5]))
        raise ContractError(f"{label} contains duplicate gene names: {preview}")
    return values


def fetch_official_contract(cache_dir: str | Path, board: str = "T1:val", timeout: int = 30) -> tuple[dict[str, Any], list[str]]:
    """Fetch and cache current official metadata/panel plus provenance.

    Network access is opt-in. The cached files are a record of a fetch, not a replacement
    for checking the current official site before submission.
    """
    cache = Path(cache_dir)
    cache.mkdir(parents=True, exist_ok=True)
    # Keep TLS verification enabled while using certifi's portable CA bundle. This
    # avoids relying on a macOS Python installation's sometimes-unconfigured store.
    verified_context = ssl.create_default_context(cafile=certifi.where())
    with urlopen(Request(OFFICIAL_PANEL_INDEX_URL, headers={"User-Agent": FETCH_USER_AGENT}), timeout=timeout, context=verified_context) as response:
        index_bytes = response.read()
    index = json.loads(index_bytes.decode("utf-8"))
    if board not in index:
        raise ContractError(f"Board {board!r} is absent from official index. Available: {', '.join(index)}")
    metadata = index[board]
    if metadata.get("task") != "T1" or board != "T1:val":
        raise ContractError("This walkthrough intentionally supports only the current T1:val contract.")
    filename = metadata.get("genes_file")
    if not filename or Path(filename).name != filename:
        raise ContractError("Official metadata supplied an unsafe or absent genes_file.")
    panel_url = OFFICIAL_PANEL_BASE_URL + filename
    with urlopen(Request(panel_url, headers={"User-Agent": FETCH_USER_AGENT}), timeout=timeout, context=verified_context) as response:
        panel_bytes = response.read()
    panel = [line.strip() for line in panel_bytes.decode("utf-8").splitlines() if line.strip()]
    _validate_board_contract(metadata, panel)
    # The official prefix is for the normalized newline-delimited gene sequence (no
    # terminal newline), so a harmless transport newline does not create a false mismatch.
    expected_prefix = metadata.get("genes_sha256")
    actual_hash = sha256("\n".join(panel).encode("utf-8")).hexdigest()
    if expected_prefix and not actual_hash.startswith(str(expected_prefix)):
        raise ContractError(
            "Downloaded panel SHA-256 does not match the prefix recorded in official metadata; "
            "refusing to cache an ambiguous contract."
        )
    index_path, panel_path = cache / "index.json", cache / filename
    index_path.write_bytes(index_bytes)
    panel_path.write_bytes(panel_bytes)
    provenance = {
        "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
        "board": board,
        "index_url": OFFICIAL_PANEL_INDEX_URL,
        "panel_url": panel_url,
        "index_sha256": sha256(index_bytes).hexdigest(),
        "panel_normalized_sha256": actual_hash,
        "panel_download_sha256": sha256(panel_bytes).hexdigest(),
        "official_panel_sha256_prefix": metadata.get("genes_sha256"),
        "note": "Local cache of a live official fetch; verify current requirements before upload.",
    }
    (cache / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    return metadata, panel


def read_cached_contract(cache_dir: str | Path, board: str = "T1:val") -> tuple[dict[str, Any], list[str]]:
    """Read a previously fetched contract without contacting the network."""
    cache = Path(cache_dir)
    index_path = cache / "index.json"
    provenance_path = cache / "provenance.json"
    try:
        index_bytes = index_path.read_bytes()
        index = json.loads(index_bytes.decode("utf-8"))
        provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ContractError(f"Could not read a complete, valid contract cache at {cache}: {exc}") from exc
    if board not in index:
        raise ContractError(f"Board {board!r} is absent from {index_path}.")
    metadata = index[board]
    filename = metadata.get("genes_file")
    if not filename or Path(filename).name != filename:
        raise ContractError("Cached metadata has an unsafe or absent genes_file.")
    panel_path = cache / filename
    try:
        panel_bytes = panel_path.read_bytes()
        panel = [line.strip() for line in panel_bytes.decode("utf-8").splitlines() if line.strip()]
    except (OSError, UnicodeDecodeError) as exc:
        raise ContractError(f"Could not read cached panel {panel_path}: {exc}") from exc
    _validate_board_contract(metadata, panel)
    normalized_hash = sha256("\n".join(panel).encode("utf-8")).hexdigest()
    expected_prefix = metadata.get("genes_sha256")
    if not expected_prefix or not normalized_hash.startswith(str(expected_prefix)):
        raise ContractError("Cached panel hash does not match the official prefix recorded in index.json.")
    checks = {
        "index_sha256": sha256(index_bytes).hexdigest(),
        "panel_download_sha256": sha256(panel_bytes).hexdigest(),
        "panel_normalized_sha256": normalized_hash,
    }
    for key, actual in checks.items():
        if provenance.get(key) != actual:
            raise ContractError(f"Cached contract integrity check failed for {key}; re-fetch the official contract.")
    if provenance.get("board") != board:
        raise ContractError("Cached provenance board does not match the requested board.")
    return metadata, panel


def _validate_board_contract(metadata: dict[str, Any], panel: Iterable[str]) -> None:
    genes = _unique_or_error(panel, "Official panel")
    if metadata.get("task") != "T1":
        raise ContractError("Expected Task 1 metadata.")
    if metadata.get("needs_coords") is not False:
        raise ContractError("Task 1 output must not require spatial coordinates.")
    expected = metadata.get("n_genes")
    if not isinstance(expected, int) or len(genes) != expected:
        raise ContractError(f"Panel length {len(genes)} disagrees with metadata n_genes={expected!r}.")
    if not isinstance(metadata.get("min_cells"), int) or not isinstance(metadata.get("max_cells"), int):
        raise ContractError("Board metadata has no usable cell-count bounds.")
    if metadata["min_cells"] < 1 or metadata["min_cells"] > metadata["max_cells"]:
        raise ContractError("Board metadata has inconsistent cell-count bounds.")


def make_synthetic_training_inputs(panel: Iterable[str] | None = None) -> tuple[ad.AnnData, ad.AnnData, list[str], dict[str, Any]]:
    """Return deterministic E8.5/E9.5-like *synthetic* inputs and a tiny fake contract.

    The fake board is deliberately named SYNTHETIC:T1:val and must never be used for upload.
    """
    genes = ["GeneA", "GeneB", "GeneC", "GeneD"] if panel is None else _unique_or_error(panel, "Synthetic panel")
    if not genes:
        raise ContractError("Synthetic panel must contain at least one gene.")
    # A and B occur in both stages. C is deliberately E9.5-only, exercising the
    # conservative unchanged-late-cell behavior used when no earlier mean exists.
    early_values = np.array([
        [1, 0, 2, 0], [3, 0, 0, 2],  # A: mean [2, 0, 1, 1]
        [0, 1, 2, 1], [2, 3, 0, 1],  # B: mean [1, 2, 1, 1]
    ], dtype=np.float32)
    late_values = np.array([
        [4, 1, 3, 1], [2, 2, 1, 3],  # A: delta [1, 1.5, 1, 1]
        [1, 4, 2, 2], [3, 2, 2, 4],  # B: delta [1, 1, 1, 2]
        [5, 0, 1, 0],                 # C: late-only, unchanged
    ], dtype=np.float32)
    if len(genes) != 4:
        # Deterministic values derived from panel position, not biological data.
        position = np.arange(len(genes), dtype=np.float32)
        early_values = np.vstack([(position % 3) + 1, (position % 4) + 1, (position % 2) + 1, (position % 5) + 1])
        late_values = np.vstack([early_values[0] + 2, early_values[1] + 1, early_values[2] + 3, early_values[3] + 1, (position % 4) + 1])
    early_x = sparse.csr_matrix(early_values)
    late_x = sparse.csr_matrix(late_values)
    early = ad.AnnData(
        early_x,
        obs=pd.DataFrame({"celltype": ["A", "A", "B", "B"]}, index=["E8.5_A0", "E8.5_A1", "E8.5_B0", "E8.5_B1"]),
        var=pd.DataFrame(index=genes),
    )
    late = ad.AnnData(
        late_x,
        obs=pd.DataFrame({"celltype": ["A", "A", "B", "B", "C"]}, index=["E9.5_A0", "E9.5_A1", "E9.5_B0", "E9.5_B1", "E9.5_C0"]),
        var=pd.DataFrame(index=genes),
    )
    contract = {"key": "SYNTHETIC:T1:val", "task": "T1", "n_genes": len(genes), "min_cells": 3, "max_cells": 5, "needs_coords": False}
    return early, late, genes, contract


def _require_same_ordered_genes(early: ad.AnnData, late: ad.AnnData, panel: Iterable[str]) -> list[str]:
    genes = _unique_or_error(panel, "Requested panel")
    for label, data in (("E8.5 input", early), ("E9.5 input", late)):
        observed = _unique_or_error(data.var_names, label)
        missing = sorted(set(genes) - set(observed))
        if missing:
            raise ContractError(f"{label} is missing {len(missing)} panel genes (e.g. {', '.join(missing[:5])}).")
    return genes


def _celltype_vector(data: ad.AnnData, label: str) -> np.ndarray:
    if "celltype" not in data.obs:
        raise ContractError(f"{label} requires obs['celltype'] for the reference per-cell-type baseline.")
    labels = data.obs["celltype"]
    if labels.isna().any():
        raise ContractError(f"{label} has missing values in obs['celltype'].")
    return labels.astype(str).to_numpy()


def _require_chunk_rows(chunk_rows: int) -> None:
    if not isinstance(chunk_rows, int) or isinstance(chunk_rows, bool) or chunk_rows < 1:
        raise ContractError("chunk_rows must be a positive integer.")


def _read_physical_rows_in_panel_order(
    x: Any,
    physical_rows: np.ndarray,
    positions: np.ndarray,
) -> Any:
    """Read increasing physical rows, then select panel columns in memory.

    ``h5py`` dense datasets reject non-increasing fancy row indexes. AnnData's
    backed sparse datasets are also most predictable with monotonic physical row
    access. Selecting columns after the bounded row read keeps a shuffled panel
    valid without relying on multi-axis HDF5 fancy indexing.
    """
    if len(physical_rows) and np.any(np.diff(physical_rows) <= 0):
        raise ContractError("Internal error: backed reads require strictly increasing physical row indexes.")
    source_rows = x[physical_rows, :]
    return source_rows[:, positions]


def _mean_rows_in_panel_order(
    data: ad.AnnData,
    row_indices: np.ndarray,
    positions: np.ndarray,
    *,
    chunk_rows: int,
) -> np.ndarray:
    """Calculate a panel mean from bounded, increasing physical row reads."""
    _require_chunk_rows(chunk_rows)
    rows = np.asarray(row_indices, dtype=np.int64)
    if rows.ndim != 1 or len(rows) == 0:
        raise ContractError("Cannot calculate a pseudobulk mean from an empty row selection.")
    ordered_rows = np.sort(rows, kind="stable")
    total = np.zeros(len(positions), dtype=np.float64)
    for start in range(0, len(ordered_rows), chunk_rows):
        source = _read_physical_rows_in_panel_order(
            data.X, ordered_rows[start : start + chunk_rows], positions,
        )
        total += np.asarray(source.sum(axis=0), dtype=np.float64).ravel()
    return total / len(ordered_rows)


def _uniform_sample_rows(
    n_available: int,
    n_cells: int,
    seed: int,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Draw a reproducible uniform subset in the generator's logical output order."""
    if not isinstance(seed, (int, np.integer)) or isinstance(seed, bool):
        raise ContractError("sampling_seed must be an integer.")
    try:
        selected = np.random.default_rng(int(seed)).choice(
            n_available, size=n_cells, replace=False,
        ).astype(np.int64, copy=False)
    except ValueError as exc:
        raise ContractError(f"Could not sample {n_cells} of {n_available} late-stage cells: {exc}") from exc
    record = {
        "method": "uniform_without_replacement",
        "random_generator": "numpy.random.default_rng",
        "seed": int(seed),
        "n_available_late_cells": int(n_available),
        "n_selected_cells": int(n_cells),
        "selected_input_row_indices": selected.tolist(),
        "output_order": "numpy Generator.choice draw order, restored after physical-order reads",
    }
    record["selected_input_row_indices_sha256"] = sha256(
        json.dumps(record["selected_input_row_indices"], separators=(",", ":")).encode("utf-8"),
    ).hexdigest()
    return selected, record


def _sparse_safe_shifted_rows(
    late: ad.AnnData,
    row_indices: np.ndarray,
    positions: np.ndarray,
    late_types: np.ndarray,
    shifts: dict[str, np.ndarray],
    chunk_rows: int = 128,
) -> sparse.csr_matrix:
    """Apply shifts in bounded blocks and restore the requested logical row order."""
    _require_chunk_rows(chunk_rows)
    requested_rows = np.asarray(row_indices, dtype=np.int64)
    if requested_rows.ndim != 1:
        raise ContractError("Internal error: row indices must be one-dimensional.")
    physical_order = np.argsort(requested_rows, kind="stable")
    physical_rows = requested_rows[physical_order]
    restore_order = np.empty(len(requested_rows), dtype=np.int64)
    restore_order[physical_order] = np.arange(len(requested_rows))
    blocks: list[sparse.csr_matrix] = []
    for start in range(0, len(physical_rows), chunk_rows):
        rows = physical_rows[start : start + chunk_rows]
        source = _read_physical_rows_in_panel_order(late.X, rows, positions)
        # Only this bounded output block becomes dense. It is necessary when a shift
        # turns implicit zeros into nonzero values; the full source is never densified.
        values = source.toarray() if sparse.issparse(source) else np.asarray(source).copy()
        row_shifts = np.vstack([shifts.get(celltype, np.zeros(len(positions), dtype=np.float32)) for celltype in late_types[rows]])
        shifted = np.maximum(values.astype(np.float32, copy=False) + row_shifts, 0, dtype=np.float32)
        blocks.append(sparse.csr_matrix(shifted, dtype=np.float32))
    in_physical_order = sparse.vstack(blocks, format="csr", dtype=np.float32)
    return in_physical_order[restore_order, :]


def _expression_problem(x: Any) -> str | None:
    """Return a concise structural reason why ``x`` cannot be Task-1 expression."""
    if x is None:
        return "is missing."
    try:
        shape = x.shape
    except AttributeError:
        return "has no shape."
    try:
        if len(shape) != 2:
            return "must be two-dimensional."
    except TypeError:
        return "must be two-dimensional."
    try:
        dtype = np.dtype(x.dtype)
    except (AttributeError, TypeError):
        return "must have a real numeric dtype."
    if np.issubdtype(dtype, np.bool_):
        return "boolean dtype is unsupported: a presence mask is not a real-valued expression matrix."
    if np.issubdtype(dtype, np.complexfloating):
        return "complex dtype is unsupported: Task 1 expression must be real-valued."
    if not np.issubdtype(dtype, np.number):
        return f"must have a real numeric dtype, not {dtype}."
    return None


def _require_supported_expression(x: Any, label: str) -> None:
    problem = _expression_problem(x)
    if problem is not None:
        raise ContractError(f"{label} {problem}")


def _require_input_expression(data: ad.AnnData, label: str) -> None:
    """Turn absent backed ``.X`` into the same concise builder error as ``None``."""
    try:
        x = data.X
    except (AttributeError, KeyError) as exc:
        raise ContractError(f"{label} is missing.") from exc
    _require_supported_expression(x, label)


def build_pseudobulk_shift_prediction(
    early: ad.AnnData,
    late: ad.AnnData,
    panel: Iterable[str],
    n_cells: int,
    contract: dict[str, Any],
    *,
    sampling_seed: int = DEFAULT_SAMPLING_SEED,
    sampling_report: dict[str, Any] | None = None,
    chunk_rows: int = 128,
) -> ad.AnnData:
    """Build a minimal output using the reference per-cell-type pseudobulk shift.

    For cell types present at E8.5 and E9.5, calculate the type-specific mean shift
    and add it to every E9.5 cell of that type, retaining within-type differences
    except where non-negative clipping is required. E9.5-only types are copied
    unchanged: guessing a missing historical shift would be an unsupported choice.
    The output has no celltype metadata because that label is not needed for a minimal
    Task 1 prediction. Late-stage cells are selected uniformly without replacement
    using a fixed seed by default; callers can capture the policy in
    ``sampling_report`` and serialize it as a sidecar. The prediction itself gets
    no sampling metadata. This is a pedagogical reference workflow, not a scorer.
    """
    panel_list = _require_same_ordered_genes(early, late, panel)
    _validate_board_contract(contract, panel_list)
    if not isinstance(n_cells, int) or not (contract["min_cells"] <= n_cells <= contract["max_cells"]):
        raise ContractError(f"n_cells must be an integer within board bounds [{contract['min_cells']}, {contract['max_cells']}].")
    if n_cells > late.n_obs:
        raise ContractError(
            f"Cannot select {n_cells} output cells from an E9.5 input containing only {late.n_obs}; "
            "choose an in-range count no larger than the available late-stage cells."
        )
    # Validate before the mean and output paths cast to float64/float32. In
    # particular, casting complex values would silently discard their imaginary part.
    _require_input_expression(early, "E8.5 input .X")
    _require_input_expression(late, "E9.5 input .X")
    _require_chunk_rows(chunk_rows)
    early_types = _celltype_vector(early, "E8.5 input")
    late_types = _celltype_vector(late, "E9.5 input")
    early_positions = pd.Index(early.var_names).get_indexer(panel_list)
    late_positions = pd.Index(late.var_names).get_indexer(panel_list)
    shifts: dict[str, np.ndarray] = {}
    for celltype in sorted(set(early_types) & set(late_types)):
        early_mean = _mean_rows_in_panel_order(
            early, np.flatnonzero(early_types == celltype), early_positions, chunk_rows=chunk_rows,
        )
        late_mean = _mean_rows_in_panel_order(
            late, np.flatnonzero(late_types == celltype), late_positions, chunk_rows=chunk_rows,
        )
        shift = (late_mean - early_mean).astype(np.float32)
        if not np.isfinite(shift).all():
            raise ContractError(f"Pseudobulk shift for celltype {celltype!r} contains non-finite values.")
        shifts[celltype] = shift
    selected, record = _uniform_sample_rows(late.n_obs, n_cells, sampling_seed)
    if sampling_report is not None:
        sampling_report.clear()
        sampling_report.update(record)
    x = _sparse_safe_shifted_rows(
        late, selected, late_positions, late_types, shifts, chunk_rows=chunk_rows,
    )
    prediction = ad.AnnData(
        X=x,
        obs=pd.DataFrame(index=[f"prediction_{i:05d}" for i in range(n_cells)]),
        var=pd.DataFrame(index=pd.Index(panel_list, name=None)),
    )
    return prediction


def _x_has_invalid_values(x: Any, chunk_rows: int = 2048) -> tuple[bool, bool]:
    """Return (has_nonfinite, has_negative) without densifying sparse matrices."""
    _require_supported_expression(x, "Prediction .X")
    if sparse.issparse(x):
        values = x.data
        return bool(not np.isfinite(values).all()), bool((values < 0).any())
    for start in range(0, x.shape[0], chunk_rows):
        block = x[start : start + chunk_rows]
        # AnnData backed sparse datasets yield a scipy sparse block on slice.
        values = block.data if sparse.issparse(block) else np.asarray(block)
        if not np.isfinite(values).all():
            return True, bool((values < 0).any())
        if (values < 0).any():
            return False, True
    return False, False


def validate_task1_prediction(prediction: ad.AnnData, panel: Iterable[str], contract: dict[str, Any]) -> ValidationReport:
    """Check the published Task-1 structural contract; never scores or mutates data."""
    report = ValidationReport()
    try:
        panel_list = _unique_or_error(panel, "Official panel")
        _validate_board_contract(contract, panel_list)
    except ContractError as exc:
        report.errors.append(str(exc))
        return report
    observed = [str(x) for x in prediction.var_names]
    duplicates = pd.Index(observed)[pd.Index(observed).duplicated()].unique().tolist()
    if duplicates:
        report.errors.append(f"Prediction var_names contains duplicate genes: {', '.join(map(str, duplicates[:5]))}.")
    if observed != panel_list:
        observed_set, panel_set = set(observed), set(panel_list)
        missing, extra = sorted(panel_set - observed_set), sorted(observed_set - panel_set)
        if not missing and not extra:
            report.errors.append("Prediction has the right gene set but wrong var_names order.")
        else:
            if missing:
                report.errors.append(f"Prediction is missing {len(missing)} required genes (e.g. {', '.join(missing[:5])}).")
            if extra:
                report.errors.append(f"Prediction has {len(extra)} unexpected genes (e.g. {', '.join(extra[:5])}).")
    if prediction.n_vars != contract["n_genes"]:
        report.errors.append(f"Expected {contract['n_genes']} variables, got {prediction.n_vars}.")
    if not (contract["min_cells"] <= prediction.n_obs <= contract["max_cells"]):
        report.errors.append(f"Expected n_obs in [{contract['min_cells']}, {contract['max_cells']}], got {prediction.n_obs}.")
    try:
        x = prediction.X
    except (AttributeError, KeyError):
        x = None
    expression_problem = _expression_problem(x)
    if expression_problem is not None:
        report.errors.append(f"Prediction .X {expression_problem}")
    else:
        nonfinite, negative = _x_has_invalid_values(x)
        if nonfinite:
            report.errors.append("Prediction .X contains NaN or infinite values.")
        if negative:
            report.errors.append("Prediction .X contains negative values.")
    if "spatial_3D" in prediction.obsm:
        report.errors.append("This companion's minimal-output policy rejects spatial_3D coordinates; confirm current server behavior separately.")
    if "celltype" in prediction.obs:
        report.warnings.append("celltype is training-only metadata; omit it from a minimal Task 1 output.")
    if expression_problem is None and x.dtype != np.float32:
        report.warnings.append(f".X dtype is {x.dtype}; float32 is the bounded recommended output dtype.")
    report.warnings.append(
        "Structural checks cannot establish the required challenge normalization; finite, non-negative raw counts "
        "or another scale can pass locally and still be unsuitable for scoring."
    )
    return report
