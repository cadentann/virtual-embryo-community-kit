from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import ssl
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
import anndata as ad
import h5py

from vec_t1_walkthrough.workflow import (
    ContractError,
    build_pseudobulk_shift_prediction,
    make_synthetic_training_inputs,
    read_cached_contract,
    validate_task1_prediction,
)
import vec_t1_walkthrough.workflow as workflow
from vec_t1_walkthrough.cli import main as cli_main


def prepared():
    early, late, panel, contract = make_synthetic_training_inputs()
    return build_pseudobulk_shift_prediction(early, late, panel, 3, contract), panel, contract


def test_synthetic_output_is_minimal_and_valid():
    prediction, panel, contract = prepared()
    report = validate_task1_prediction(prediction, panel, contract)
    assert report.ok
    assert prediction.X.dtype == np.float32
    assert sparse.isspmatrix_csr(prediction.X)
    assert list(prediction.var_names) == panel
    assert "celltype" not in prediction.obs
    assert "spatial_3D" not in prediction.obsm


def test_per_celltype_shift_has_hand_checked_values_and_retains_heterogeneity():
    early, late, panel, contract = make_synthetic_training_inputs()
    prediction = build_pseudobulk_shift_prediction(early, late, panel, 5, contract)
    actual = prediction.X.toarray()
    expected = np.array([
        [5, 2.5, 4, 2],  # E9.5 A0 + (mean_A_late - mean_A_early) = +[1, 1.5, 1, 1]
        [3, 3.5, 2, 4],
        [2, 5, 3, 4],    # E9.5 B0 + [1, 1, 1, 2]
        [4, 3, 3, 6],
        [5, 0, 1, 0],    # E9.5-only C is deliberately unchanged
    ], dtype=np.float32)
    selected = np.random.default_rng(0).choice(late.n_obs, size=5, replace=False)
    np.testing.assert_allclose(actual, expected[selected])
    late_values = late.X.toarray()
    # Within-type A variation remains exactly after the same additive shift.
    a_positions = np.flatnonzero(late.obs["celltype"].to_numpy()[selected] == "A")
    np.testing.assert_allclose(
        actual[a_positions[0]] - actual[a_positions[1]],
        late_values[selected[a_positions[0]]] - late_values[selected[a_positions[1]]],
    )
    # Output strips training-only labels even though the baseline uses them internally.
    assert "celltype" not in prediction.obs


def test_baseline_requires_celltype_in_each_reference_input():
    early, late, panel, contract = make_synthetic_training_inputs()
    del late.obs["celltype"]
    with pytest.raises(ContractError, match=r"requires obs\['celltype'\]"):
        build_pseudobulk_shift_prediction(early, late, panel, 3, contract)


@pytest.mark.parametrize(
    ("mutation", "fragment"),
    [
        (lambda a: a[:, [1, 0, 2, 3]].copy(), "wrong var_names order"),
        (lambda a: _rename(a, ["GeneA", "GeneB", "GeneC", "Extra"]), "missing 1 required"),
        (lambda a: _rename(a, ["GeneA", "GeneB", "GeneC", "GeneC"]), "duplicate"),
        (lambda a: _nan(a), "NaN or infinite"),
        (lambda a: _negative(a), "negative"),
    ],
)
def test_adversarial_outputs_are_rejected(mutation, fragment):
    prediction, panel, contract = prepared()
    report = validate_task1_prediction(mutation(prediction), panel, contract)
    assert not report.ok
    assert any(fragment in message for message in report.errors)


def _rename(adata, names):
    result = adata.copy()
    result.var_names = names
    return result


def _nan(adata):
    result = adata.copy()
    result.X.data[0] = np.nan
    return result


def _negative(adata):
    result = adata.copy()
    result.X.data[0] = -1
    return result


def test_extra_gene_is_reported_alongside_missing_gene():
    prediction, panel, contract = prepared()
    report = validate_task1_prediction(_rename(prediction, ["GeneA", "GeneB", "GeneC", "Extra"]), panel, contract)
    assert any("unexpected genes" in message for message in report.errors)


def test_dense_matrix_is_checked_without_changing_it():
    prediction, panel, contract = prepared()
    prediction = prediction.copy()
    prediction.X = prediction.X.toarray()
    before = prediction.X.copy()
    report = validate_task1_prediction(prediction, panel, contract)
    assert report.ok
    np.testing.assert_array_equal(prediction.X, before)


def test_backed_h5ad_is_checked_without_forcing_full_matrix(tmp_path):
    prediction, panel, contract = prepared()
    path = tmp_path / "prediction.h5ad"
    prediction.write_h5ad(path)
    backed = ad.read_h5ad(path, backed="r")
    report = validate_task1_prediction(backed, panel, contract)
    assert report.ok


@pytest.mark.parametrize(
    ("name", "matrix", "fragment"),
    [
        ("complex", np.full((3, 4), 1 + 2j, dtype=np.complex64), "complex dtype is unsupported"),
        ("string", np.full((3, 4), "x", dtype="U1"), "real numeric dtype"),
        ("object", np.full((3, 4), "x", dtype=object), "real numeric dtype"),
        ("bool", np.ones((3, 4), dtype=bool), "boolean dtype is unsupported"),
    ],
)
def test_validator_rejects_unsupported_real_anndata_expression_dtypes(name, matrix, fragment):
    prediction, panel, contract = prepared()
    prediction.X = matrix
    report = validate_task1_prediction(prediction, panel, contract)
    assert not report.ok, name
    assert any(fragment in error for error in report.errors)


def test_validator_rejects_missing_x_without_throwing():
    prediction, panel, contract = prepared()
    prediction.X = None
    report = validate_task1_prediction(prediction, panel, contract)
    assert not report.ok
    assert report.errors == ["Prediction .X is missing."]


def test_validator_rejects_non_two_dimensional_proxy_without_throwing():
    _, panel, contract = prepared()

    class OneDimensionalPrediction:
        X = np.array([1.0, 2.0, 3.0])
        var_names = pd.Index(panel)
        n_vars = len(panel)
        n_obs = 3
        obsm = {}
        obs = pd.DataFrame(index=["one", "two", "three"])

    report = validate_task1_prediction(OneDimensionalPrediction(), panel, contract)
    assert not report.ok
    assert report.errors == ["Prediction .X must be two-dimensional."]


def test_backed_complex_h5ad_is_rejected_without_loading_or_casting_full_matrix(tmp_path):
    prediction, panel, contract = prepared()
    prediction.X = prediction.X.astype(np.complex64)
    prediction.X.data[0] = 1 + 2j
    path = tmp_path / "complex_prediction.h5ad"
    prediction.write_h5ad(path)
    backed = ad.read_h5ad(path, backed="r")
    try:
        report = validate_task1_prediction(backed, panel, contract)
    finally:
        backed.file.close()
    assert not report.ok
    assert any("complex dtype is unsupported" in error for error in report.errors)


def test_backed_string_h5ad_is_rejected_without_throwing(tmp_path):
    prediction, panel, contract = prepared()
    prediction.X = np.full(prediction.shape, "x", dtype="U1")
    path = tmp_path / "string_prediction.h5ad"
    prediction.write_h5ad(path)
    backed = ad.read_h5ad(path, backed="r")
    try:
        report = validate_task1_prediction(backed, panel, contract)
    finally:
        backed.file.close()
    assert not report.ok
    assert any("real numeric dtype" in error for error in report.errors)


@pytest.mark.parametrize("stage", ["early", "late"])
def test_builder_rejects_complex_input_before_float_cast(stage):
    early, late, panel, contract = make_synthetic_training_inputs()
    data = early if stage == "early" else late
    data.X = data.X.astype(np.complex64)
    data.X.data[0] = 1 + 2j
    with pytest.raises(ContractError, match="complex dtype is unsupported"):
        build_pseudobulk_shift_prediction(early, late, panel, 3, contract)


def test_builder_rejects_boolean_input_under_explicit_policy():
    early, late, panel, contract = make_synthetic_training_inputs()
    late.X = np.ones(late.shape, dtype=bool)
    with pytest.raises(ContractError, match="boolean dtype is unsupported"):
        build_pseudobulk_shift_prediction(early, late, panel, 3, contract)


def test_builder_rejects_missing_input_x_before_working():
    early, late, panel, contract = make_synthetic_training_inputs()
    early.X = None
    with pytest.raises(ContractError, match=r"E8.5 input \.X is missing"):
        build_pseudobulk_shift_prediction(early, late, panel, 3, contract)


def _write_verified_cache(cache: Path, panel: list[str], min_cells: int, max_cells: int) -> None:
    cache.mkdir()
    panel_bytes = ("\n".join(panel) + "\n").encode("utf-8")
    normalized_hash = sha256("\n".join(panel).encode("utf-8")).hexdigest()
    index = {"T1:val": {
        "key": "T1:val", "task": "T1", "n_genes": len(panel), "genes_file": "panel.txt",
        "genes_sha256": normalized_hash[:16], "needs_coords": False,
        "min_cells": min_cells, "max_cells": max_cells,
    }}
    index_bytes = json.dumps(index).encode("utf-8")
    (cache / "index.json").write_bytes(index_bytes)
    (cache / "panel.txt").write_bytes(panel_bytes)
    (cache / "provenance.json").write_text(json.dumps({
        "board": "T1:val",
        "index_sha256": sha256(index_bytes).hexdigest(),
        "panel_download_sha256": sha256(panel_bytes).hexdigest(),
        "panel_normalized_sha256": normalized_hash,
    }), encoding="utf-8")


def test_validate_cli_rejects_complex_h5ad_cleanly_and_closes_handle(tmp_path, capsys):
    prediction, panel, contract = prepared()
    prediction.X = prediction.X.astype(np.complex64)
    prediction.X.data[0] = 1 + 2j
    prediction_path = tmp_path / "complex_prediction.h5ad"
    prediction.write_h5ad(prediction_path)
    cache = tmp_path / "cache"
    _write_verified_cache(cache, panel, contract["min_cells"], contract["max_cells"])

    assert cli_main(["validate", str(prediction_path), "--cache", str(cache)]) == 2
    captured = capsys.readouterr()
    assert "ERROR: Prediction .X complex dtype is unsupported" in captured.out
    assert "Traceback" not in captured.out + captured.err
    with h5py.File(prediction_path, "r+") as handle:
        assert "X" in handle


def test_validate_cli_rejects_missing_x_h5ad_cleanly_and_closes_handle(tmp_path, capsys):
    prediction, panel, contract = prepared()
    prediction.X = None
    prediction_path = tmp_path / "missing_x_prediction.h5ad"
    prediction.write_h5ad(prediction_path)
    cache = tmp_path / "cache"
    _write_verified_cache(cache, panel, contract["min_cells"], contract["max_cells"])

    assert cli_main(["validate", str(prediction_path), "--cache", str(cache)]) == 2
    captured = capsys.readouterr()
    assert "ERROR: Prediction .X is missing." in captured.out
    assert "Traceback" not in captured.out + captured.err
    with h5py.File(prediction_path, "r+") as handle:
        assert "X" not in handle


def test_build_cli_rejects_complex_input_before_writing_or_leaking_handles(tmp_path, capsys):
    early, late, panel, contract = make_synthetic_training_inputs()
    early.X = early.X.astype(np.complex64)
    early.X.data[0] = 1 + 2j
    early_path, late_path = tmp_path / "early_complex.h5ad", tmp_path / "late_valid.h5ad"
    early.write_h5ad(early_path)
    late.write_h5ad(late_path)
    cache = tmp_path / "cache"
    _write_verified_cache(cache, panel, contract["min_cells"], contract["max_cells"])
    output = tmp_path / "must_not_exist.h5ad"

    assert cli_main([
        "build", "--early", str(early_path), "--late", str(late_path),
        "--cache", str(cache), "--n-cells", "3", "--output", str(output),
    ]) == 2
    captured = capsys.readouterr()
    assert "ERROR: E8.5 input .X complex dtype is unsupported" in captured.err
    assert "Traceback" not in captured.out + captured.err
    assert not output.exists()
    assert not Path(str(output) + ".selection.json").exists()
    for source in (early_path, late_path):
        with h5py.File(source, "r+") as handle:
            assert "X" in handle


def test_builder_rejects_missing_input_gene():
    early, late, panel, contract = make_synthetic_training_inputs()
    late = late[:, :3].copy()
    with pytest.raises(ContractError, match="missing 1 panel genes"):
        build_pseudobulk_shift_prediction(early, late, panel, 3, contract)


def test_cell_bounds_are_read_from_contract_not_hardcoded():
    early, late, panel, contract = make_synthetic_training_inputs()
    with pytest.raises(ContractError, match=r"within board bounds \[3, 5\]"):
        build_pseudobulk_shift_prediction(early, late, panel, 2, contract)


def test_builder_rejects_more_cells_than_late_input_contains():
    early, late, panel, contract = make_synthetic_training_inputs()
    contract = {**contract, "max_cells": 10}
    with pytest.raises(ContractError, match="only 5"):
        build_pseudobulk_shift_prediction(early, late, panel, 6, contract)


def test_validation_does_not_mutate_source():
    prediction, panel, contract = prepared()
    before = prediction.X.copy()
    validate_task1_prediction(prediction, panel, contract)
    assert (prediction.X != before).nnz == 0


def test_live_contract_fetch_is_mocked_and_records_normalized_hash(monkeypatch, tmp_path):
    panel_bytes = b"GeneA\nGeneB\n"  # terminal newline differs from normalized hash on purpose
    normalized = sha256(b"GeneA\nGeneB").hexdigest()
    index = {"T1:val": {"key": "T1:val", "task": "T1", "n_genes": 2, "genes_file": "panel.txt", "genes_sha256": normalized[:16], "needs_coords": False, "min_cells": 1, "max_cells": 2}}

    class Response:
        def __init__(self, payload): self.payload = payload
        def read(self): return self.payload
        def __enter__(self): return self
        def __exit__(self, *_): return False

    payloads, contexts, requests = [json.dumps(index).encode(), panel_bytes], [], []
    def fake_urlopen(*args, **kwargs):
        contexts.append(kwargs["context"])
        requests.append(args[0])
        return Response(payloads.pop(0))
    monkeypatch.setattr(workflow, "urlopen", fake_urlopen)
    metadata, panel = workflow.fetch_official_contract(tmp_path)
    provenance = json.loads((tmp_path / "provenance.json").read_text())
    assert metadata["key"] == "T1:val"
    assert panel == ["GeneA", "GeneB"]
    assert provenance["panel_normalized_sha256"] == normalized
    assert len(contexts) == 2
    assert contexts[0] is contexts[1]
    assert isinstance(contexts[0], ssl.SSLContext)
    assert contexts[0].verify_mode == ssl.CERT_REQUIRED
    assert contexts[0].check_hostname is True
    assert all(request.get_header("User-agent") == workflow.FETCH_USER_AGENT for request in requests)
    cached_metadata, cached_panel = read_cached_contract(tmp_path)
    assert cached_metadata["key"] == "T1:val"
    assert cached_panel == panel

    (tmp_path / "panel.txt").write_text("OtherA\nOtherB\n", encoding="utf-8")
    with pytest.raises(ContractError, match="hash"):
        read_cached_contract(tmp_path)


def test_build_cli_runs_end_to_end_from_local_backed_inputs(monkeypatch, tmp_path):
    panel_bytes = b"GeneA\nGeneB\n"
    normalized = sha256(b"GeneA\nGeneB").hexdigest()
    index = {"T1:val": {"key": "T1:val", "task": "T1", "n_genes": 2, "genes_file": "panel.txt", "genes_sha256": normalized[:16], "needs_coords": False, "min_cells": 1, "max_cells": 2}}

    class Response:
        def __init__(self, payload): self.payload = payload
        def read(self): return self.payload
        def __enter__(self): return self
        def __exit__(self, *_): return False

    payloads = [json.dumps(index).encode(), panel_bytes]
    monkeypatch.setattr(workflow, "urlopen", lambda *args, **kwargs: Response(payloads.pop(0)))
    cache = tmp_path / "cache"
    workflow.fetch_official_contract(cache)
    early, late, _, _ = make_synthetic_training_inputs(["GeneA", "GeneB"])
    early_path, late_path = tmp_path / "E8.5.h5ad", tmp_path / "E9.5.h5ad"
    early.write_h5ad(early_path)
    late.write_h5ad(late_path)
    output = tmp_path / "candidate.h5ad"
    assert cli_main([
        "build", "--early", str(early_path), "--late", str(late_path),
        "--cache", str(cache), "--n-cells", "2", "--output", str(output),
    ]) == 0
    built = ad.read_h5ad(output)
    assert built.shape == (2, 2)
    assert list(built.var_names) == ["GeneA", "GeneB"]
    assert "celltype" not in built.obs
    selection = json.loads(Path(str(output) + ".selection.json").read_text())
    assert selection == {
        "schema_version": 1,
        "selection_policy": "uniform_without_replacement",
        "sampling_seed": 0,
        "available_cells": 5,
        "selected_cells": 2,
        "selected_physical_indices_sha256": selection["selected_physical_indices_sha256"],
    }
    assert len(selection["selected_physical_indices_sha256"]) == 64
    assert "sampling" not in built.uns


def test_public_cli_accepts_an_arbitrary_small_fake_ordered_panel(tmp_path):
    panel_path = tmp_path / "fake_panel.txt"
    panel_path.write_text("zeta\nalpha\nmu\n")
    output = tmp_path / "synthetic_output.h5ad"
    assert cli_main(["demo", "--panel", str(panel_path), "--n-cells", "5", "--output", str(output)]) == 0
    result = ad.read_h5ad(output)
    assert list(result.var_names) == ["zeta", "alpha", "mu"]
    assert result.shape == (5, 3)
    sidecar = json.loads(Path(str(output) + ".selection.json").read_text())
    assert sidecar["schema_version"] == 1
    assert sidecar["selection_policy"] == "uniform_without_replacement"
    assert sidecar["sampling_seed"] == 0
    assert "selected_input_row_indices" not in sidecar
    assert len(sidecar["selected_physical_indices_sha256"]) == 64
    assert "sampling" not in result.uns


def test_demo_cli_refuses_to_overwrite_existing_output(tmp_path, capsys):
    output = tmp_path / "keep_me.h5ad"
    output.write_bytes(b"sentinel")
    assert cli_main(["demo", "--output", str(output)]) == 2
    assert output.read_bytes() == b"sentinel"
    assert "Refusing to overwrite" in capsys.readouterr().err


def test_demo_cli_refuses_to_overwrite_existing_selection_sidecar(tmp_path, capsys):
    output = tmp_path / "new_output.h5ad"
    sidecar = Path(str(output) + ".selection.json")
    sidecar.write_text("sentinel", encoding="utf-8")
    assert cli_main(["demo", "--output", str(output)]) == 2
    assert not output.exists()
    assert sidecar.read_text(encoding="utf-8") == "sentinel"
    assert "Refusing to overwrite" in capsys.readouterr().err


def test_generated_t1_like_thousands_gene_fixture_stays_sparse_and_ordered():
    panel = [f"synthetic_gene_{i:05d}" for i in range(2_048)]
    early, late, returned_panel, contract = make_synthetic_training_inputs(panel)
    assert sparse.issparse(early.X)
    assert sparse.issparse(late.X)
    assert "celltype" in early.obs and "celltype" in late.obs
    assert "spatial_3D" not in early.obsm and "spatial_3D" not in late.obsm
    prediction = build_pseudobulk_shift_prediction(early, late, returned_panel, 3, contract)
    report = validate_task1_prediction(prediction, returned_panel, contract)
    assert report.ok
    assert sparse.isspmatrix_csr(prediction.X)
    assert list(prediction.var_names) == panel
    assert prediction.shape == (3, 2_048)


def test_notebook_cells_have_nonempty_unique_nbformat_ids():
    notebook = json.loads((Path(__file__).parents[1] / "notebooks" / "task1_baseline_to_submission.ipynb").read_text())
    ids = [cell.get("id") for cell in notebook["cells"]]
    assert all(isinstance(cell_id, str) and cell_id for cell_id in ids)
    assert len(ids) == len(set(ids))


def _shuffled_inputs(matrix_format: str):
    """Small independent fixture: shuffled observations, variables, and panel."""
    panel = ["g2", "g0", "g3", "g1"]
    early_genes = ["g0", "g1", "g2", "g3"]
    late_genes = ["g3", "g1", "g0", "g2"]
    early_values = np.array([
        [1, 4, 3, 2], [3, 2, 1, 4],  # A
        [5, 4, 3, 1], [2, 3, 4, 5],  # B
    ], dtype=np.float32)
    late_values = np.array([
        [6, 3, 4, 2],  # B
        [3, 5, 2, 4],  # A
        [7, 1, 5, 3],  # B
        [2, 6, 1, 5],  # A
        [4, 2, 8, 1],  # C, late-only
        [5, 4, 6, 2],  # B
    ], dtype=np.float32)
    if matrix_format == "dense":
        early_x, late_x = early_values, late_values
    elif matrix_format == "csr":
        early_x, late_x = sparse.csr_matrix(early_values), sparse.csr_matrix(late_values)
    elif matrix_format == "csc":
        early_x, late_x = sparse.csc_matrix(early_values), sparse.csc_matrix(late_values)
    else:  # pragma: no cover - test-only guard
        raise AssertionError(matrix_format)
    early = ad.AnnData(
        early_x,
        obs=pd.DataFrame({"celltype": ["A", "A", "B", "B"]}, index=["e4", "e1", "e3", "e2"]),
        var=pd.DataFrame(index=early_genes),
    )
    late = ad.AnnData(
        late_x,
        obs=pd.DataFrame({"celltype": ["B", "A", "B", "A", "C", "B"]}, index=["m_cell", "a_cell", "z_cell", "b_cell", "c_cell", "y_cell"]),
        var=pd.DataFrame(index=late_genes),
    )
    contract = {"key": "SYNTHETIC:T1:val", "task": "T1", "n_genes": 4, "min_cells": 4, "max_cells": 4, "needs_coords": False}
    return early, late, panel, contract, early_values, late_values, early_genes, late_genes


def _independent_reference(early_values, late_values, early_genes, late_genes, panel, early_types, late_types, selected):
    """Hand reference kept separate from workflow internals and backed access code."""
    early_cols = [early_genes.index(gene) for gene in panel]
    late_cols = [late_genes.index(gene) for gene in panel]
    output = late_values[selected][:, late_cols].astype(np.float32, copy=True)
    for celltype in sorted(set(early_types) & set(late_types)):
        delta = (
            late_values[np.asarray(late_types) == celltype][:, late_cols].mean(axis=0)
            - early_values[np.asarray(early_types) == celltype][:, early_cols].mean(axis=0)
        )
        output[np.asarray(late_types)[selected] == celltype] += delta
    return np.maximum(output, 0)


@pytest.mark.parametrize("matrix_format", ["dense", "csr", "csc"])
def test_actual_anndata_in_memory_and_backed_shuffled_inputs_match_hand_reference(tmp_path, matrix_format):
    early, late, panel, contract, early_values, late_values, early_genes, late_genes = _shuffled_inputs(matrix_format)
    seed = 23
    selected = np.random.default_rng(seed).choice(late.n_obs, size=4, replace=False)
    expected = _independent_reference(
        early_values, late_values, early_genes, late_genes, panel,
        early.obs["celltype"].tolist(), late.obs["celltype"].tolist(), selected,
    )

    in_memory_record: dict[str, object] = {}
    in_memory = build_pseudobulk_shift_prediction(
        early, late, panel, 4, contract, sampling_seed=seed,
        sampling_report=in_memory_record, chunk_rows=2,
    )
    np.testing.assert_allclose(in_memory.X.toarray(), expected, rtol=1e-6, atol=2e-7)
    assert in_memory_record["selected_input_row_indices"] == selected.tolist()
    assert "sampling" not in in_memory.uns

    early_path, late_path = tmp_path / f"early_{matrix_format}.h5ad", tmp_path / f"late_{matrix_format}.h5ad"
    early.write_h5ad(early_path)
    late.write_h5ad(late_path)
    backed_early = ad.read_h5ad(early_path, backed="r")
    backed_late = ad.read_h5ad(late_path, backed="r")
    try:
        source_contract = {**contract, "min_cells": late.n_obs, "max_cells": late.n_obs}
        input_report = validate_task1_prediction(backed_late, late_genes, source_contract)
        assert input_report.ok
        backed_record: dict[str, object] = {}
        backed = build_pseudobulk_shift_prediction(
            backed_early, backed_late, panel, 4, contract, sampling_seed=seed,
            sampling_report=backed_record, chunk_rows=2,
        )
    finally:
        backed_early.file.close()
        backed_late.file.close()
    np.testing.assert_allclose(backed.X.toarray(), expected, rtol=1e-6, atol=2e-7)
    assert backed_record == in_memory_record


def test_backed_dense_unsorted_rows_restore_requested_order_across_multiple_chunks(tmp_path):
    _, late, panel, _, _, late_values, _, late_genes = _shuffled_inputs("dense")
    late_path = tmp_path / "late_dense.h5ad"
    late.write_h5ad(late_path)
    backed_late = ad.read_h5ad(late_path, backed="r")
    requested_rows = np.array([4, 1, 5, 0], dtype=np.int64)
    positions = pd.Index(late_genes).get_indexer(panel)
    try:
        actual = workflow._sparse_safe_shifted_rows(
            backed_late, requested_rows, positions, late.obs["celltype"].to_numpy(), {}, chunk_rows=2,
        )
    finally:
        backed_late.file.close()
    np.testing.assert_allclose(actual.toarray(), late_values[requested_rows][:, positions])


def test_fixed_seed_sampling_is_reproducible_and_identifier_prefixes_do_not_select_a_population():
    n_per_type = 1_000
    early = ad.AnnData(
        np.array([[1.0], [2.0]], dtype=np.float32),
        obs=pd.DataFrame({"celltype": ["A", "B"]}, index=["early_a", "early_b"]),
        var=pd.DataFrame(index=["g"]),
    )
    late_values = np.concatenate([
        np.ones((n_per_type, 1), dtype=np.float32),
        np.full((n_per_type, 1), 2.0, dtype=np.float32),
    ])
    types = np.array(["A"] * n_per_type + ["B"] * n_per_type)
    late = ad.AnnData(
        late_values,
        obs=pd.DataFrame({"celltype": types}, index=[f"A_{i:04d}" for i in range(n_per_type)] + [f"B_{i:04d}" for i in range(n_per_type)]),
        var=pd.DataFrame(index=["g"]),
    )
    contract = {"key": "SYNTHETIC:T1:val", "task": "T1", "n_genes": 1, "min_cells": 1_000, "max_cells": 1_000, "needs_coords": False}
    first_record: dict[str, object] = {}
    first = build_pseudobulk_shift_prediction(early, late, ["g"], 1_000, contract, sampling_report=first_record)
    second_record: dict[str, object] = {}
    second = build_pseudobulk_shift_prediction(early, late, ["g"], 1_000, contract, sampling_report=second_record)
    np.testing.assert_array_equal(first.X.toarray(), second.X.toarray())
    assert first_record == second_record
    selected = np.asarray(first_record["selected_input_row_indices"])
    assert not np.array_equal(selected, np.arange(1_000))  # guards against lexical-prefix selection
    observed_types = types[selected]
    assert set(observed_types) == {"A", "B"}

    renamed = late.copy()
    renamed.obs_names = [f"renamed_{i:04d}" for i in range(renamed.n_obs)]
    renamed_record: dict[str, object] = {}
    renamed_prediction = build_pseudobulk_shift_prediction(
        early, renamed, ["g"], 1_000, contract, sampling_report=renamed_record,
    )
    np.testing.assert_array_equal(first.X.toarray(), renamed_prediction.X.toarray())
    assert renamed_record == first_record
