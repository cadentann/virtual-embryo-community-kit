"""Command line entry points for the walkthrough."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from urllib.error import URLError

import anndata as ad

from .workflow import (
    ContractError, build_pseudobulk_shift_prediction, fetch_official_contract,
    make_synthetic_training_inputs, read_cached_contract, validate_task1_prediction,
)


def _sampling_sidecar_path(output: Path) -> Path:
    return Path(str(output) + ".selection.json")


def _require_new_output_paths(output: Path) -> Path:
    sidecar = _sampling_sidecar_path(output)
    existing = [path for path in (output, sidecar) if path.exists()]
    if existing:
        raise ContractError("Refusing to overwrite existing output or sampling sidecar: " + ", ".join(map(str, existing)))
    return sidecar


def _write_sampling_sidecar(path: Path, sampling: dict[str, object]) -> None:
    """Write the public selection record without source identifiers or AnnData metadata."""
    payload = {
        "schema_version": 1,
        "selection_policy": sampling["method"],
        "sampling_seed": sampling["seed"],
        "available_cells": sampling["n_available_late_cells"],
        "selected_cells": sampling["n_selected_cells"],
        "selected_physical_indices_sha256": sampling["selected_input_row_indices_sha256"],
    }
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def _demo(output: Path, panel_path: Path | None, n_cells: int, sampling_seed: int) -> int:
    sidecar = _require_new_output_paths(output)
    panel = None if panel_path is None else [line.strip() for line in panel_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    early, late, panel, contract = make_synthetic_training_inputs(panel)
    sampling: dict[str, object] = {}
    pred = build_pseudobulk_shift_prediction(
        early, late, panel, n_cells=n_cells, contract=contract,
        sampling_seed=sampling_seed, sampling_report=sampling,
    )
    report = validate_task1_prediction(pred, panel, contract)
    report.require_ok()
    output.parent.mkdir(parents=True, exist_ok=True)
    pred.write_h5ad(output)
    _write_sampling_sidecar(sidecar, sampling)
    print(f"Wrote synthetic, submission-shaped demonstration: {output}")
    print(f"Wrote sampling sidecar: {sidecar}")
    print("It is not challenge data and must not be uploaded: it intentionally fails the live board contract.")
    return 0


def _validate(path: Path, cache: Path) -> int:
    contract, panel = read_cached_contract(cache)
    pred = ad.read_h5ad(path, backed="r")
    try:
        report = validate_task1_prediction(pred, panel, contract)
    finally:
        pred.file.close()
    for issue in report.errors:
        print("ERROR:", issue)
    for issue in report.warnings:
        print("WARNING:", issue)
    print("PASS" if report.ok else "FAIL", "— structural checks only; official server remains authoritative.")
    return 0 if report.ok else 2


def _build(early_path: Path, late_path: Path, cache: Path, output: Path, n_cells: int, sampling_seed: int) -> int:
    sidecar = _require_new_output_paths(output)
    contract, panel = read_cached_contract(cache)
    early = late = None
    sampling: dict[str, object] = {}
    try:
        early = ad.read_h5ad(early_path, backed="r")
        late = ad.read_h5ad(late_path, backed="r")
        prediction = build_pseudobulk_shift_prediction(
            early, late, panel, n_cells=n_cells, contract=contract,
            sampling_seed=sampling_seed, sampling_report=sampling,
        )
    finally:
        for data in (early, late):
            if data is not None:
                data.file.close()
    report = validate_task1_prediction(prediction, panel, contract)
    report.require_ok()
    output.parent.mkdir(parents=True, exist_ok=True)
    prediction.write_h5ad(output)
    _write_sampling_sidecar(sidecar, sampling)
    print(f"Wrote Task-1-shaped candidate from local E8.5/E9.5 inputs: {output}")
    print(f"Wrote sampling sidecar: {sidecar}")
    print("Reopen and run `validate` before upload; structural checks do not verify normalization or model quality.")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Synthetic-only Task 1 submission walkthrough")
    sub = parser.add_subparsers(dest="command", required=True)
    p_demo = sub.add_parser("demo", help="write an offline synthetic output")
    p_demo.add_argument("--output", type=Path, default=Path("artifacts/synthetic_t1_prediction.h5ad"))
    p_demo.add_argument("--panel", type=Path, help="optional small fake ordered panel; never pass a challenge panel for this synthetic demo")
    p_demo.add_argument("--n-cells", type=int, default=3, help="synthetic output cells (contract range: 3–5)")
    p_demo.add_argument("--seed", type=int, default=0, dest="sampling_seed", help="fixed seed for uniform sampling without replacement")
    p_fetch = sub.add_parser("fetch-contract", help="opt in to cache live official T1:val metadata/panel")
    p_fetch.add_argument("--cache", type=Path, default=Path("official_contract_cache"))
    p_validate = sub.add_parser("validate", help="validate a candidate against a previously fetched contract")
    p_validate.add_argument("prediction", type=Path)
    p_validate.add_argument("--cache", type=Path, default=Path("official_contract_cache"))
    p_build = sub.add_parser("build", help="build from local released E8.5/E9.5 files using a verified cached contract")
    p_build.add_argument("--early", type=Path, required=True, help="local released Task-1 E8.5 .h5ad")
    p_build.add_argument("--late", type=Path, required=True, help="local released Task-1 E9.5 .h5ad")
    p_build.add_argument("--n-cells", type=int, required=True, help="output cells; must satisfy the cached board bounds")
    p_build.add_argument("--output", type=Path, required=True, help="new output path; existing files are never overwritten")
    p_build.add_argument("--cache", type=Path, default=Path("official_contract_cache"))
    p_build.add_argument("--seed", type=int, default=0, dest="sampling_seed", help="fixed seed for uniform sampling without replacement")
    args = parser.parse_args(argv)
    try:
        if args.command == "demo":
            return _demo(args.output, args.panel, args.n_cells, args.sampling_seed)
        if args.command == "fetch-contract":
            contract, panel = fetch_official_contract(args.cache)
            print(f"Cached {contract['key']} with {len(panel)} ordered genes in {args.cache}.")
            return 0
        if args.command == "build":
            return _build(args.early, args.late, args.cache, args.output, args.n_cells, args.sampling_seed)
        return _validate(args.prediction, args.cache)
    except (ContractError, OSError, URLError, ValueError, TypeError, AttributeError, KeyError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        if args.command in {"fetch-contract", "validate", "build"}:
            print("No contract-backed result was accepted. Check the error and re-fetch the official cache if needed.", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
