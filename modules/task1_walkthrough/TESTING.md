# Testing

Run in a fresh environment:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e '.[test]'
vec-t1-walkthrough demo --output artifacts/synthetic_t1_prediction.h5ad
pytest -q
```

The preferred test installation is `pip install -e '.[test]'`, which keeps test-only `pytest` separate from runtime dependencies. `requirements.txt` is a convenient runtime-plus-pytest alternative; the Jupyter UI/executor remains the explicit `.[notebook]` extra.

Local Jupyter is the only supported notebook environment. Colab is not assessed or supported in this candidate.

The closeout G1 run on 2026-09-11 used Python 3.11.0, `anndata 0.12.19`, `certifi 2026.7.22`, `numpy 2.4.4`, `pandas 2.3.3`, `scipy 1.17.1`, and `pytest 9.1.1`; all 42 unit tests passed. Final frozen-bundle acceptance and its exact dependency versions are recorded separately in the private evidence bundle. The declared upper bounds are intentional: `anndata>=0.10,<0.13`, `certifi>=2024,<2028`, `numpy>=1.26,<3`, `pandas>=2,<3`, `scipy>=1.11,<2`, and test `pytest>=8,<10`.

The tests use only synthetic data. They cover hand-calculated per-cell-type shifts, retained within-type variation, E9.5-only-type unchanged behavior, a generated 2,048-gene Task-1-like sparse fixture, a minimal sparse float32 output, exact panel order, missing/extra/duplicate genes, NaN and negative values, dense and sparse `.X`, metadata-derived cell bounds, requests larger than the available late-stage cell pool, output-overwrite refusal, missing `celltype`, missing input genes, and non-mutation by validation. They also write and reopen actual AnnData dense, CSR, and CSC inputs in backed mode; compare shuffled observations/panel order and multi-chunk reads to an independent hand reference; exercise an unsorted dense-backed row request with order restoration; and show fixed-seed uniform sampling is reproducible and not controlled by identifier prefixes. G1 regressions reject complex, string, object, boolean, missing, and non-2D expression inputs without throwing; reject complex/boolean builder inputs before real casts; and exercise clean CLI build/validate failures with real `.h5ad` files and closed handles. Boolean matrices are deliberately invalid: a presence mask is not a real-valued expression matrix. Per-cell-type aggregation and shift application read increasing physical rows in bounded blocks; only a bounded output block is densified before returning CSR. The builder writes sampling policy/seed only to a JSON sidecar, never prediction metadata. The validator examines sparse `.data` without densifying it; dense inputs are scanned in row chunks.

Not covered: organizer server behavior, a production challenge-data run, current network availability in every environment, the standalone scorer, or any hidden target. The notebook's code paths are equivalent to the exercised library calls. Notebook execution uses the optional `pip install -e '.[notebook]'` runtime rather than burdening CLI-only users with Jupyter. Re-fetch and inspect the official contract immediately before an actual submission.
