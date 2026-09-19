# Synthetic example report

**Synthetic-only example.** This readable example was generated from the
deterministic `synthetic_dense_t1_like.h5ad` and
`synthetic_panel_reordered.txt` fixtures made by
`examples/make_synthetic_fixtures.py`. Paths below are deliberately made
package-relative; this file contains no challenge data or private local paths.

## File identity

- Path: `examples/generated/synthetic_dense_t1_like.h5ad`
- Size: 12,800 bytes
- SHA-256 local fingerprint: `a6b15ee96e14bd0847919a1935415de31c68cb31800b9686bf5743dba69cf73b`
- Note: This is a local fingerprint, not an authenticated checksum.

## AnnData structure

- n_obs: 3
- n_vars: 4
- X storage: `dense_dataset`
- X dtype: `float32`
- obs columns: `celltype`
- var columns: none
- obsm, layers, and uns keys: none
- Axis lengths match X: `True`

## Expression QC

- Finite: `True`
- Negative values detected: `False`
- Finite min / max: `0.0 / 6.0`
- Sparsity estimate: `0.5`

## Cell-type and spatial summary

| Synthetic label | Count | Fraction |
| --- | ---: | ---: |
| type_a | 2 | 0.6667 |
| type_b | 1 | 0.3333 |

No `obsm/spatial_3D` field is present in this dense synthetic fixture.

## Gene-panel comparison

- Panel path: `examples/generated/synthetic_panel_reordered.txt`
- Exact ordered match: `False`
- Same unique set, different order: `True`
- Missing from file / extra in file: `0 / 0`
- Duplicate file names / duplicate panel names: `0 / 0`
- Result: `reorder_plan.csv` can be written as a source-index plan only; the
  source fixture remains unchanged.

## Scope

This is a read-only structural/QC example, not a score, authorization, or
normalization guarantee. Finite non-negative values can still be raw counts or
the wrong scale. The synthetic names, values, labels, and dimensions above are
fabricated solely for this example.
