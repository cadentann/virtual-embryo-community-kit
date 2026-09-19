# Changelog

## G1 input-type hardening — 2026-09-11

- Reject missing, non-2D, boolean, complex, string, and object expression matrices through validation reports before numeric operations; builder inputs receive the same real-numeric check before float casts.
- Close backed validation/build handles on both success and expected failure, and add real AnnData/HDF5 API and CLI regressions for these boundaries.

## Targeted repair — 2026-09-11

- Replaced lexical-prefix cell selection with fixed-seed uniform sampling without replacement. A collision-protected `<output>.selection.json` records policy, seed, counts, and a digest of selected integer row indices, never prediction metadata or cell identifiers.
- Made dense and sparse backed input reads HDF5-safe by reading increasing physical row blocks and restoring requested output order; expanded real AnnData dense/CSR/CSC, shuffled-order, multi-chunk, hand-reference, and sampling-bias regressions.

## Post-review hardening — 2026-09-10

- Added the local E8.5/E9.5 `build` path, cache-integrity verification, normalization/phase/scorer distinctions, output collision protection, notebook runtime instructions, and a full-dimension synthetic stress benchmark.
- Expanded the unit suite with adversarial regression coverage.

## 0.1.0 — 2026-09-10

- Initial private release candidate.
- Added an offline synthetic Task-1 walkthrough, opt-in official contract fetch/cache, strict structural preflight, sparse-safe pseudobulk aggregation, CLI, notebook, and adversarial tests.
- No organizer notebook, challenge matrix, score, or held-out target is redistributed.
