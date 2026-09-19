# Known limitations

- This is an onboarding and format-preflight resource, not a competitive method.
- The reference shift needs `obs["celltype"]` in both released training-stage inputs. It only computes shifts for types present in both; E9.5-only types remain unchanged rather than receiving an invented historical shift. This is a conservative implementation detail, not a biological interpretation.
- The bundled notebook uses a fake `SYNTHETIC:T1:val` contract with four genes. It intentionally fails the live board contract and must not be uploaded.
- Selecting a board-bounded number of late-stage rows is output preparation, not part of the pseudobulk-shift method. This package defaults to a reproducible uniform sample without replacement (`--seed 0`) and writes its policy, seed, counts, and a digest of selected integer row indices to `<output>.selection.json`, not the prediction. A different seed deliberately produces a different sample.
- `fetch-contract` requires live access to the official panel service. Its cache can become stale; inspect `provenance.json` and re-fetch before submission.
- A passing check is not organizer approval, data-eligibility approval, an official format acceptance, or a score.
- The local checker deliberately accepts only two-dimensional, real numeric expression matrices. Boolean presence masks, complex matrices, and string/object encodings are rejected rather than coerced; this type policy is a conservative local preflight boundary, not an assertion that it exhaustively mirrors server behavior.
- The structural checker does not reimplement the published standalone local scorer (`veckit`) or official server evaluation and does not assess biological plausibility.
- A full-size output can still require substantial memory, especially when additive shifts turn many implicit sparse zeros into explicit values. This package keeps each conversion to bounded row blocks and never densifies the full input matrix, but it accumulates the CSR output in memory and cannot remove the intrinsic cost. A deliberately near-dense 1,000 × 32,285 synthetic stress case used about 767 MiB peak RAM and wrote a 248 MiB file on the tested macOS host; budget at least 1 GiB RAM and 300 MiB free disk, with more headroom for real inputs.
- This resource makes no claim about cell-type biology or E10.5/E12.5 expression. It never accesses held-out truth.
