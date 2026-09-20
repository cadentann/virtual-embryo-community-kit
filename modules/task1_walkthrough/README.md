# Task 1: baseline → current-contract submission walkthrough

This is a small, synthetic-only bridge between a pseudobulk-shift baseline demonstration and the **current Task 1 submission contract** of the [Virtual Embryo Challenge](https://virtualembryo.ai/challenge). It is designed for a participant who can run a notebook but does not yet know which parts of a demo are illustrative and which parts are required for a real submission.

It creates and validates a *submission-shaped* `.h5ad` file. It does **not** produce a competitive model, score a prediction, obtain withheld truth, replace the published standalone local scorer (`veckit`) or official challenge server, or decide data eligibility.

## The important distinction

| Step | What this resource does | What it cannot do |
| --- | --- | --- |
| Offline demo | Builds a deterministic synthetic `.h5ad` with Task-1-like structure | Uses no challenge data; it intentionally fails the live board contract and must not be uploaded |
| Real prediction | Shows how to fetch the current board metadata/panel, align data, and fail safely | Does not infer E10.5 truth |
| Local score | Explains the published standalone local scorer (`veckit`) when a permissible target is available | Cannot locally score currently withheld validation or never-distributed test truth |
| Leaderboard | Leaves score evaluation to the official service | Cannot reproduce server validation without released target data |

As checked on 2026-09-20, Task 1 uses E8.5 and E9.5 released training stages, E10.5 validation, and E12.5 test; the Task-1 release has no E9.25. E10.5 validation truth is currently withheld and evaluated through the leaderboard; the official schedule says it will be released at the final phase on 2026-10-20. E12.5 test truth is not distributed. The current Task-1 contract has 32,285 ordered genes, no spatial coordinates, and a `T1:val` board. The machine-readable `T1:val` index and data-page prose specify 1,000–5,118 cells, while evaluation prose says at least 1,000 with no cap. This companion enforces the dated index as a conservative local preflight policy; it does not resolve the official disagreement, and the upload service/organizers remain authoritative. Check the [Task 1 page](https://virtualembryo.ai/challenge/tasks/temporal), [timeline](https://virtualembryo.ai/challenge/timeline), [data contract](https://virtualembryo.ai/challenge/data), [panel index](https://virtualembryo.ai/challenge/panels/index.json), and [submission page](https://virtualembryo.ai/challenge/account/submissions) again before upload: requirements can change.

## What this corrects without republishing the baseline notebook

The pedagogical method is credited to the [published `pseudobulk_shift` baseline and linked organizer notebook](https://virtualembryo.ai/challenge/baselines). This is original companion code/prose, not a redistributed organizer notebook or a new predictive method.

| Published demo assumption captured 2026-09-10 | Current-contract bridge here |
| --- | --- |
| T1 real-use recipe names E9.25 | Uses only released E8.5 and E9.5 inputs |
| Example local command names an E10.5 target | Requires no withheld target; server evaluation stays separate |
| Gene intersection is silent | Fetches the current board panel and fails on missing, extra, duplicate, or reordered output genes |
| All E9.5 cells are copied | Uses current board bounds and a fixed-seed uniform sample without replacement |
| Whole sparse matrices may be densified | Aggregates sparsely and reads/applies shifts in bounded physical row blocks, restoring the sampled order |

## Five-minute offline quickstart

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e '.[test]'
vec-t1-walkthrough demo --output artifacts/synthetic_t1_prediction.h5ad
pytest -q
```

Expected terminal output includes:

```text
Wrote synthetic, submission-shaped demonstration: artifacts/synthetic_t1_prediction.h5ad
It is not challenge data and must not be uploaded: it uses SYNTHETIC:T1:val.
```

The demo refuses to overwrite an existing output. Choose a new path or move the old artifact first.

The executable notebook [`notebooks/task1_baseline_to_submission.ipynb`](notebooks/task1_baseline_to_submission.ipynb) performs the same offline demonstration and is intentionally safe to run without network access. To run it in Jupyter, install the optional notebook runtime and start Jupyter from this project directory:

```bash
pip install -e '.[notebook]'
jupyter lab
```

Local Jupyter is the only supported notebook path. Google Colab is not assessed or supported.

For a non-interactive reproducibility run from the project directory:

```bash
jupyter nbconvert --to notebook --execute notebooks/task1_baseline_to_submission.ipynb \
  --output /tmp/task1_baseline_to_submission.executed.ipynb --ExecutePreprocessor.timeout=180
```

For a black-box format test with an arbitrary **small fake** ordered panel (never a challenge panel), use:

```bash
printf 'zeta\nalpha\nmu\n' > fake_panel.txt
vec-t1-walkthrough demo --panel fake_panel.txt --n-cells 5 --output artifacts/fake_panel_demo.h5ad
```

The output uses precisely that order and a `SYNTHETIC:T1:val` contract. It is useful for checking integration plumbing, not for challenge submission.

## Moving from the demo to a real Task-1 output

1. Train or otherwise build a prediction using permitted inputs. The walkthrough's reference pseudobulk shift calculates a separate E8.5→E9.5 mean-expression delta for each `obs["celltype"]` shared by both stages, adds that delta to each E9.5 cell of its type, and clips at zero. It therefore retains within-type E9.5 variation rather than replacing every cell with one global profile. For an E9.5-only type, it conservatively keeps the E9.5 cell unchanged; this is documented behavior, not a biological inference. To meet the board's output shape, the example selects `n_cells` uniformly without replacement with fixed seed `0`; it does not privilege lexical cell identifiers. This output-preparation choice is separate from the published pseudobulk-shift method. For a submission, `n_obs` is a sample-size choice subject to the current board bounds, not a measurement or prediction of embryo size. This is not presented as a strong model.
2. Fetch **current** official metadata and the ordered board panel. This is an explicit network step and saves the URLs, timestamp, and SHA-256 provenance locally.

   ```bash
   vec-t1-walkthrough fetch-contract --cache official_contract_cache
   ```

3. Build a candidate from your **local, permitted copies** of the released Task-1 E8.5 and E9.5 files. The command reads them in backed mode, uses exactly the cached panel order and current cell bounds, writes a new path, and refuses to overwrite an existing file:

   ```bash
   vec-t1-walkthrough build \
     --early /path/to/released/E8.5_RNA.h5ad \
     --late /path/to/released/E9.5_RNA.h5ad \
     --cache official_contract_cache \
     --n-cells 1000 \
     --seed 0 \
     --output artifacts/task1_candidate.h5ad
   ```

   Do not assume a released training `.h5ad` has the upload order. The helper rejects duplicate input/panel genes and missing required input genes. It uses `min_cells`/`max_cells` from the fetched metadata rather than a hardcoded limit. The command writes `artifacts/task1_candidate.h5ad.selection.json`, recording the uniform-without-replacement policy, seed, available/selected counts, and `selected_physical_indices_sha256`, the SHA-256 digest of selected integer input-row indices. It deliberately omits cell identifiers and keeps the record out of the prediction's AnnData metadata. Choose and record a different integer seed only when you intend a different sample. The 1,000–5,118 bounds are the dated index-derived companion policy described above, not a resolution of the evaluation prose's “at least 1,000”/“no cap” wording.
4. The command writes a minimal `AnnData`: a two-dimensional, finite, non-negative `float32` `.X`, ordered `var_names`, and only required metadata. The local checker accepts only real numeric expression dtypes; it rejects boolean presence masks, complex values, and string/object matrices rather than coercing them. Its stricter minimal-output policy omits and rejects `spatial_3D`, while current server behavior remains authoritative; it also omits training-only `celltype` metadata.
5. Reopen the written file through the structural preflight:

   ```bash
   vec-t1-walkthrough validate path/to/prediction.h5ad --cache official_contract_cache
   ```

   `PASS` means only that this companion tool's local structural checks passed. It does **not** establish that `.X` uses the required challenge normalization: raw counts or another finite, non-negative scale can pass locally and still be unsuitable for scoring. Confirm preprocessing against the current [data contract](https://virtualembryo.ai/challenge/data). Task 1 does not require or use spatial coordinates; this companion deliberately rejects `spatial_3D` to keep its minimal output, while current upload prose says `obsm` is ignored. That stricter companion policy is not asserted to mirror server acceptance. The official upload checks and challenge server evaluation, including the current contract, are authoritative.
6. Keep three operations distinct. This package provides a **local structural preflight**. The published standalone local scorer (`veckit`) requires participant-supplied, permissible target/reference files; for a training-stage smoke test, its documented pattern is:

   ```bash
   veckit --task T1 \
     --input prediction_for_E9.5.h5ad \
     --target /path/to/released/E9.5_RNA.h5ad \
     --reference /path/to/released/E8.5_RNA.h5ad
   ```

   Record the installed `veckit` version/commit because its public implementation is evolving. This does not preview a leaderboard score. As of 2026-09-20, withheld E10.5 validation truth is evaluated by the **official challenge server**; E12.5 test truth is never distributed. See the official [evaluation](https://virtualembryo.ai/challenge/evaluation), [baseline/scorer references](https://virtualembryo.ai/challenge/baselines), and [rules](https://virtualembryo.ai/challenge/rules).

## Why gene order is a release gate

`var_names` is an ordered contract, not just a set. The validator distinguishes the following failures:

- same genes but wrong order;
- missing required genes;
- unexpected genes;
- duplicate genes;
- variable count mismatch.

It also rejects out-of-range cell counts, non-finite values, negative values, missing/non-2D/non-real-numeric `.X`, and `spatial_3D` under this companion's stricter minimal-output policy. Boolean `.X` is deliberately rejected because a presence mask is not a real-valued expression matrix. Current official evaluation prose says `obsm` is ignored, so the `spatial_3D` rejection is not presented as a server rejection. These local choices are not a claim that the server rejects every same condition. It does not silently reorder an output. Any optional reorder facility in an official utility is not a substitute for generating the correct ordered output yourself.

## Real data and permissions

This release contains only deterministic synthetic matrices and no challenge `.h5ad`, target, or extracted biological content. The dataset is unpublished. Its use in a publication or public presentation requires prior approval from Dr. Neil Chi and the data-generating team. Outside-Challenge use, reuse, publication, presentation, or redistribution also carries the acknowledgment and credit requirements in the current [Terms](https://virtualembryo.ai/challenge/terms). Ask the organizers about uses not clearly covered, including redistribution of derived artifacts. This is a factual summary, not legal advice. External public data may be permitted subject to the current rules and disclosure requirements, but this tool is not an eligibility oracle.

## Layout

```text
src/vec_t1_walkthrough/  utility and CLI
notebooks/               runnable synthetic tutorial
tests/                   structural/adversarial tests
```

Read [TESTING.md](TESTING.md), [KNOWN_LIMITATIONS.md](KNOWN_LIMITATIONS.md), [CHANGELOG.md](CHANGELOG.md), and the repository-level [license](../../LICENSE) and [attribution](../../ATTRIBUTION.md) before adapting this work.
