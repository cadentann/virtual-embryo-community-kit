# Start here: synthetic A → B quickstart

From a fresh clone:

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e './modules/task1_walkthrough[test,notebook]'
python -m pip install -e './modules/data_inspector[test,notebook]'

mkdir -p local_outputs
vec-t1-walkthrough demo --output local_outputs/task1_synthetic_prediction.h5ad --seed 20260911
vec-inspect local_outputs/task1_synthetic_prediction.h5ad --report-dir local_outputs/inspection
```

The A command also writes `local_outputs/task1_synthetic_prediction.h5ad.selection.json`, recording the fixed seed, uniform-without-replacement policy, counts, and a source-local selection fingerprint. The prediction and reports are synthetic and are ignored by Git. The B reports are Markdown and JSON; inspect the terminal `attention` section as well as the files.

For local notebooks, install the relevant `.[notebook]` extra, then start Jupyter from the module directory. This release supports local Jupyter only; Google Colab is **not assessed or supported**. No public URL or Colab bootstrap is included.

In the activated shell, starting at the repository root:

```bash
cd modules/task1_walkthrough
python -m jupyter lab notebooks/task1_baseline_to_submission.ipynb
```

Jupyter stays running in that shell. In a **new terminal**, navigate to the
repository root (the folder containing this file), activate the environment
created above, and run the public checks:

```bash
. .venv/bin/activate
python -m pytest -q modules/task1_walkthrough/tests modules/data_inspector/tests tests
```

For the inspector notebook, start from the repository root in an activated
shell, then run:

```bash
cd modules/data_inspector
python -m jupyter lab notebooks/vec_data_inspector_demo.ipynb
```

The quickstart installation above already includes the notebook extras. The
inspector notebook generates synthetic fixtures and reports under that
module's `reports/` directory.

The optional live-contract command is separate and networked:

```bash
vec-t1-walkthrough fetch-contract --cache local_outputs/official_contract_cache
```

It fetches current public metadata/panel files and writes provenance. Do not use it for the offline synthetic quickstart, and do not treat local structural validation as a server score or eligibility decision.
