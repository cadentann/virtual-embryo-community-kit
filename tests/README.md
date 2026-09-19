# Public integration tests

From the repository root, create and activate a fresh virtual environment, install both modules with their test and notebook extras, then run:

```bash
python -m pip install -e './modules/task1_walkthrough[test,notebook]' -e './modules/data_inspector[test,notebook]'
PYTHONDONTWRITEBYTECODE=1 python -m pytest -p no:cacheprovider -q tests
```

The tests invoke the documented public CLIs with the active interpreter, run the synthetic A → B path, validate notebook hygiene, and scan distributable source content for private/data residue. Local packaging caches and documented synthetic-output directories created after installation/quickstart are ignored by this runtime check; the separate pre-ZIP release gate rejects those directories before shipment. The tests do not contact challenge services or use challenge data. The bounded large sparse sentinel belongs to the private verification run and is generated at test time there; no large artifact is committed.

The current public tests do not require `tests/blackbox_config.json`. To extend the portable private-style harness, copy `tests/blackbox_config.example.json` to `tests/blackbox_config.json`; that local file is ignored by Git.
