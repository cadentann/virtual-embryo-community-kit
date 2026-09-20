# VEC First-Submission and Read-Only Preflight Kit

A companion community resource for the Virtual Embryo Challenge, with a synthetic Task-1 baseline-to-submission walkthrough and a read-only `.h5ad` preflight inspector. It is for participants preparing their first prediction or checking the structure of a local file.

This is an independent resource, not an organizer-endorsed tool.

Contribution A preserves and credits the published per-cell-type `pseudobulk_shift` teaching method; it is not a new model. Contribution B performs structural and descriptive QC; it is not the official scorer, a biological interpretation engine, an eligibility oracle, or an upload-service replacement. The modules are independently installable and do not import one another.

## Quickstart

Follow [START_HERE.md](START_HERE.md) for the synthetic A → B workflow. It uses only repository-relative paths and makes no challenge-data request.

## Scope and safety

A structural pass does not establish normalization, biological quality, official acceptance, eligibility, or score. No challenge data, hidden truth, cached panel snapshot, real-data report, or organizer notebook is bundled. Keep reports made from real permitted files private unless current terms and permissions allow sharing.

Official requirements can change. Recheck the [Task 1 page](https://virtualembryo.ai/challenge/tasks/temporal), [data/panel guidance](https://virtualembryo.ai/challenge/data), [panel index](https://virtualembryo.ai/challenge/panels/index.json), [submission guidance](https://virtualembryo.ai/challenge/account/submissions), and [rules/terms](https://virtualembryo.ai/challenge/rules) before use. The official pages, machine-readable metadata, and ordered T1 panel were checked on 2026-09-20. No downloaded official source snapshot is bundled here.

The supported notebook route is local Jupyter from each module directory after its editable install. Colab is not assessed or supported. The current T1 machine index and data-page prose say 1,000–5,118 cells; the evaluation prose says “at least 1,000” and “no cap.” The walkthrough treats the machine index as a dated local preflight input, not as a resolution of that official-page disagreement.

Original material that Caden Tan controls is licensed under [Apache License 2.0](LICENSE); see [NOTICE](NOTICE) and [ATTRIBUTION.md](ATTRIBUTION.md) for scope and provenance. Also see [CONTRIBUTING.md](CONTRIBUTING.md), module limitations, and [tests/README.md](tests/README.md).
