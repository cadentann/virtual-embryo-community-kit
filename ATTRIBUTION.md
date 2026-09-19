# Attribution

This project is a companion resource for the [Virtual Embryo Challenge](https://virtualembryo.ai/challenge). The organizers and official sources remain authoritative for challenge rules, data, panels, submission format, and terms: [Task 1](https://virtualembryo.ai/challenge/tasks/temporal), [data](https://virtualembryo.ai/challenge/data), [panels](https://virtualembryo.ai/challenge/panels/index.json), [submissions](https://virtualembryo.ai/challenge/account/submissions), [rules](https://virtualembryo.ai/challenge/rules), [terms](https://virtualembryo.ai/challenge/terms), and [community](https://virtualembryo.ai/challenge/community). Official pages were checked on 2026-09-11; requirements may change.

Contribution A credits the published organizer `pseudobulk_shift` baseline as its pedagogical method. The organizer notebook itself is not redistributed, and A does not claim a new predictive method.

The public [`aristoteleo/veckit`](https://github.com/aristoteleo/veckit) project is credited as the existing standalone local scoring implementation. No `veckit` code, sample data, scorer, or metrics are bundled here. [`cellgeni/h5ad-cli`](https://github.com/cellgeni/h5ad-cli) and [`vals/anndata-design-inspector`](https://github.com/vals/anndata-design-inspector) are related public tools; B has a narrower VEC-oriented, non-mutating QC/report scope.

Subject to final similarity/provenance review, the code, prose, notebooks, synthetic fixture generator, and committed synthetic example reports in this candidate were created for this project. Challenge data, official panels/snapshots, mirrored source pages, organizer material, and reports derived from real challenge files are excluded and are not covered by the proposed project license.

Runtime dependencies include Python, NumPy, SciPy, pandas, h5py, AnnData, certifi, pytest, and optional Jupyter tooling. Their own project metadata and licenses govern those dependencies; this repository does not relicense them.
