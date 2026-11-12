# API reference

TODO: generate from the docstrings with mkdocstrings. Outline:

## Loading

- `vf.load(path, locus, features, *, seed=0, run_operations=True)` → `Repertoires`
- `vf.feature_blocks()`, `vf.FEATURE_BLOCKS`, `vf.FeatureBlock`

## `Repertoires` (steps 1–4)

- `overview()`, `with_metadata(*exprs)`, `qc(...)` → `QCResult`, `design(...)` → `StudyDesign`
- Attributes: `metadata`, `X`, `feature_names`, `features`, `clonotypes`

## `StudyDesign` (steps 5–11)

- `summary()`, `samples()`
- `balance()` → `BalanceResult`
- `pca(n_components, block=None)` → `PCAResult` (`.block_shares`: which blocks drive each PC)
- `variance(...)` → `VarianceResult`
- `baseline(...)` → `BaselineResult`
- `within_batch(base, ...)` → table
- `shap(base, target)` → `ShapResult`

## Results and output

- Every result has a `table`, and shows as it in a notebook
- `vf.save(folder, **results)`: tables as TSV

## Plots

- `vf.plots.qc`, `balance`, `pca`, `pca_by_batch`, `variance`, `baseline`, `roc_pr`, `within_batch`,
  `shap`, `shap_beeswarm`

## Statistics

- `vairrify.stat_utils`: `prepare_design`, `permanova`, `pvca`, `dispersion`, `variation_partitioning`,
  `cramers_v`, `metadata_associations`
