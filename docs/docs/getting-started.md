# Getting started

From clonotype files to a first audit, following the [walkthrough notebooks](walkthrough.md).

## 1. Ingest the data (notebook 0)

- What you need: a folder of clonotype tables (AIRR, MiXCR, Adaptive, 10x, ...), and a way to map each
  file to a repertoire and a patient (e.g. a regular expression on the file name)
- `tcrio`'s `DatasetIngester` writes a dataset in a common format
- Metadata: patient-level (age, sex, ancestry) and repertoire-level (batch, disease status, CV fold)
- Check the filtering report: which rows were dropped, and why
- TODO: link to tcrio's documentation for readers and mappers

## 2. Check sample overlap and contamination (notebook 1, optional)

- Shared clonotypes between pairs of repertoires; unrelated people rarely share many of their top clones
- In the manuscript (Rawat et al. T1D cohort): a pair whose top 250 clones shared ≥ 50 clones, with a
  positive Spearman correlation of their duplicate counts, was treated as contamination and removed
  (6 pairs, 12 patients)
- TODO: tcrio's `OverlapAnalyzer` uses `top_n=250` but flags at 40 shared clones by default; align with
  the manuscript's 50, or explain the difference

## 3. Run the audit (notebook 2)

- Choose the dataset, the locus and the features (see [Feature blocks](feature-blocks.md))
- Inspect and clean the metadata
- QC, then the study design: label, batch, covariates, controls
- Balance, PCA, sources of variation, baseline classifier, within-batch evaluation, SHAP
- To compare features (e.g. V genes with depth), rerun with another choice

## 4. Outputs

- Every step returns a table, and `vf.plots` draws it
- `vf.save(folder, qc=qc, ...)` writes the tables as TSV
- TODO: a suggested folder layout per dataset, locus, label and features
