# Walkthrough notebooks

The `notebooks/` folder of the repository walks through a full audit, step by step. Each decision sits
in a cell marked ⚙.

| Notebook | What it does | Needs |
|---|---|---|
| `0_ingest_data` | clonotype files and metadata → a tcrio dataset; filtering report | clonotype tables, metadata |
| `1_sample_overlap_and_contamination` | shared clonotypes between repertoires, possible contamination | an ingested dataset |
| `2f_audit` (TODO: rename) | the audit: QC, design, balance, PCA, variance, baseline, within batch, SHAP | an ingested dataset |

- TODO: an example dataset readers can download, and the expected run time
- TODO: render the notebooks here with mkdocs-jupyter
