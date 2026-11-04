# Walkthrough notebooks

These notebooks run a full audit, from clonotype files to the shortcut analysis. Read them here, with
their outputs, or download them to run on your own data.

| Notebook | Download | What it does | Needs |
|---|---|---|---|
| [0. Ingest the data](nbs/0_ingest_data.ipynb) | [⬇ .ipynb](download/0_ingest_data.ipynb){: download="0_ingest_data.ipynb" } | standardization and filtering of clonotype tables and metedata | clonotype tables, metadata |
| [1. Sample overlap and contamination](nbs/1_sample_overlap_and_contamination.ipynb) | [⬇ .ipynb](download/1_sample_overlap_and_contamination.ipynb){: download="1_sample_overlap_and_contamination.ipynb" } | shared clonotypes between repertoires, possible contamination | an ingested dataset |
| [2. The vAIRRify audit](nbs/2_shortcut_learning.ipynb) | [⬇ .ipynb](download/2_shortcut_learning.ipynb){: download="2_shortcut_learning.ipynb" } | shortcut learning analysis on general repertoire properties | an ingested dataset |

## Advanced examples

The audit of notebook 2 with other choices. Start with the walkthrough above.

| Notebook | Download | What changes |
|---|---|---|
| [B cells (IGH)](nbs/2_shortcut_learning_bcr.ipynb) | [⬇ .ipynb](download/2_shortcut_learning_bcr.ipynb){: download="2_shortcut_learning_bcr.ipynb" } | the same dataset, on the B-cell heavy chain instead of TRB, with V and J genes |
| [Another dataset and label](nbs/2_shortcut_learning_vlasova.ipynb) | [⬇ .ipynb](download/2_shortcut_learning_vlasova.ipynb){: download="2_shortcut_learning_vlasova.ipynb" } | the `vlasova_2026` dataset, with COVID status as the label and V genes, J genes and depth combined |

## Running a notebook

1. Install vairrify (see [Installation](installation.md)) and JupyterLab: `pip install jupyterlab`.
2. Download the notebook and open it with `jupyter lab`.
3. Point them to your own datasets and choose your parameters
