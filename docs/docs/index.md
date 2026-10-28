# vAIRRify

Machine-learning models trained on adaptive immune receptor repertoires (AIRR-seq) can reach high
accuracy by recognising technical variation, such as the sequencing batch, depth or library
preparation, instead of the biology they are meant to learn: **shortcut learning**. vAIRRify audits a
dataset for such shortcuts, before a model is trained or while its results are evaluated.

It implements the analyses of the first part of *Disentangling Biological Signal from Technical Bias
in Disease Diagnostics from Adaptive Immune Receptor Repertoires* (see [Citation](#citation)).

Why it matters, from the manuscript:

- In the Mal-ID training data (Zaslavsky et al.), a random forest on simple repertoire properties
  (repertoire size, V and J usage, V-J pairing, VDJ recombination statistics) predicted each disease
  as well as, or better than, three published classifiers (Table 1)
- Those classifiers performed worse within batches than across them, a sign of shortcut learning (Fig. 4)
- They also transferred poorly to independent COVID-19, lupus and T1D cohorts

## What it checks

| Question | Analysis | Manuscript |
|---|---|---|
| Are the label and the batch tied together? | repertoires per batch and label, Cramér's V, mixed- vs single-label batches | Fig. S1 |
| What drives variation in the repertoire properties? | PCA, PVCA, PERMANOVA | Fig. 2 |
| Can simple repertoire properties predict the label? | random-forest baseline classifier | Table 1, Fig. 3A |
| Does that signal hold once the batch is fixed? | within-batch evaluation | within-batch evaluation |
| Which features carry it? | SHAP | Fig. 3B |

## Install

```bash
pip install vairrify          # TODO: once on PyPI; until then see Installation
```

## Quickstart

```python
import vairrify as vf

reps = vf.load("data/processed/my_dataset", locus="TRB", features="V genes")
qc = reps.qc()
design = reps.design(label="disease", batch="study_name", covariates=["age", "sex"], exclude=qc.failed)

design.variance()                    # what drives V-gene usage
base = design.baseline()             # can V-gene usage predict the label?
design.within_batch(base)            # ... also within mixed-label batches?
```

## Where next

- [Installation](installation.md) and [Getting started](getting-started.md): from clonotype files to a first audit
- [Concepts](concepts.md): shortcut learning, variable roles, feature blocks
- [The audit, step by step](guide/qc.md): what each step answers and how to read it
- [Walkthrough notebooks](walkthrough.md)

## Citation

Stadelmaier J.\*, Van Deuren V.\*, Greiff V., Kowarik M., Mariotti-Ferrandiz E., Minervina A. A.,
Nahnsen S., Pogorelyy M. V., Rawat P., Sandve G. K., Semerci A., Stervbo U., Thomas P. G., Yaari G.,
Gabernet G.§, Meysman P.§ *Disentangling Biological Signal from Technical Bias in Disease Diagnostics
from Adaptive Immune Receptor Repertoires.* TODO: journal or preprint, year, DOI.
(\* shared first authors, § shared senior authors)

vAIRRify is MIT-licensed.
