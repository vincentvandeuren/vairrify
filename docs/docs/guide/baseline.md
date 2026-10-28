# Baseline classifier (Table 1, Fig. 3A)

`base = design.baseline(cv_column=..., min_class_size=10, cv_repeats=3)` · `vf.plots.baseline(base)`

## What it answers

Can a simple classifier predict the label from the features alone? If it can, a downstream model
can take the same route.

## How it works

- A random forest predicts each label class against the rest, scored out of fold (AUC and average
  precision)
- `cv_column`: the downstream model's own folds. Without it: 3 random folds that keep patients
  together and are stratified on batch × label, repeated `cv_repeats` times.
- The batch is predicted too (a positive control), as are the controls (should sit at 0.5)

## Reading the result

- A high AUC is a shortcut risk, not evidence of biology
- TODO: Table 1 comparison for zaslavsky, run with vAIRRify (general repertoire properties, HLA, both)

## In the manuscript (Table 1, Fig. 3)

- Random forest, 200 trees, `class_weight="balanced_subsample"`, no tuning (scikit-learn 1.8.0);
  the same 3-fold cross-validation splits as the original Mal-ID study; ROC AUC and average precision
  on the out-of-fold scores
- TCR only, which showed the stronger biases; 676 features of general repertoire properties
- Out-of-fold AUC per disease (one against the rest):

| Disease | General properties | HLA | General + HLA | Mal-ID models 1 / 2 / 3 (reported) |
|---|---|---|---|---|
| HIV | 0.964 | 0.894 | 0.963 | 0.93 / 0.74 / 0.93 |
| T1D | 0.974 | 0.960 | 0.985 | 0.93 / 0.72 / 0.94 |
| Lupus | 0.949 | 0.732 | 0.955 | 0.92 / 0.81 / 0.92 |
| COVID-19 | 0.958 | 0.884 | 0.958 | 0.91 / 0.89 / 0.91 |
| Influenza | 0.999 | 0.960 | 1.000 | 0.97 / 0.83 / 0.95 |
| Healthy | 0.957 | 0.826 | 0.959 | 0.91 / 0.80 / 0.90 |

- Influenza: AUC 0.999, AP 0.987; repertoire size was the most important feature. Influenza
  repertoires were significantly larger than the rest (Mann-Whitney U, p = 9.68e-20), reflecting
  sequencing depth or cell numbers in that study.
- HLA features alone reached AUC 0.960 for influenza
- HIV from HLA: the most predictive alleles were tied to ancestry (European-associated alleles
  negatively, HLA-A\*23:01, enriched in African populations, positively), none with known HIV
  associations: a geographic shortcut

## Caveats

- Forest settings differ from the manuscript's (400 trees, `max_features=0.05`, `max_depth=10`; pass
  `rf=` to change them)
- Classes with fewer than `min_class_size` repertoires (or rest) are skipped
