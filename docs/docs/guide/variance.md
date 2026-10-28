# Sources of variation (Fig. 2A)

`var = design.variance(permutations=999, pvca_bins={"age": [0, 10, 20, ...]})` · `vf.plots.variance(var)`

## What it answers

How much of the variation in the features does each variable explain: the batch, the label, the
covariates?

## How it works

- **PVCA**: variance components of a mixed model fitted to the leading PCs, as in
  `pvca::pvcaBatchAssess`. Continuous variables are binned (`pvca_bins`, or `pvca_n_bins` equal groups).
- **PERMANOVA**, as `vegan::adonis2`: each variable's share of the total variation (R²), **alone** and
  **adjusted** for all others, with a permutation p-value

## Reading the result

- A label that is large alone but small adjusted shows up mainly through what it is tied to, usually
  the batch
- Controls should sit near zero
- PVCA shares and PERMANOVA R² answer different questions and don't add up
- Zaslavsky (TRB, V genes): the batch takes 64% of the variance (PVCA) and 39% of the total variation
  alone (PERMANOVA R²), 12% adjusted. With depth (`log10_clonotypes`) as a technical variable, the
  batch's PVCA share is 52%.

## In the manuscript (Fig. 2A)

- PVCA in R (the `pvca` package), with age binned in decades (0–9, 10–19, ..., 70+)
- The study source (batch) was the largest contributor among all available metadata: 14% of the
  variance in BCR V-gene usage, 53% in TCR V-gene usage
- This despite every sample being PBMC and processed with the same library preparation protocol

## Differences from the manuscript

- vAIRRify's PVCA is a Python port of `pvca::pvcaBatchAssess`, without its `vsn2` normalisation
- PERMANOVA isn't in the manuscript
- TODO: check which variables Fig. 2A included (with depth, zaslavsky's batch share is 52%, close to
  the manuscript's 53%)

## Caveats

- PVCA needs at least two features
- PERMANOVA's permutations treat a patient's repertoires as independent
- With depth features, `log10_clonotypes` as a variable is circular
