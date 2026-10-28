# PCA (Fig. 2B)

`pca = design.pca(n_components=10)` · `vf.plots.pca(pca, x="PC1", y="PC2")` · `vf.plots.pca_by_batch(pca)`

## What it answers

Do repertoires cluster by batch rather than by label?

## How it works

- Features standardised (for gene usage: log1p counts, L2-normalised, as in the manuscript)
- `pca.explained`: variance explained per component; `x` and `y` pick the components to plot
- Loadings panel: the 12 features loading most on the chosen components, and in orange the V genes
  whose usage varies for germline reasons (see [Germline variation](../germline.md))
- Per-batch panels: one batch coloured over all others in grey

## Reading the result

- Batches forming their own clusters means batch drives the features
- Features that load on these PCs and also carry the baseline classifier (SHAP) point to a batch shortcut
- Zaslavsky (TRB, V genes): PC1 explains 24% and PC2 14%. TODO: describe the clusters.

## In the manuscript (Fig. 2B)

- V-gene counts were log1p-transformed, L2-normalised and scaled, then PCA (scikit-learn), for TCR beta
  and BCR heavy chain
- Several V genes drive the batch-related variation, consistent with technical differences in
  amplification and sequencing efficiency
- TCR repertoires from paediatric cohorts clustered at high PC1. Age and batch were confounded, so
  their contributions couldn't be separated.
- In batches holding both healthy and diseased individuals, such as the HIV study, repertoires
  clustered together regardless of status: the skew elsewhere is likely batch-dependent, not disease-dependent
- Recommendation: run such a projection on every dataset, to catch outliers and clear biases

## Caveats

- Needs at least two features
