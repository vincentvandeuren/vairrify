# PCA (Fig. 2B)

`pca = design.pca(n_components=10)` · `vf.plots.pca(pca, x="PC1", y="PC2")` · `vf.plots.pca_by_batch(pca)`
· per block: `design.pca(n_components=10, block="V genes")`

## What it answers

Do repertoires cluster by batch rather than by label?

## How it works

- Features standardised (for gene usage: log1p counts, L2-normalised, as in the manuscript)
- `pca.explained`: variance explained per component; `x` and `y` pick the components to plot
- Loadings panel: the 12 features loading most on the chosen components, and in orange the V genes
  whose usage varies for germline reasons (see [Germline variation](../germline.md))
- Per-batch panels: one batch coloured over all others in grey
- Several blocks: each block's standardised features are weighted by 1 / √(its number of features that
  vary in the analysed repertoires), so every block has the same total variance. Features that don't
  vary (e.g. HLA alleles nobody in the cohort carries) don't count. That doesn't make them count equally on the leading PCs: a
  block whose features vary together puts its variance on one axis and can take a whole PC, while a block
  of many independent features spreads it over many and hardly shows
    - The axis labels name the blocks driving each PC, e.g. "PC1 (27%) – 90% depth": each block's
      share of the PC's squared loadings. `pca.block_shares` has them for every PC
    - A third panel splits each PC's explained variance by block (a bar per PC, the plotted two in
      bold): a segment is the variance that PC captures from that block. Each block starts with the same
      share of the total, so the legend gives how much of each block's own variance the PCs shown
      capture. Zaslavsky (V genes + VJ bias + depth): PC1-10 hold all of depth, 69% of V genes, 9% of
      VJ bias
    - Block colours are muted, so they don't read as the label colours of the PCA panel beside them
    - `design.pca(..., block=b)`: a PCA of one block alone, as in the manuscript's V-gene PCA

## Reading the result

- Batches forming their own clusters means batch drives the features
- Features that load on these PCs and also carry the baseline classifier (SHAP) point to a batch shortcut
- Zaslavsky (TRB, V genes): PC1 explains 24% and PC2 14%. TODO: describe the clusters.
- Zaslavsky (TRB, V genes + VJ bias + depth): joint PC1 is 90% depth, whose five features move largely
  together (73% on its own PC1); VJ bias makes up at most 7% of the first four joint PCs, so only its
  own PCA shows it

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
- With several blocks, read the joint PCA together with the per-block ones: a block absent from the
  leading joint PCs can still separate batches on its own
