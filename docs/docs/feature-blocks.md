# Feature blocks

A feature block is a per-repertoire table computed by tcrio. `vf.feature_blocks()` lists them.
*Technical* blocks are the manuscript's general repertoire properties, which technical bias shapes;
*biological* blocks are confounders.

| Block | Kind | Source (tcrio) | Features |
|---|---|---|---|
| V genes | technical | `GeneCountsSummary` | share of clonotypes using each V gene |
| J genes | technical | `GeneCountsSummary` | share of clonotypes using each J gene |
| VJ bias | technical | `GeneCountsSummary` | log ratio of each V-J pairing to independent V and J usage |
| VDJ stats | technical | `VdjStatisticsSummary` | mean, sd, skew, kurtosis of trimming, P and N nucleotides |
| repertoire size | technical | repertoire counts | number of unique clonotypes |
| depth | technical | `RarefactionReport` | reads; expected clonotypes at 20k / 48k / 200k / 500k reads |
| HLA | biological | `RepertoireHlaInference` | inferred HLA carriage (score > 0.5) |

## Notes per block

- **V genes, J genes**: PCA and the variance analysis use the manuscript's profile (log1p counts,
  L2-normalised). TRBV6-2 and TRBV6-3 are identical; tcrio reports both as TRBV6-2.
- **V genes with germline variation** (TRBV4-3, TRBV6-2/6-3, TRBV7-2, TRBV7-3, TRBV28, TRBV30):
  deletion polymorphisms and linked alleles, so their variation needn't be technical. Marked in the PCA
  plot; evidence on [Germline variation](germline.md).
- **depth**: past a repertoire's own reads the rarefaction curve is extrapolated. TODO: how to read
  that, and the share of repertoires reaching each depth.
- **HLA**: dichotomised at 0.5, as in the manuscript, so repertoire size doesn't leak in. On the
  zaslavsky data the dichotomised scores predict much less than Table 1's HLA column (influenza
  AUC 0.54 vs 0.96); the raw scores reproduce it. TODO: check which the table used.
- Only TRB is supported by `GeneCountsSummary` and the HLA model for now

## In the manuscript (Methods, "Baseline machine learning")

- **Repertoire size**: the total number of unique clonotypes
- **V and J usage**: the frequency of clonotypes using each gene
- **V-J pairing**: the log-fold change of each V-J pair relative to what independent V and J usage predicts
- **VDJ statistics**: germline-encoded nucleotides identified from the V and J alleles (for D, the
  longest possible match); from these, exonuclease deletions (trims) and untemplated (N) insertions at
  each junction boundary and their totals; then the mean, standard deviation, skewness and kurtosis
- **HLA**: inferred from the TCR repertoire with the method of ElAbd et al. (HLA2TCR), dichotomised
  (> 0.5) "to prevent repertoire size from influencing HLA genotype prediction scores"
- Duplicate clonotypes (same `v_call`, CDR3 amino acids and `j_call`) were merged; V-gene counts use a
  single clonal representative
- The general repertoire properties came to 676 features

## Adding a block

- TODO: example of a `FeatureBlock` reading another tcrio result, registered in `FEATURE_BLOCKS`
