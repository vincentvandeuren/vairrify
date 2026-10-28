# Limitations

- **Repeated measures**: a patient's repertoires all stay. The baseline's random folds keep them
  together, but PERMANOVA's permutations treat them as independent. TODO: restricted permutations.
- **Labels**: one categorical label per design. Multi-class labels are predicted one class against
  the rest; continuous labels aren't supported (bin them with `with_metadata`).
- **Mixed-label batches**: the within-batch evaluation needs batches holding both a class and the
  rest. Merged datasets often have few.
- **Extrapolated depth**: clonotypes at depths beyond a repertoire's own reads come from tcrio's
  extrapolation, not from measured data
- **Loci**: V/J gene usage and HLA inference are TRB-only in tcrio for now
- **HLA**: inferred from the repertoire, so the scores may carry dataset biases themselves
- **Technical vs biological**: some V genes vary for germline reasons (see
  [Germline variation](germline.md)); a predictive feature isn't technical by definition

## Limitations the manuscript discusses

- Treating broad repertoire differences between diseases (e.g. V-gene usage) as technical is an
  assumption: they can carry disease signal too, e.g. TRBV9+ T cells likely underlie ankylosing spondylitis
- Some features labelled technical may be biological bystanders, such as HLA associations of a
  disease; such features should still predict within batches
- Inferred HLA scores come from a model trained on another dataset; they may reflect dataset biases,
  and repertoire properties besides HLA
- DNA- and RNA-based protocols sample different molecules, so a failure to transfer between them can be
  biological as well as technical
- Independent cohorts differ in depth and patient characteristics (e.g. mild vs hospitalised COVID-19),
  which complicates any transfer
