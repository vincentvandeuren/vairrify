# Balance over batches (Fig. S1)

`balance = design.balance()` · `vf.plots.balance(balance)`

## What it answers

Is the label tied to the batch? If so, recognising the batch is enough to predict the label.

## How it works

- Repertoires per batch and label level
- Bias-corrected Cramér's V between batch and label: **balanced** (< 0.1), **imbalanced** (0.1–0.5),
  **confounded** (≥ 0.5)
- Each batch is **mixed-label** (more than one label level) or **single-label**

## Reading the result

- Only mixed-label batches allow a [within-batch evaluation](within-batch.md)
- Zaslavsky (TRB, disease): confounded (Cramér's V 0.93); 4 of 11 batches are mixed-label: the HIV
  study, the T1D biobank and the two adult lupus batches

## In the manuscript (Fig. S1)

- Repertoires of each disease class per batch, in the Mal-ID data
- Identifiers for library preparation or sequencing runs weren't available, so each study was
  treated as one batch
- Only HIV, T1D and lupus had mixed-label batches, so only those could be evaluated within batch;
  COVID-19, influenza and healthy controls came from single-label batches
- The T1D mixed-label batch held only three healthy repertoires
- Recommendation: process samples of the different labels in the same batches
