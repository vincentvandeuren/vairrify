# Within-batch evaluation and SHAP (Fig. 3B)

`within = design.within_batch(base, min_group_class_size=5)` · `vf.plots.within_batch(within)`
`shap = design.shap(base, target=None, top=15)` · `vf.plots.shap(shap)` ·
`vf.plots.shap_beeswarm(shap, by="value")` · `vf.plots.shap_beeswarm(shap, by="batch")`

## What it answers

- Within batch: does the signal hold once the batch is fixed?
- SHAP: which features carry it?

## How it works

- **All batches**: the out-of-fold AUC from the baseline
- **Within batch**: the same scores, judged inside mixed-label batches that hold at least
  `min_group_class_size` of both the class and the rest
- **SHAP**: mean |SHAP| over out-of-fold repertoires, for one label class (by default the most
  recoverable one)
- Bar colour: positively associated with the class (pink) or negatively (blue), from the sign of the
  correlation between a feature's value and its SHAP value
- **Beeswarms**: the repertoires behind the bars, for the same `top` features. One dot per out-of-fold
  repertoire and feature, at its SHAP value: how much that feature changed the repertoire's predicted
  probability of the class. Dots pile up where values are similar; their height within a row means
  nothing else
    - `by="value"`: coloured by the feature's value in that repertoire (pink high, blue low)
    - `by="batch"`: coloured by the repertoire's batch, triangles for repertoires of the class. With
      more batches than colours, the batches whose mean prediction is furthest from the average keep a
      colour
    - Last row: per repertoire, the summed SHAP values of all the features not shown (`others=False`
      hides it)

## Reading the result

- A high AUC across batches that drops within batches points to a batch shortcut
- Blank cells: no batch holds both the class and the rest
- Zaslavsky (TRB, V genes): HIV reaches AUC 0.96 across batches but 0.65 within the HIV study
- Batch beeswarm: one colour on the positive side of most rows points to a batch that carries the
  signal. Zaslavsky (TRB, V genes): the repertoires pushed towards influenza are the Flu vaccine UPenn
  2021 batch. Whether that is influenza or the batch, this plot can't tell; the within-batch AUC can
- With many features (e.g. VJ bias), the summed row can outweigh the features shown: the model then
  relies on many small differences rather than a few features

## In the manuscript (Fig. 4)

- The idea: a classifier that uses batch-associated signal as a shortcut looks better across batches
  than it is; a within-batch evaluation is invariant to such shortcuts, and closer to a clinical
  setting where disease and batch aren't aligned
- The manuscript applied it to the published Mal-ID models; vAIRRify applies it to the baseline classifier
- For HIV, T1D and lupus, and for BCR and TCR alike, within-batch AUCs were lower than across-batch ones:
    - lupus (BCR + TCR metamodel): 0.75 and 0.59 within the two mixed-label batches, against 0.98
      across batches
    - HIV (TCR metamodel): 0.90 across batches, 0.61 within the mixed-label HIV study
- Healthy repertoires from mixed-label batches scored higher for the disease than healthy repertoires
  from healthy-only batches

## SHAP in the manuscript (Fig. 3B, S2, S3)

- `TreeExplainer` (shap 0.48.0) on out-of-fold predictions
- Influenza: repertoire size first; TRBV18 enriched and TRBV3-1 depleted relative to other datasets,
  "likely due to differences in sequencing library preparation"
- In vAIRRify's zaslavsky run, TRBV3-1 comes out *positively* associated with influenza: after tcrio
  ingestion it's present in all 37 influenza repertoires but in only 223 of 542 overall (e.g. 11 of 95
  HIV, 13 of 58 COVID-19). TODO: check the ingestion.

## Caveats

- Few batches are mixed-label in many merged datasets; with small classes, lower `min_group_class_size`
