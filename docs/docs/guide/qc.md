# QC: repertoires that shouldn't be analysed

`qc = reps.qc(min_clonotypes=1000, max_zero=10, outlier_iqr=3)` · `vf.plots.qc(qc)`

## What it answers

Which repertoires look like failed libraries? A few of them can dominate everything downstream.

## The checks

| Flag | Rule | Failure? |
|---|---|---|
| low depth | fewer than `min_clonotypes` clonotypes | yes |
| missing genes | more than `max_zero` genes never observed (gene-usage blocks only) | yes |
| atypical profile | distance to the typical profile more than `outlier_iqr` IQRs above the upper quartile | no, a reason to look |

## Reading the result

- `qc.failed`, `qc.atypical`, `qc.flagged`: repertoire IDs, to pass as `exclude=` to the design
- The plot: depth against distance to the typical profile, with both thresholds
- TODO: example of an atypical repertoire that is biology (e.g. lymphopenia)

## Caveats

- Low depth and atypical profiles can be biology; check where exclusions fall (per batch, per label)
