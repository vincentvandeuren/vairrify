# Study design: variables and samples

`design = reps.design(label=..., batch=..., other_technical=[...], covariates=[...], controls=[...], exclude=qc.failed)`

## What it fixes

- The variables and their roles (see [Concepts](../concepts.md))
- The repertoires every later step uses

## Settings

- `label`, `batch`: categorical; cast with `with_metadata` if stored as numbers
- `covariate_missing`: `"level"` (an `"unknown"` group) or `"drop"`
- `min_level_size`: covariate and control levels smaller than this merge into `"other"`

## Reading the result

- `design.summary()`: each variable's role and kind, and how PVCA treats it
- `design.samples()`: how many repertoires and patients each step uses

## Caveats

- TODO: nested technical variables (e.g. run within experiment)
