# Concepts

## Shortcut learning in AIRR-seq

- Technical factors (sequencing run, library prep, depth, cell input) shape which receptors are
  observed, and how many
- When those factors line up with a label, a model can predict the label from them: a shortcut that
  won't transfer to a clinic

### In the manuscript (Fig. 1)

- The observed repertoire is shaped by the disease-specific effect plus three factors: genetic
  background and exposure history (biological), and technical factors (which receptors, and how many,
  are sequenced)
- Technical sources it names:
    - primer hybridisation and amplification efficiency, which distort V-gene usage and clone
      frequencies (and vary with GC content)
    - transcript stability in RNA-based protocols
    - the number of cells recovered and the sequencing depth, which affect clonal abundance, diversity
      and the detection of rare V genes
    - sample preparation (e.g. the technician), and cell composition and demographics differing between batches
- Biological shortcuts:
    - HLA genotype: for diseases with HLA-associated risk (e.g. T1D) a model can learn risk alleles
      instead of disease-associated receptors
    - infections with an uneven geography (HIV, CMV), together with geographic HLA frequencies: a
      model can recognise where a subject comes from rather than the disease
- Shortcuts in AIRR data are harder to spot than in imaging: a feature like V-gene frequency can be
  shaped both by confounding and by the disease, and the disease-associated receptors that would serve
  as ground truth are usually unknown

## Repertoire-shaping factors and variable roles

Every metadata variable used gets one role:

| Role | What it is | Example |
|---|---|---|
| `label` | what a downstream model would predict (categorical) | disease, serostatus |
| `batch` | the main technical grouping | study, run, library prep |
| `other_technical` | further technical factors | experiment, depth (`log10_clonotypes`) |
| `covariates` | biological factors that may shape the repertoire | age, sex, ancestry |
| `controls` | should explain nothing | CV fold, computational split, `random_control` |

- A variable's kind follows its dtype: numbers are continuous, everything else categorical
- One design looks at one label; for a second label, make a second design

## Features

- The features come from one [feature block](feature-blocks.md), or several side by side
- One feature selection per run: QC, PCA, variance and the baseline classifier all use it
- With several blocks, each counts equally in QC, PCA and the variance analysis, however many
  features it has

## Samples

- Analysed: repertoires with every feature, a label and a batch, not excluded by QC
- Repeated measures: a patient's repertoires all stay. The baseline's random folds keep them
  together; PERMANOVA's permutations treat them as independent.
- PCA and the variance analysis also drop repertoires missing a continuous variable
- TODO: `covariate_missing`, `min_level_size`

## Recommendations from the manuscript

For study design:

- Process samples of the different labels in the same batches, so batch and disease aren't confounded
- Match age, sex and ancestry between the disease and control groups

For analysis (what vAIRRify automates):

- Project the data, e.g. a PCA of V-gene usage, to find the largest sources of variation and outliers
- Train a "basic" classifier on general repertoire properties to reveal biases missed by eye
- In datasets integrated from several studies, evaluate on subsets of consistently generated
  repertoires (within batch)

Beyond vAIRRify:

- Validate receptor-level models with well-characterised disease-associated receptors; their
  frequency also makes a useful baseline
