"""Steps 5-11: the study design, a fixed choice of variables and samples, and the analyses on it.

The variables are the manuscript's repertoire-shaping factors (Fig. 1): technical factors (batch,
depth, ...), the label a model would predict (disease status), biological covariates (age, sex,
ancestry, ...), and controls that should explain nothing. A variable's kind follows its dtype:
numbers are continuous, everything else categorical.

Samples: repertoires with every feature and a label and batch value, not excluded by QC. A
patient's repertoires all stay: the baseline classifier's random folds keep them together, but
PERMANOVA's permutations treat them as independent. PCA and variance also drop repertoires missing
a continuous variable, and with `covariate_missing="drop"` those missing a categorical covariate.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Mapping, Optional, Sequence

import numpy as np
import polars as pl
from scipy import stats

from . import baseline as _baseline
from . import stat_utils as su
from . import variance as _variance
from .repertoires import Repertoires, kind
from .results import Result


@dataclass(frozen=True, repr=False)
class BalanceResult(Result):
    """`table`: repertoires per batch and label level, and whether each batch holds mixed labels."""
    batch: str
    label: str
    cramers_v: float
    chi2_p: float
    design: str          # balanced / imbalanced / confounded


@dataclass(frozen=True)
class StudyDesign:
    repertoires: Repertoires
    roles: Mapping[str, str]             # variable -> technical / label / covariate / control
    kinds: Mapping[str, str]             # variable -> categorical / continuous
    label: str
    batch: Optional[str]
    metadata: pl.DataFrame               # all repertoires; rare covariate / control levels merged into "other"
    keep: np.ndarray                     # the analysed repertoires
    excluded: np.ndarray                 # excluded by QC
    structure: su.Design = field(repr=False)  # PCA / variance: their profiles (standardised, blocks weighted) and metadata

    # --- 5. variables / 6. samples ----------------------------------------------------------

    @classmethod
    def build(cls, reps: Repertoires, *, label, batch, other_technical, covariates, controls, exclude,
              covariate_missing, min_level_size) -> "StudyDesign":
        md = reps.metadata
        roles = {}
        for role, group in [("technical", ([batch] if batch else []) + other_technical),
                            ("label", [label]), ("covariate", covariates), ("control", controls + ["random_control"])]:
            for v in group:
                if v in roles:
                    raise ValueError(f"{v!r} is listed as both {roles[v]} and {role}")
                if v not in md.columns:
                    raise ValueError(f"{v!r} is not a metadata column; see overview()")
                roles[v] = role
        kinds = {v: kind(md[v]) for v in roles}
        for v in [label, batch]:
            if v and kinds[v] != "categorical":
                raise ValueError(f"{v!r} must be categorical; cast it with with_metadata(pl.col({v!r}).cast(pl.String))")

        excluded = np.isin(reps.ids, list(exclude)) & reps.present
        keep = reps.present & ~excluded & md[label].is_not_null().to_numpy()
        if batch:
            keep &= md[batch].is_not_null().to_numpy()

        # merge rare covariate / control levels into "other"
        analysis_md = md
        for v in [c for c in roles if kinds[c] == "categorical" and roles[c] in ("covariate", "control")]:
            counts = md.filter(pl.Series(keep)).group_by(v).len()
            rare = counts.filter((pl.col("len") < min_level_size) & pl.col(v).is_not_null())[v].to_list()
            if rare:
                print(f"{v}: merging {len(rare)} level(s) with < {min_level_size} repertoires into 'other': {rare}")
                analysis_md = analysis_md.with_columns(
                    pl.when(pl.col(v).is_in(rare)).then(pl.lit("other")).otherwise(pl.col(v).cast(pl.String)).alias(v))

        structure, dropped = _variance.prepare(reps.profile()[keep], analysis_md.filter(pl.Series(keep)),
                                               [v for v in roles if kinds[v] == "continuous"],
                                               [v for v in roles if kinds[v] == "categorical"], covariate_missing)
        for v in dropped:
            print(f"{v!r} left out of PCA / variance: a single level left")
        structure = structure.with_X(su.standardize(structure.X) * reps.weights(structure.X))
        return cls(reps, roles, kinds, label, batch, analysis_md, keep, excluded, structure)

    @property
    def variables(self) -> list[str]:
        """The variables in the variance analysis (all of them, unless one was left with a single level)."""
        return [v for v in self.roles if v in self.structure.variables]

    def summary(self) -> pl.DataFrame:
        """Each variable's role and kind, and how PVCA treats it."""
        md = self.repertoires.metadata
        return pl.DataFrame({
            "variable": list(self.roles),
            "role": list(self.roles.values()),
            "kind": [self.kinds[v] for v in self.roles],
            "distinct": [md[v].drop_nulls().n_unique() for v in self.roles],
            "nulls": [md[v].null_count() for v in self.roles],
            "PVCA": ["binned" if self.kinds[v] == "continuous" else "as is" for v in self.roles],
        })

    def samples(self) -> pl.DataFrame:
        """How many repertoires each analysis uses, and where the rest went."""
        reps = self.repertoires
        present = reps.present
        return pl.DataFrame([
            ("all repertoires", reps.metadata.height),
            (f"with {reps.locus} data for every feature", int(present.sum())),
            ("not excluded by QC", int((present & ~self.excluded).sum())),
            (f"analysed: has {self.label!r}" + (f" and {self.batch!r}" if self.batch else ""), int(self.keep.sum())),
            ("  of which patients", self.metadata.filter(pl.Series(self.keep))["patient_id"].n_unique()),
            ("PCA / variance: no missing continuous values", self.structure.n),
        ], schema=["step", "repertoires"], orient="row")

    # --- 7. balance -------------------------------------------------------------------------

    def balance(self) -> BalanceResult:
        """Repertoires per batch and label level (manuscript Fig. S1), and how tied label and batch are
        (bias-corrected Cramér's V): balanced < 0.1 <= imbalanced < 0.5 <= confounded. Only mixed-label
        batches, holding more than one label level, allow a within-batch evaluation."""
        if not self.batch:
            raise ValueError("no batch set")
        md = self.metadata.filter(pl.Series(self.keep))
        wide = md.group_by(self.batch, self.label).len().pivot(on=self.label, index=self.batch, values="len").fill_null(0)
        levels = sorted(c for c in wide.columns if c != self.batch)
        counts = wide.select(levels).to_numpy()
        table = (wide.select(self.batch, *levels)
                 .with_columns(pl.Series("batch type", np.where((counts > 0).sum(axis=1) > 1, "mixed-label", "single-label")))
                 .sort("batch type", self.batch))
        v = su.cramers_v(md[self.batch].to_numpy(), md[self.label].to_numpy())
        p = stats.chi2_contingency(counts)[1] if min(counts.shape) > 1 else np.nan
        return BalanceResult(table, self.batch, self.label, v, p,
                             "balanced" if v < 0.1 else ("imbalanced" if v < 0.5 else "confounded"))

    # --- 8-9. structure ---------------------------------------------------------------------

    def pca(self, n_components: int = 10, block: Optional[str] = None):
        """Manuscript Fig. 2B. By default of all loaded blocks together, each counting equally; `block`
        restricts it to one of them. `.block_shares` shows which blocks drive each PC."""
        return _variance.pca(self, n_components, block)

    def variance(self, *, permutations: int = 999, pvca_bins: Optional[Mapping[str, Sequence[float]]] = None,
                 pvca_n_bins: int = 4):
        """Manuscript Fig. 2A (PVCA, with at least two features) and PERMANOVA."""
        return _variance.variance(self, permutations=permutations, pvca_bins=pvca_bins or {}, pvca_n_bins=pvca_n_bins)

    # --- 10-11. baseline classifier, within batch -------------------------------------------

    def baseline(self, *, cv_column: Optional[str] = None, rf=None, min_class_size: int = 10, cv_repeats: int = 3):
        """Manuscript Table 1."""
        return _baseline.baseline(self, cv_column=cv_column, rf=rf, min_class_size=min_class_size, cv_repeats=cv_repeats)

    def within_batch(self, base, *, min_group_class_size: int = 5) -> pl.DataFrame:
        return _baseline.within_batch(self, base, min_group_class_size=min_group_class_size)

    def shap(self, base, target: Optional[str] = None, top: int = 15):
        """Mean |SHAP| for one label class (`target`); by default the most recoverable one. Draw it with
        `plots.shap`, and the per-repertoire values with `plots.shap_beeswarm`."""
        return _baseline.shap_importance(self, base, target=target, top=top)
