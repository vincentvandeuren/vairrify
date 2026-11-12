"""Steps 1-4: load a tcrio dataset, look at its metadata, clean it up, flag failed repertoires.

`Repertoires` holds one row per repertoire: its metadata, its depth, and the features every
analysis uses, from one feature block or several side by side (e.g. the manuscript's general
repertoire properties). It never changes once made: `with_metadata` returns a new one, so
re-running a notebook cell can't leave stale state behind. `design()` fixes the variables and the
samples and hands over to `StudyDesign` for the analyses.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import polars as pl
from tcr_io.dataset import Dataset

from .features import FEATURE_BLOCKS
from .qc import QCResult, flag_repertoires


def kind(s: pl.Series) -> str:
    # Numbers are continuous, everything else categorical; cast in with_metadata to change it.
    return "continuous" if s.dtype.is_numeric() else "categorical"


@dataclass(frozen=True)
class Repertoires:
    name: str
    locus: str
    features: list[str]                     # the feature blocks every analysis uses
    metadata: pl.DataFrame                  # one row per repertoire, in dataset order
    clonotypes: np.ndarray                  # unique clonotypes per repertoire; NaN without this locus
    X: np.ndarray                           # repertoires x features, the blocks side by side
    feature_names: list[str]
    block_of: np.ndarray                    # each feature's block
    dataset: Dataset = field(repr=False)
    seed: int = 0

    # --- 1. load ----------------------------------------------------------------------------

    @classmethod
    def load(cls, path: str | Path, locus: str, features: str | Sequence[str], *,
             seed: int = 0, run_operations: bool = True) -> "Repertoires":
        """Metadata, depth and the feature blocks, one row per repertoire, running a block's tcrio
        operation where it is missing or out of date. Two metadata columns are added: `log10_clonotypes` (depth; list
        it in `other_technical` to use it) and `random_control` (three random groups, a negative control)."""
        features = [features] if isinstance(features, str) else list(features)
        if unknown := [b for b in features if b not in FEATURE_BLOCKS]:
            raise ValueError(f"unknown feature block(s) {unknown}; see vf.feature_blocks()")
        path = Path(path)
        ds = Dataset(path).select_locus(locus)
        rep_ids = ds.repertoire_meta.select("repertoire_id")
        metadata = (rep_ids.join(ds.full_repertoire_meta, on="repertoire_id", how="left")
                    .join(ds.full_patient_meta, on="patient_id", how="left")
                    .drop(pl.selectors.by_dtype(pl.List)))
        clonotypes = rep_ids.join(ds.repertoire_counts, on="repertoire_id", how="left")["n_clonotypes"].cast(pl.Float64).to_numpy()
        with np.errstate(divide="ignore", invalid="ignore"):
            metadata = metadata.with_columns(
                pl.Series("log10_clonotypes", np.log10(clonotypes)),
                pl.Series("random_control", np.random.default_rng(seed).choice(["r1", "r2", "r3"], metadata.height)),
            )
        tables = [rep_ids.join(FEATURE_BLOCKS[b].read(ds, run_operations=run_operations), on="repertoire_id", how="left")
                  .drop("repertoire_id") for b in features]
        reps = cls(path.name, locus, features, metadata, clonotypes, np.hstack([t.to_numpy() for t in tables]),
                   [c for t in tables for c in t.columns], np.repeat(features, [t.width for t in tables]), ds, seed)
        print(f"{reps.name} {locus}: {metadata.height} repertoires, {int(reps.present.sum())} with "
              f"{' + '.join(features)} ({len(reps.feature_names)} features), {metadata.width} metadata columns")
        return reps

    @property
    def ids(self) -> np.ndarray:
        return self.metadata["repertoire_id"].to_numpy()

    @property
    def present(self) -> np.ndarray:
        """Repertoires with every feature for this locus."""
        return np.isfinite(self.X).all(axis=1) & np.isfinite(self.clonotypes) & (self.clonotypes > 0)

    def profile(self) -> np.ndarray:
        """The features as QC, PCA and the variance analysis see them: each block's profile, e.g. the
        manuscript's V-gene profile (log1p counts, L2-normalised)."""
        parts = []
        for b in self.features:
            X, profile = self.X[:, self.block_of == b], FEATURE_BLOCKS[b].profile
            parts.append(profile(X, self.clonotypes) if profile else X)
        return np.hstack(parts)

    def weights(self, X: np.ndarray) -> np.ndarray:
        """Per feature, 1 / sqrt(how many features of its block vary in `X`, the repertoires analysed): applied to
        standardised features, every block counts equally. A feature that doesn't vary (e.g. an HLA allele nobody
        carries) adds nothing, so it mustn't count towards its block's size either; it gets weight 0."""
        varies = np.nanstd(X, axis=0) > 0
        counts = {b: varies[self.block_of == b].sum() for b in self.features}
        return np.array([1 / np.sqrt(counts[b]) if counts[b] and v else 0.0 for b, v in zip(self.block_of, varies)])

    # --- 2. inspect -------------------------------------------------------------------------

    def overview(self) -> pl.DataFrame:
        """Every metadata column: its kind, distinct and missing values, and a few examples."""
        md = self.metadata
        return pl.DataFrame([
            {
                "column": c,
                "kind": kind(md[c]),
                "dtype": str(md[c].dtype),
                "distinct": md[c].drop_nulls().n_unique(),
                "nulls": md[c].null_count(),
                "examples": ", ".join(str(x) for x in md[c].drop_nulls().unique(maintain_order=True).head(4).to_list()),
                "varies within patient": (bool(md.group_by("patient_id").agg(pl.col(c).n_unique())[c].max() > 1)
                                          if c != "patient_id" else None),
            }
            for c in md.columns
        ])

    # --- 3. clean up ------------------------------------------------------------------------

    def with_metadata(self, *exprs: pl.Expr, **named: pl.Expr) -> "Repertoires":
        """Add or replace metadata columns, as polars `with_columns`."""
        return replace(self, metadata=self.metadata.with_columns(*exprs, **named))

    # --- 4. QC ------------------------------------------------------------------------------

    def qc(self, *, min_clonotypes: int = 1000, max_zero: int = 10, outlier_iqr: float = 3) -> QCResult:
        return flag_repertoires(self, min_clonotypes=min_clonotypes, max_zero=max_zero, outlier_iqr=outlier_iqr)

    # --- 5-6. variables and samples ---------------------------------------------------------

    def design(self, *, label: str, batch: Optional[str], other_technical: Sequence[str] = (),
               covariates: Sequence[str] = (), controls: Sequence[str] = (), exclude: Sequence[str] = (),
               covariate_missing: str = "level", min_level_size: int = 5):
        """Fix the variables and the samples. See `StudyDesign.build`."""
        from .design import StudyDesign

        return StudyDesign.build(
            self, label=label, batch=batch, other_technical=list(other_technical), covariates=list(covariates),
            controls=list(controls), exclude=list(exclude), covariate_missing=covariate_missing,
            min_level_size=min_level_size,
        )
