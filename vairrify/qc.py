"""Step 4: repertoires that shouldn't be analysed.

A few failed libraries can dominate everything downstream. Each repertoire is checked for:

- low depth: fewer than `min_clonotypes` clonotypes
- missing genes: more than `max_zero` genes never observed, in the gene-usage blocks
- atypical profile: more than `outlier_iqr` IQRs above the upper quartile of the distance to the
  typical profile (root-mean-square robust z-score over the features, each block counting equally)

The first two are failures. An atypical profile is only a reason to look: it can be biology.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
import polars as pl

from .features import FEATURE_BLOCKS, gene_usage_profile
from .results import Result

if TYPE_CHECKING:
    from .repertoires import Repertoires

FAILURES = ("low depth", "missing genes")


@dataclass(frozen=True, repr=False)
class QCResult(Result):
    cutoff: float
    min_clonotypes: int

    def _ids(self, expr: pl.Expr) -> list[str]:
        return self.table.filter(expr)["repertoire_id"].to_list()

    @property
    def failed(self) -> list[str]:
        return self._ids(pl.col("flags").str.contains("|".join(FAILURES)))

    @property
    def atypical(self) -> list[str]:
        return self._ids(pl.col("flags").str.contains("atypical"))

    @property
    def flagged(self) -> list[str]:
        return self._ids(pl.col("flags") != "")


def flag_repertoires(reps: Repertoires, *, min_clonotypes: int, max_zero: int, outlier_iqr: float) -> QCResult:
    has = reps.present
    profile = np.where(has[:, None], reps.profile(), np.nan)
    typical = np.nanmedian(profile, axis=0)
    spread = 1.4826 * np.nanmedian(np.abs(profile - typical), axis=0)
    spread = np.where(spread > 0, spread, np.nanstd(profile, axis=0) + 1e-12)
    z2 = ((profile - typical) / spread) ** 2
    distance = np.sqrt((z2 * reps.weights(profile) ** 2).sum(axis=1) / len(reps.features))
    q1, q3 = np.nanpercentile(distance, [25, 75])
    cutoff = q3 + outlier_iqr * (q3 - q1)
    genes = np.isin(reps.block_of, [b for b in reps.features if FEATURE_BLOCKS[b].profile is gene_usage_profile])
    zeros = (reps.X[:, genes] == 0).sum(axis=1) if genes.any() else np.full(len(has), None)

    table = pl.DataFrame({  # NaN -> null: polars orders NaN above every number, so NaN > cutoff would be true
        "repertoire_id": reps.ids,
        "clonotypes": pl.Series(reps.clonotypes).fill_nan(None),
        "genes at zero": pl.Series(zeros, dtype=pl.Int64),
        "distance": pl.Series(distance).fill_nan(None),
    }).with_columns(
        pl.concat_str([
            pl.when(pl.col("clonotypes") < min_clonotypes).then(pl.lit("low depth")),
            pl.when(pl.col("genes at zero") > max_zero).then(pl.lit("missing genes")),
            pl.when(pl.col("distance") > cutoff).then(pl.lit("atypical profile")),
        ], separator=", ", ignore_nulls=True).alias("flags")
    ).filter(pl.Series(has))  # repertoires without data for this locus drop out at the design step anyway
    return QCResult(table.sort("distance", descending=True, nulls_last=True), float(cutoff), min_clonotypes)
