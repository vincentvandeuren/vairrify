"""Steps 8-9: what drives variation in the features?

- PCA (manuscript Fig. 2B): the leading components, and the features that load most on them.
- PVCA (Fig. 2A): variance components averaged over the leading PCs. Continuous variables are
  binned for it, with the given edges or into equal-sized groups.
- PERMANOVA: each variable's share of the total variation, alone and adjusted for all others.
  A label that is large alone but small adjusted shows up mainly through what it is tied to,
  usually the batch. Controls should sit near zero.

PCA and PVCA need at least two features. PVCA shares and PERMANOVA R² answer different questions
and don't add up.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Mapping, Optional, Sequence

import numpy as np
import polars as pl
from sklearn.decomposition import PCA

from . import stat_utils as su
from .results import Result

if TYPE_CHECKING:
    from .design import StudyDesign


@dataclass(frozen=True, repr=False)
class PCAResult(Result):
    """`table`: each feature's block and loadings on PC1, PC2, ...; `scores` the repertoires' metadata and PCs."""
    scores: pl.DataFrame
    explained: np.ndarray
    label: str
    batch: Optional[str]
    block: Optional[str] = None   # the one block analysed; None = all loaded blocks together
    block_totals: Optional[Mapping[str, float]] = None  # each block's share of the total variance (1 / blocks by design)

    @property
    def block_shares(self) -> pl.DataFrame:
        """Per block, its share of each PC's loadings (squared loadings; each PC's shares sum to 1). With
        several blocks: which blocks drive a PC."""
        pcs = [c for c in self.table.columns if c.startswith("PC")]
        return self.table.group_by("block", maintain_order=True).agg((pl.col(pcs) ** 2).sum())


@dataclass(frozen=True, repr=False)
class VarianceResult(Result):
    """`table`: per variable, its PVCA share and PERMANOVA R² alone / adjusted."""
    pvca: Optional[pl.DataFrame]   # effect, proportion (interactions summed), role; None with a single feature
    n_pcs: Optional[int]
    model_r2: float


def pca(design: StudyDesign, n_components: int, block: Optional[str] = None) -> PCAResult:
    # The design's standardised profile; one block's columns of it for a per-block PCA.
    reps = design.repertoires
    if block is not None and block not in reps.features:
        raise ValueError(f"{block!r} is not a loaded block; loaded: {reps.features}")
    cols = np.ones(len(reps.feature_names), dtype=bool) if block is None else reps.block_of == block
    X = design.structure.X[:, cols]
    if X.shape[1] < 2:
        raise ValueError("PCA needs at least two features")
    fit = PCA(n_components=min(n_components, X.shape[1])).fit(X)
    s = fit.transform(X)
    pcs = [f"PC{i + 1}" for i in range(fit.n_components_)]
    scores = design.structure.metadata.with_columns(*[pl.Series(pc, s[:, i]) for i, pc in enumerate(pcs)])
    loadings = pl.DataFrame({"feature": np.asarray(reps.feature_names)[cols], "block": reps.block_of[cols],
                             **{pc: fit.components_[i] for i, pc in enumerate(pcs)}})
    var = X.var(axis=0, ddof=1)  # as PCA's explained_variance_ratio_
    totals = {str(b): float(var[reps.block_of[cols] == b].sum() / var.sum()) for b in dict.fromkeys(reps.block_of[cols])}
    return PCAResult(loadings, scores, fit.explained_variance_ratio_, design.label, design.batch, block, totals)


def prepare(X, md, continuous, categorical, missing) -> tuple[su.Design, list[str]]:
    # prepare_design, dropping any categorical variable left with a single level.
    # TODO: an option on stat_utils.prepare_design would replace matching on the error text.
    categorical, dropped = list(categorical), []
    while True:
        try:
            return su.prepare_design(X, md, continuous=continuous, categorical=categorical, missing=missing), dropped
        except ValueError as err:
            single = [c for c in categorical if f"{c!r} has fewer than two levels" in str(err)]
            if not single:
                raise
            categorical.remove(single[0])
            dropped.append(single[0])


def _binned(values: np.ndarray, edges: Sequence[float]) -> np.ndarray:
    edges = np.asarray(edges, dtype=float)
    idx = np.clip(np.searchsorted(edges, values, side="right") - 1, 0, len(edges) - 2)
    names = [f"{lo:.3g}+" if np.isinf(hi) else f"{lo:.3g}-{hi:.3g}" for lo, hi in zip(edges[:-1], edges[1:])]
    return np.array(names)[idx]


def variance(design: StudyDesign, *, permutations: int, pvca_bins: Mapping[str, Sequence[float]],
             pvca_n_bins: int) -> VarianceResult:
    structure, variables, roles = design.structure, design.variables, dict(design.roles)
    permanova = su.permanova(structure, terms=variables, permutations=permutations, seed=design.repertoires.seed)
    table = (permanova.filter(pl.col("test").is_in(["alone", "marginal"])).pivot(on="test", index="term", values=["r2", "p_value"])
             .select(pl.col("term").alias("variable"), pl.col("r2_alone").alias("R² alone"),
                     pl.col("r2_marginal").alias("R² adjusted"), pl.col("p_value_marginal").alias("p adjusted")))

    pvca_df, n_pcs = None, None
    if structure.X.shape[1] > 1:
        pvca_md = structure.metadata  # PVCA models every variable as a random effect, so continuous ones are binned
        for v in [v for v in variables if design.kinds[v] == "continuous"]:
            x = pvca_md[v].to_numpy()
            edges = pvca_bins.get(v, np.unique(np.quantile(x, np.linspace(0, 1, pvca_n_bins + 1))))
            pvca_md = pvca_md.with_columns(pl.Series(v, _binned(x, edges)))
        pvca = su.pvca(prepare(structure.X, pvca_md, [], variables, "drop")[0])
        n_pcs = pvca.n_pcs
        table = table.join(pvca.summary.select(pl.col("effect").alias("variable"), pl.col("proportion").alias("PVCA share")),
                           on="variable", how="left")
        pvca_df = (pvca.summary.with_columns(
                       pl.when(pl.col("effect").str.contains(":")).then(pl.lit("interactions")).otherwise(pl.col("effect")).alias("effect"))
                   .group_by("effect").agg(pl.col("proportion").sum())
                   .with_columns(pl.col("effect").replace_strict(roles, default="residual / interactions").alias("role"))
                   .sort("proportion", descending=True))
    table = (table.with_columns(pl.col("variable").replace_strict(roles).alias("role"),
                                pl.col("variable").replace_strict(dict(design.kinds)).alias("kind"))
             .sort("R² alone", descending=True)
             .select("variable", "role", "kind", *(["PVCA share"] if pvca_df is not None else []), "R² alone", "R² adjusted", "p adjusted"))
    return VarianceResult(table, pvca_df, n_pcs, permanova.filter(pl.col("term") == "model")["r2"][0])
