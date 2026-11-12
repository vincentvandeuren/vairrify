"""Feature blocks: per-repertoire tables the audit can use, and where each comes from in tcrio.

A block is one tcrio result read as one row per repertoire. Its kind follows the manuscript:
"technical" blocks are the general repertoire properties that technical bias shapes (repertoire
size, gene usage, VDJ recombination statistics); "biological" blocks are confounders such as
inferred HLA. An analysis uses one block, or several side by side.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Optional
import warnings

import numpy as np
import polars as pl
from tcr_io import operations as ops
from tcr_io.dataset import Dataset
from tcr_io.structure import Artifact


def gene_usage_profile(freq: np.ndarray, clonotypes: np.ndarray) -> np.ndarray:
    # The manuscript's gene-usage profile: log1p of the counts, L2-normalised per repertoire.
    log_counts = np.log1p(freq * clonotypes[:, None])
    with np.errstate(invalid="ignore", divide="ignore"):
        return log_counts / np.linalg.norm(log_counts, axis=1, keepdims=True)


@dataclass(frozen=True)
class FeatureBlock:
    name: str
    kind: str                                         # "technical" or "biological"
    description: str
    read_table: Callable[[Dataset], pl.DataFrame]     # repertoire_id + one column per feature
    artifact: Optional[Artifact] = None               # the tcrio result read; None = read_table alone
    profile: Optional[Callable[[np.ndarray, np.ndarray], np.ndarray]] = None  # what PCA and variance see

    @property
    def operation(self) -> Optional[type]:
        return self.artifact.owner if self.artifact is not None else None

    def supports(self, locus: str) -> bool:
        return self.operation is None or locus in self.operation.supported_loci

    def read(self, ds: Dataset, *, run_operations: bool = True) -> pl.DataFrame:
        if not self.supports(ds.locus):
            raise ValueError(f"{self.name!r} comes from {self.operation.__name__}, which does not run on {ds.locus}")
        if self.operation is not None and run_operations:
            # tcrio reruns only a missing or stale result (e.g. after a version bump), and warns otherwise
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", message=".*already complete.*")
                ds.run_operation(self.operation())
        elif self.artifact is not None and not ds.get_result(self.artifact).exists():
            raise LookupError(f"{self.name!r} needs {self.operation.__name__}; run it, or pass run_operations=True")
        table = self.read_table(ds)
        return table.select("repertoire_id", pl.exclude("repertoire_id").cast(pl.Float64))


def _result(artifact: Artifact, transform: Callable[[pl.DataFrame], pl.DataFrame] = lambda t: t):
    return lambda ds: transform(ds.get_result(artifact).read())


# V genes whose usage varies between people for germline reasons, so their variation needn't be technical:
# deletion polymorphisms (Omer et al. 2022, Genome Med 14:2) and their linked alleles, confirmed as bimodal
# usage in three large cohorts (see the docs, "Germline variation"). tcrio stores the identical
# TRBV6-2 / TRBV6-3 pair as TRBV6-2.
GERMLINE_V_GENES = ("TRBV4-3", "TRBV6-2", "TRBV6-3", "TRBV7-2", "TRBV7-3", "TRBV28", "TRBV30")

RAREFACTION_DEPTHS =(20_000, 48_000, 200_000, 500_000)  # on tcrio's rarefaction grid


def _depth(ds: Dataset) -> pl.DataFrame:
    # Reads, and the expected unique clonotypes at fixed read depths: interpolated in log depth
    # between the rarefaction grid's points, extrapolated past a repertoire's own reads.
    curves = ds.get_result(ops.RarefactionReport.rarefaction_curves).read().sort("subsampling_depth")
    rows = []
    for (rep,), c in curves.group_by("repertoire_id"):
        clonotypes = np.interp(np.log(RAREFACTION_DEPTHS), np.log(c["subsampling_depth"]), c["expected_richness"])
        rows.append({"repertoire_id": rep, "reads": c["total_abundance"][0],
                     **{f"clonotypes at {d // 1000}k reads": n for d, n in zip(RAREFACTION_DEPTHS, clonotypes)}})
    return pl.DataFrame(rows)


FEATURE_BLOCKS: dict[str, FeatureBlock] = {b.name: b for b in (
    FeatureBlock("V genes", "technical", "share of clonotypes using each V gene",
                 _result(ops.GeneCountsSummary.v_counts), ops.GeneCountsSummary.v_counts, gene_usage_profile),
    FeatureBlock("J genes", "technical", "share of clonotypes using each J gene",
                 _result(ops.GeneCountsSummary.j_counts), ops.GeneCountsSummary.j_counts, gene_usage_profile),
    FeatureBlock("VJ bias", "technical", "log ratio of each V-J pairing to what independent V and J usage predicts",
                 _result(ops.GeneCountsSummary.vj_bias), ops.GeneCountsSummary.vj_bias),
    FeatureBlock("VDJ stats", "technical", "mean, sd, skew and kurtosis of trimming, P and N nucleotides",
                 _result(ops.VdjStatisticsSummary.summary, lambda t: t.select(pl.exclude("^.*_sum$"))),
                 ops.VdjStatisticsSummary.summary),
    FeatureBlock("depth", "technical", "reads, and expected unique clonotypes at 20k / 48k / 200k / 500k reads (rarefaction)",
                 _depth, ops.RarefactionReport.rarefaction_curves),
    # dichotomised as in the manuscript, so repertoire size does not leak in through the scores
    FeatureBlock("HLA", "biological", "inferred HLA carriage (score > 0.5)",
                 _result(ops.RepertoireHlaInference.inferred_hla,
                         lambda t: t.select("repertoire_id", (pl.exclude("repertoire_id") > 0.5).cast(pl.Float64))),
                 ops.RepertoireHlaInference.inferred_hla),
)}


def feature_blocks() -> pl.DataFrame:
    """The blocks a feature set can be built from."""
    return pl.DataFrame([
        {"block": b.name, "kind": b.kind, "source": b.operation.__name__ if b.operation else "read_table",
         "description": b.description}
        for b in FEATURE_BLOCKS.values()
    ])
