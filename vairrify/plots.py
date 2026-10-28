"""Figures for the audit steps. Each takes a step's result and returns the matplotlib figure."""

from __future__ import annotations

import math
from typing import Sequence

import matplotlib.pyplot as plt
import polars as pl
import seaborn as sns

from .baseline import BaselineResult, ShapResult
from .design import BalanceResult
from .features import GERMLINE_V_GENES
from .qc import FAILURES, QCResult
from .variance import PCAResult, VarianceResult

BLUE, ORANGE, AQUA, PINK, GREY = "#2a78d6", "#eb6834", "#1baf7a", "#e87ba4", "#a3a29d"
CATEGORICAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]  # validated order
SEQUENTIAL = sns.light_palette(BLUE, as_cmap=True)
DIVERGING = sns.blend_palette([ORANGE, "#f1f0ed", BLUE], as_cmap=True)  # AUC around 0.5
ROLE_COLORS = {"technical": BLUE, "label": ORANGE, "covariate": AQUA, "control": GREY, "residual / interactions": "#d6d5d0"}


def qc(result: QCResult):
    status = pl.when(pl.col("flags").str.contains("|".join(FAILURES))).then(pl.lit("failed")) \
        .when(pl.col("flags") != "").then(pl.lit("atypical")).otherwise(pl.lit("ok"))
    df = result.table.with_columns(pl.col("clonotypes").log10().alias("log10 clonotypes"), status.alias("status"))
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    sns.scatterplot(df, x="log10 clonotypes", y="distance", hue="status", hue_order=["ok", "atypical", "failed"],
                    palette={"ok": "#c9c8c3", "atypical": BLUE, "failed": ORANGE}, s=18, linewidth=0, ax=ax)
    ax.axhline(result.cutoff, color="0.4", linestyle="--", linewidth=1)
    ax.axvline(math.log10(result.min_clonotypes), color="0.4", linestyle=":", linewidth=1)
    ax.set(ylabel="distance to typical profile", title="QC: depth vs atypical profile")
    ax.legend(title="", loc="upper right")
    fig.tight_layout()
    return fig


def balance(result: BalanceResult):
    """Repertoires per batch and label level (manuscript Fig. S1); single-label batches are marked."""
    counts = result.table.with_columns(
        pl.when(pl.col("batch type") == "single-label").then(pl.format("{} (single-label)", result.batch))
        .otherwise(pl.col(result.batch)).alias(result.batch)).drop("batch type").to_pandas().set_index(result.batch)
    mixed = int((result.table["batch type"] == "mixed-label").sum())
    fig, ax = plt.subplots(figsize=(0.75 * counts.shape[1] + 6, 0.4 * counts.shape[0] + 2))
    sns.heatmap(counts, annot=True, fmt="d", cmap=SEQUENTIAL, linewidths=2, linecolor="white",
                cbar_kws={"label": "repertoires"}, ax=ax)
    ax.set(xlabel=result.label, ylabel="", title=f"{result.label} per {result.batch}: {result.design} "
           f"(Cramér's V {result.cramers_v:.2f}), {mixed} of {result.table.height} batches mixed-label")
    fig.tight_layout()
    return fig


def _colour(scores: pl.DataFrame, label: str):
    # Label levels beyond the validated palette fold into "other".
    levels = scores[label].value_counts(sort=True)[label].to_list()
    shown = levels if len(levels) <= len(CATEGORICAL) else levels[: len(CATEGORICAL) - 1]
    palette = dict(zip(shown, CATEGORICAL)) | ({"other": GREY} if len(shown) < len(levels) else {})
    return palette, pl.when(pl.col(label).is_in(shown)).then(pl.col(label)).otherwise(pl.lit("other"))


def _explained(result: PCAResult, pc: str) -> str:
    return f"{pc} ({result.explained[int(pc.removeprefix('PC')) - 1]:.0%})"


def pca(result: PCAResult, x: str = "PC1", y: str = "PC2", highlight: Sequence[str] = GERMLINE_V_GENES):
    """Manuscript Fig. 2B: two of the leading PCs, and the features driving them. Features in
    `highlight` (by default V genes whose usage varies for germline reasons) are marked in orange."""
    palette, colour = _colour(result.scores, result.label)
    pcs = result.scores.with_columns(colour.alias("colour"))
    top = result.table.with_columns((pl.col(x) ** 2 + pl.col(y) ** 2).alias("length")).sort("length", descending=True).head(12)
    marked = result.table.filter(pl.col("feature").is_in(list(highlight)))
    fig, (ax, ax_load) = plt.subplots(1, 2, figsize=(14, 5.5), width_ratios=[1.25, 1])
    sns.scatterplot(pcs, x=x, y=y, hue="colour", hue_order=list(palette), palette=palette,
                    s=14, linewidth=0, alpha=0.8, ax=ax)
    ax.set(xlabel=_explained(result, x), ylabel=_explained(result, y), title=f"PCA, coloured by {result.label}")
    ax.legend(title="", loc="upper left", bbox_to_anchor=(1.01, 1))
    sns.scatterplot(result.table, x=x, y=y, color="#c9c8c3", s=14, linewidth=0, ax=ax_load)
    sns.scatterplot(top, x=x, y=y, color=BLUE, s=24, linewidth=0, label="top 12", ax=ax_load)
    if marked.height:
        sns.scatterplot(marked, x=x, y=y, color=ORANGE, s=30, linewidth=0, label="germline variation", ax=ax_load)
    for feature, fx, fy in pl.concat([top.select("feature", x, y), marked.select("feature", x, y)]).unique("feature").iter_rows():
        ax_load.annotate(feature, (fx, fy), xytext=(3, 3), textcoords="offset points", fontsize=8)
    ax_load.axhline(0, color="0.6", linewidth=0.8)
    ax_load.axvline(0, color="0.6", linewidth=0.8)
    ax_load.set(xlabel=f"{x} loading", ylabel=f"{y} loading", title=f"Features driving {x} and {y}")
    ax_load.legend(title="", loc="upper left", bbox_to_anchor=(1.01, 1), fontsize=8)
    fig.tight_layout()
    return fig


def pca_by_batch(result: PCAResult, x: str = "PC1", y: str = "PC2", ncols: int = 4):
    """One panel per batch: its repertoires coloured over all the others in grey."""
    if not result.batch:
        return None
    palette, colour = _colour(result.scores, result.label)
    pcs = result.scores.with_columns(colour.alias("colour"))
    batches = pcs[result.batch].value_counts(sort=True)[result.batch].to_list()
    nrows = math.ceil(len(batches) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(15, 3.3 * nrows), sharex=True, sharey=True, squeeze=False)
    for ax, b in zip(axes.flat, batches):
        sns.scatterplot(pcs, x=x, y=y, color="#e3e2dd", s=8, linewidth=0, ax=ax)
        sns.scatterplot(pcs.filter(pl.col(result.batch) == b), x=x, y=y, hue="colour", hue_order=list(palette),
                        palette=palette, s=16, linewidth=0, legend=False, ax=ax)
        ax.set_title(f"{b} (n = {pcs.filter(pl.col(result.batch) == b).height})", fontsize=9)
    for ax in axes.flat[len(batches):]:
        ax.set_visible(False)
    handles = [plt.Line2D([], [], marker="o", linestyle="", color=c, label=l) for l, c in palette.items()]
    fig.legend(handles=handles, title=result.label, loc="upper left", bbox_to_anchor=(1.0, 0.98))
    fig.suptitle(f"{x} vs {y} per {result.batch}")
    fig.tight_layout()
    return fig


def variance(result: VarianceResult):
    """Manuscript Fig. 2A (PVCA), next to PERMANOVA R² alone and adjusted."""
    fig, (ax_pvca, ax_perm) = plt.subplots(1, 2, figsize=(14, 0.45 * result.table.height + 2.5))
    if result.pvca is not None:
        continuous = result.table.filter(pl.col("kind") == "continuous")["variable"].to_list()
        pvca = result.pvca.with_columns(
            pl.when(pl.col("effect").is_in(continuous)).then(pl.format("{} (binned)", "effect")).otherwise(pl.col("effect")).alias("effect"))
        sns.barplot(pvca, y="effect", x="proportion", hue="role", palette=ROLE_COLORS, dodge=False, ax=ax_pvca)
        ax_pvca.set(xlabel="share of variance", ylabel="", title=f"PVCA ({result.n_pcs} PCs)")
        ax_pvca.legend(title="", loc="lower right")
    else:
        ax_pvca.set_axis_off()
        ax_pvca.set_title("no PVCA for a single feature")
    perm = (result.table.unpivot(index=["variable", "role"], on=["R² alone", "R² adjusted"], variable_name="test", value_name="r2")
            .with_columns(pl.format("{} ({})", "variable", "role").alias("variable")))
    order = perm.filter(pl.col("test") == "R² alone").sort("r2", descending=True)["variable"].to_list()
    sns.barplot(perm, y="variable", x="r2", hue="test", order=order, hue_order=["R² alone", "R² adjusted"],
                palette={"R² alone": "#c9c8c3", "R² adjusted": "#52514e"}, ax=ax_perm)  # neutrals: blue/orange mean roles on the left
    ax_perm.set(xlabel="R² (share of total variation)", ylabel="", title=f"PERMANOVA (all variables: R² = {result.model_r2:.2f})")
    ax_perm.legend(title="", loc="lower right")
    fig.tight_layout()
    return fig


def baseline(result: BaselineResult):
    """Manuscript Table 1 as a figure: out-of-fold AUC per class."""
    df = result.table.with_columns(pl.format("{}: {}", "variable", "class").alias("target"))
    fig, ax = plt.subplots(figsize=(8, 0.3 * df.height + 1.6))
    sns.barplot(df, y="target", x="auc", hue="role", palette=ROLE_COLORS, dodge=False, ax=ax)
    ax.axvline(0.5, color="0.3", linewidth=1, linestyle="--")
    ax.set(xlim=(0, 1), xlabel="out-of-fold ROC AUC (0.5 = chance)", ylabel="",
           title=f"Baseline classifier ({'folds: ' + result.cv_column if result.cv_column else 'random folds'})")
    ax.legend(title="", loc="upper left", bbox_to_anchor=(1.01, 1))
    fig.tight_layout()
    return fig


def within_batch(table: pl.DataFrame):
    """AUC with all batches combined vs within mixed-label batches (blank: no batch holds both classes)."""
    wide = table.select("class", "all batches", "within batch").to_pandas().set_index("class")
    fig, ax = plt.subplots(figsize=(6.5, 0.45 * len(wide) + 1.8))
    sns.heatmap(wide, annot=True, fmt=".2f", cmap=DIVERGING, center=0.5, vmin=0, vmax=1, linewidths=2, linecolor="white",
                cbar_kws={"label": "AUC"}, ax=ax)
    ax.grid(False)
    ax.set_facecolor("#f7f7f5")
    ax.set(xlabel="", ylabel="", title="Within- vs across-batch evaluation")
    fig.tight_layout()
    return fig


def shap(result: ShapResult):
    """Manuscript Fig. 3B: the features carrying a label class, coloured by whether they are positively
    associated with it (a higher value pushes towards the class, pink) or negatively (blue)."""
    more, less = f"positively associated with {result.target}", f"negatively associated with {result.target}"
    df = result.table.with_columns(pl.col("association").replace({"positive": more, "negative": less}))
    fig, ax = plt.subplots(figsize=(7.5, 0.3 * df.height + 1.5))
    sns.barplot(df, y="feature", x="mean_abs_shap", hue="association", palette={more: PINK, less: BLUE, "none": GREY},
                dodge=False, ax=ax)
    ax.set(xlabel="mean |SHAP| (out of fold)", ylabel="",
           title=f"{result.label} = {result.target} vs rest (AUC {result.auc:.2f})")
    ax.legend(title="", loc="lower right")
    fig.tight_layout()
    return fig
