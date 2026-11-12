"""Figures for the audit steps. Each takes a step's result and returns the matplotlib figure."""

from __future__ import annotations

import math
from typing import Optional, Sequence

import colorcet as cc
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np
import polars as pl
import seaborn as sns
from sklearn.metrics import (
    average_precision_score,
    precision_recall_curve,
    roc_auc_score,
    roc_curve,
)

from .baseline import BaselineResult, ShapResult
from .design import BalanceResult
from .features import GERMLINE_V_GENES
from .qc import FAILURES, QCResult
from .variance import PCAResult, VarianceResult

BLUE, ORANGE, AQUA, PINK, GREY = "#2a78d6", "#eb6834", "#1baf7a", "#e87ba4", "#a3a29d"
CATEGORICAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]  # validated order
SEQUENTIAL = sns.light_palette(BLUE, as_cmap=True)
DIVERGING = sns.blend_palette([ORANGE, "#f1f0ed", BLUE], as_cmap=True)  # AUC around 0.5
FEATURE_VALUE = sns.blend_palette([BLUE, "#d6d5d0", PINK], as_cmap=True)  # a feature's low to high values
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
    # "PC1 (37%)"; with several blocks also the blocks driving it, e.g. "PC1 (37%) – 91% V genes".
    label = f"{pc} ({result.explained[int(pc.removeprefix('PC')) - 1]:.0%})"
    shares = result.block_shares.select("block", pc).sort(pc, descending=True)
    if shares.height > 1:
        main = [(b, s) for b, s in shares.iter_rows() if s >= 0.2] or shares.head(1).rows()
        label += " – " + ", ".join(f"{s:.0%} {b}" for b, s in main)
    return label


def _of_block(result: PCAResult) -> str:
    return f" of {result.block}" if result.block else ""


# Muted versions of the palette's last six hues (red, violet, green, magenta, yellow, aqua), so blocks don't read as
# the label colours of the PCA panel beside them: OKLCH chroma x 0.55, lightness alternating dark / light, so touching
# segments (blocks are coloured by stack position) stay apart. Validated as adjacent pairs: CVD dE >= 16, normal >= 21.
BLOCK_COLOURS = ["#c4726c", "#5c599e", "#74a86f", "#8d5669", "#d2a970", "#478367"]


def _block_bars(result: PCAResult, x: str, y: str, ax):
    # Each PC's share of the total variance, split by block in proportion to the block's share of the PC's loadings:
    # a segment is the variance the PC captures from that block.
    shares = result.block_shares
    pcs = [c for c in shares.columns if c.startswith("PC")]
    bottom, handles = np.zeros(len(pcs)), []
    for i, (block, *share) in enumerate(shares.iter_rows()):
        colour = BLOCK_COLOURS[i % len(BLOCK_COLOURS)]
        part = np.array(share) * result.explained[: len(pcs)]
        captured = part.sum() / result.block_totals[block]  # of the block's own variance
        handles.append(ax.bar(range(len(pcs)), part, bottom=bottom, color=colour, width=0.75,
                              edgecolor="black", linewidth=0.5, label=f"{block}: {captured:.0%}"))
        bottom += part
    ax.set_xticks(range(len(pcs)), pcs, fontsize=8)
    for tick in ax.get_xticklabels():
        tick.set_fontweight("bold" if tick.get_text() in (x, y) else "normal")
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{100 * v:g}%"))  # 2.5%, 5%, ...: no rounded ticks
    ax.spines[["top", "right"]].set_visible(False)
    ax.set(ylabel="variance explained", title="Variance per PC, by block")
    ax.legend(handles=handles[::-1], title=f"block: its variance in PC1-{len(pcs)}", loc="upper right", fontsize=8,
              title_fontsize=8, frameon=False)


def pca(result: PCAResult, x: str = "PC1", y: str = "PC2", highlight: Sequence[str] = GERMLINE_V_GENES):
    """Manuscript Fig. 2B: two of the leading PCs, and the features driving them. Features in
    `highlight` (by default V genes whose usage varies for germline reasons) are marked in orange.
    With several blocks, a third panel splits each PC's explained variance by block, and gives how much of
    each block's own variance the PCs capture."""
    palette, colour = _colour(result.scores, result.label)
    pcs = result.scores.with_columns(colour.alias("colour"))
    top = result.table.with_columns((pl.col(x) ** 2 + pl.col(y) ** 2).alias("length")).sort("length", descending=True).head(12)
    marked = result.table.filter(pl.col("feature").is_in(list(highlight)))
    blocks = result.block_shares.height > 1
    fig, axes = plt.subplots(1, 3 if blocks else 2, figsize=(20 if blocks else 14, 5.5),
                             width_ratios=[1.25, 1, 0.9] if blocks else [1.25, 1])
    ax, ax_load = axes[:2]
    if blocks:
        _block_bars(result, x, y, axes[2])
    sns.scatterplot(pcs, x=x, y=y, hue="colour", hue_order=list(palette), palette=palette,
                    s=14, linewidth=0, alpha=0.8, ax=ax)
    ax.set(xlabel=_explained(result, x), ylabel=_explained(result, y), title=f"PCA{_of_block(result)}, coloured by {result.label}")
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
    fig.suptitle(f"{x} vs {y}{_of_block(result)} per {result.batch}")
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


def roc_pr(result: BaselineResult, variable: Optional[str] = None, classes: Optional[Sequence[str]] = None):
    """Out-of-fold ROC and precision-recall curves of the baseline classifier: one per class of `variable` (by
    default the label), each against the rest; `classes` picks some. Colours from colorcet's glasbey_dark, which
    stays distinct for many classes. Chance (thin, dotted) is the diagonal (ROC) and, in the class's colour, the
    class's share of the repertoires (PR). With `cv_repeats` > 1 the scores are averaged over the repeats, as for
    the table, so the curves' AUC and AP match it."""
    if variable is None:
        labels = result.table.filter(pl.col("role") == "label")["variable"]
        if labels.is_empty():
            raise ValueError("no label class was predicted; see base.table")
        variable = labels[0]
    fits = {c: f for (v, c), f in result.fits.items() if v == variable and (classes is None or c in classes)}
    if not fits:
        raise ValueError(f"no predicted class of {variable!r} matches; see base.table")
    fig, (ax_roc, ax_pr) = plt.subplots(1, 2, figsize=(11, 6 + 0.12 * len(fits)), layout="constrained")
    handles = []
    for i, (c, fit) in enumerate(fits.items()):  # colour by class, in a fixed order, not by rank
        colour = cc.glasbey_dark[i % len(cc.glasbey_dark)]
        ok = ~np.isnan(fit.scores)
        y, s = fit.y[ok], fit.scores[ok]
        fpr, tpr, _ = roc_curve(y, s)
        handles += ax_roc.plot(fpr, tpr, color=colour, linewidth=2,
                               label=f"{c}: AUC {roc_auc_score(y, s):.2f}, AP {average_precision_score(y, s):.2f}"
                                     f" (chance {y.mean():.2f})")
        precision, recall, _ = precision_recall_curve(y, s)
        ax_pr.step(recall, precision, where="post", color=colour, linewidth=2)
        ax_pr.axhline(y.mean(), color=colour, linewidth=0.8, linestyle=":")
    handles += ax_roc.plot([0, 1], [0, 1], color="0.5", linewidth=0.8, linestyle=":",
                           label="chance: the diagonal (ROC), the class's share (PR, in its colour)")
    ax_roc.set(xlabel="false positive rate", ylabel="true positive rate", title="ROC")
    ax_pr.set(xlabel="recall", ylabel="precision", title="Precision-recall")
    for ax in (ax_roc, ax_pr):
        ax.set(xlim=(0, 1), ylim=(0, 1.02), aspect="equal")
    folds = f"folds: {result.cv_column}" if result.cv_column else "random folds"
    fig.suptitle(f"Baseline classifier, {variable}: each class against the rest, out of fold ({folds})")
    fig.legend(handles=handles, loc="outside lower center", ncol=2 if len(handles) > 5 else 1, fontsize=8,
               frameon=False)
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


def _swarm(x: np.ndarray, width: float = 0.4, nbins: int = 100, seed: int = 0) -> np.ndarray:
    # Vertical offsets that pile up dots with similar x, as in shap's beeswarm, at most `width` from the row.
    bins = np.round(nbins * (x - x.min()) / (np.ptp(x) or 1.0))
    order = np.lexsort((np.random.default_rng(seed).random(len(x)), bins))  # by bin, random within one
    offsets, layer, last = np.zeros(len(x)), 0, None
    for i in order:
        layer = layer + 1 if bins[i] == last else 0
        offsets[i], last = math.ceil(layer / 2) * (1 if layer % 2 else -1), bins[i]
    return offsets * width / max(np.abs(offsets).max(), 1)


def _scaled(v: np.ndarray) -> np.ndarray:
    # A feature's values on 0-1 between its 5th and 95th percentile, so a few outliers don't wash out the colours.
    lo, hi = np.percentile(v, [5, 95])
    if hi <= lo:
        lo, hi = v.min(), v.max()
    return np.clip((v - lo) / (hi - lo), 0, 1) if hi > lo else np.full(len(v), 0.5)


def _dot(colour, label: str, marker: str = "o"):
    return plt.Line2D([], [], marker=marker, linestyle="", color=colour, label=label)


def _batch_palette(result: ShapResult) -> dict:
    # Batches beyond the validated palette fold into "other". The batches the model treats most differently keep
    # a colour: those whose mean summed SHAP value (= their mean prediction) is furthest from the average, so a
    # small batch that carries the signal still stands out.
    total = result.values.sum(axis=1) + result.others
    ranked = (pl.DataFrame({"batch": result.batch, "total": total})
              .group_by("batch").agg((pl.col("total").mean() - float(total.mean())).abs().alias("distance"))
              .sort("distance", "batch", descending=[True, False])["batch"].to_list())
    return dict(zip(ranked if len(ranked) <= len(CATEGORICAL) else ranked[: len(CATEGORICAL) - 1], CATEGORICAL))


def shap_beeswarm(result: ShapResult, by: str = "value", others: bool = True):
    """The repertoires behind `shap`: one dot per repertoire and feature, at its out-of-fold SHAP value (how
    much that feature moved the repertoire's predicted probability of the class). Dots pile up where values
    are similar; their height within a row means nothing else. Same features, same order as `shap`.

    by="value": coloured by the feature's value in that repertoire (pink high, blue low).
    by="batch": coloured by the repertoire's batch, triangles for repertoires of the class: which batches
    carry the signal. With more batches than colours, those with the largest SHAP values keep a colour.
    others: a last row with, per repertoire, the summed SHAP values of all the features not shown (grey
    when coloured by value). Show more features with `design.shap(base, top=...)`."""
    if by not in ("value", "batch"):
        raise ValueError('by must be "value" or "batch"')
    if by == "batch" and result.batch is None:
        raise ValueError('the design has no batch; use by="value"')
    rows = [(f, result.values[:, j], result.data[:, j]) for j, f in enumerate(result.table["feature"])]
    if others and result.n_others:
        rows.append((f"sum of {result.n_others} other features", result.others, None))
    order = np.random.default_rng(0).permutation(len(result.in_target))  # random draw order: no batch always on top
    if by == "batch":
        palette = _batch_palette(result)
        batch_colour = np.array([palette.get(b, GREY) for b in result.batch], dtype=object)
    fig, ax = plt.subplots(figsize=(11 if by == "batch" else 8, 0.35 * len(rows) + 1.5), layout="constrained")
    for i, (_, x, v) in enumerate(rows):
        row = len(rows) - 1 - i  # the first feature at the top
        y = row + _swarm(x)
        ax.axhline(row, color="#ebeae6", linewidth=0.6, zorder=0)
        if by == "value":
            ax.scatter(x[order], y[order], color=FEATURE_VALUE(_scaled(v))[order] if v is not None else GREY,
                       s=10, linewidths=0)
        else:
            for marker, size, member in (("o", 10, ~result.in_target), ("^", 20, result.in_target)):
                idx = order[member[order]]
                ax.scatter(x[idx], y[idx], color=list(batch_colour[idx]), marker=marker, s=size, linewidths=0, alpha=0.85)
    ax.axvline(0, color="0.6", linewidth=0.8)
    ax.set_yticks(range(len(rows)), [name for name, _, _ in rows][::-1])
    ax.set_ylim(-0.7, len(rows) - 0.3)
    ax.set(xlabel=f"SHAP value: change in predicted probability of {result.target} (out of fold)", ylabel="",
           title=f"{result.label} = {result.target} vs rest (AUC {result.auc:.2f})")
    if by == "value":
        cbar = fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(0, 1), cmap=FEATURE_VALUE), ax=ax,
                            ticks=[0, 1], aspect=40, pad=0.02)
        cbar.ax.set_yticklabels(["low", "high"])
        cbar.set_label("feature value")
        cbar.outline.set_visible(False)
    else:
        n_other = len(set(result.batch)) - len(palette)
        batches = [_dot(c, b) for b, c in palette.items()] + ([_dot(GREY, f"other ({n_other} batches)")] if n_other else [])
        fig.legend(handles=batches, title="batch", loc="outside right upper", fontsize=8)
        fig.legend(handles=[_dot("0.35", result.target, "^"), _dot("0.35", "rest")], title=result.label,
                   loc="outside right lower", fontsize=8)
    return fig
