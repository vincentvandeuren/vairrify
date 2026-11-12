"""Steps 10-11: the baseline classifier, and whether its signal holds within batches.

Baseline classifier (manuscript Table 1 / Fig. 3). A random forest predicts each label class against
the rest from the features, scored out of fold. With `cv_column` the folds are the downstream
model's own. Otherwise 3 random folds keep each patient's repertoires together and are stratified
on batch x label (so every fold keeps each batch's label mix), repeated `cv_repeats` times and
averaged. Batch levels (a positive control) and the controls (should sit at 0.5) are predicted too.

Within batch. The manuscript's within-batch evaluation, applied to the baseline: per label class,
the AUC with all batches combined next to the AUC judged inside the mixed-label batches that hold
enough of both classes.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Mapping, Optional

import numpy as np
import polars as pl
from sklearn.base import clone
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import GroupKFold, StratifiedGroupKFold

from .results import Result

if TYPE_CHECKING:
    from .design import StudyDesign

ROLE_ORDER = {"label": 0, "technical": 1, "control": 2}


def default_forest() -> RandomForestClassifier:
    return RandomForestClassifier(random_state=42, n_estimators=400, max_features=0.05, max_depth=10,
                                  n_jobs=-1, class_weight="balanced_subsample")


@dataclass(frozen=True)
class Fit:
    y: np.ndarray        # is the class
    splits: list         # the (first) CV split
    scores: np.ndarray   # out-of-fold P(class)


@dataclass(frozen=True, repr=False)
class BaselineResult(Result):
    """`table`: out-of-fold AUC and average precision per class."""
    rows: np.ndarray                       # the repertoires predicted, as positions in repertoires.metadata
    fits: Mapping[tuple[str, str], Fit]    # (variable, class) -> fit
    rf: RandomForestClassifier
    cv_column: Optional[str]


def oof_scores(rf, X: np.ndarray, y: np.ndarray, splits) -> np.ndarray:
    # Out-of-fold P(class); NaN where a training fold lacks either class.
    scores = np.full(len(y), np.nan)
    for train, test in splits:
        if 0 < y[train].sum() < len(train):
            scores[test] = clone(rf).fit(X[train], y[train]).predict_proba(X[test])[:, 1]
    return scores


def auc(y: np.ndarray, scores: np.ndarray) -> float:
    ok = ~np.isnan(scores)
    return roc_auc_score(y[ok], scores[ok]) if 0 < y[ok].sum() < ok.sum() else np.nan


def average_precision(y: np.ndarray, scores: np.ndarray) -> float:
    ok = ~np.isnan(scores)
    return average_precision_score(y[ok], scores[ok]) if 0 < y[ok].sum() < ok.sum() else np.nan


def _splits(design: StudyDesign, rows: np.ndarray, y: np.ndarray, cv_column: Optional[str], repeat: int) -> list:
    # The downstream model's folds, or random folds that keep patients together and each batch's label mix.
    md = design.repertoires.metadata
    if cv_column:
        groups = md[cv_column].cast(pl.String).to_numpy()[rows]
        return list(GroupKFold(n_splits=len(np.unique(groups))).split(rows, y, groups))
    strata = y.astype(str)
    if design.batch:
        strata = np.char.add(md[design.batch].cast(pl.String).to_numpy()[rows].astype(str), strata)
    seed = design.repertoires.seed + repeat
    return list(StratifiedGroupKFold(n_splits=3, shuffle=True, random_state=seed).split(rows, strata, md["patient_id"].to_numpy()[rows]))


def baseline(design: StudyDesign, *, cv_column: Optional[str], rf, min_class_size: int, cv_repeats: int) -> BaselineResult:
    reps = design.repertoires
    md = reps.metadata
    if cv_column and cv_column not in md.columns:
        raise ValueError(f"cv_column {cv_column!r} is not a metadata column; see overview(), or set it to None")
    rf = rf if rf is not None else default_forest()
    rows = np.flatnonzero(design.keep & (md[cv_column].is_not_null().to_numpy() if cv_column else True))
    X = reps.X[rows]
    targets = [v for v in design.roles
               if design.kinds[v] == "categorical" and design.roles[v] != "covariate" and v != cv_column]
    fits, out = {}, []
    for v in targets:
        values = design.metadata[v].cast(pl.String).fill_null("unknown").to_numpy()[rows]
        levels, counts = np.unique(values, return_counts=True)
        classes = [c for c, k in zip(levels, counts) if min(k, len(values) - k) >= min_class_size and c != "unknown"]
        for c in classes[:1] if len(levels) == 2 else classes:  # a binary variable needs one row
            y = values == c
            splits = [_splits(design, rows, y, cv_column, r) for r in range(1 if cv_column else cv_repeats)]
            runs = [oof_scores(rf, X, y, s) for s in splits]
            with np.errstate(invalid="ignore"):
                scores = np.nanmean(runs, axis=0) if len(runs) > 1 else runs[0]
            fits[(v, c)] = Fit(y, splits[0], scores)
            out.append({"variable": v, "role": design.roles[v], "class": c, "n": len(y), "positives": int(y.sum()),
                        "auc": auc(y, scores), "ap": average_precision(y, scores)})
    table = pl.DataFrame(out).sort(pl.col("role").replace_strict(ROLE_ORDER), pl.col("auc"), descending=[False, True])
    return BaselineResult(table, rows, fits, rf, cv_column)


def within_batch(design: StudyDesign, base: BaselineResult, *, min_group_class_size: int) -> pl.DataFrame:
    """Per label class: AUC with all batches combined, and inside mixed-label batches."""
    if not design.batch:
        raise ValueError("no batch set")
    g = design.repertoires.metadata[design.batch].cast(pl.String).to_numpy()[base.rows]
    out = []
    for (v, c), fit in base.fits.items():
        if v != design.label:
            continue
        y = fit.y
        mixed = [b for b in np.unique(g) if min(y[g == b].sum(), (~y[g == b]).sum()) >= min_group_class_size]
        parts = [((g == b).sum(), auc(y[g == b], fit.scores[g == b])) for b in mixed]
        out.append({"class": c, "all batches": auc(y, fit.scores),
                    "within batch": float(np.average([a for _, a in parts], weights=[n for n, _ in parts])) if parts else None,
                    "mixed-label batches": len(mixed), "batches with class": len(np.unique(g[y]))})
    return pl.DataFrame(out)


@dataclass(frozen=True, repr=False)
class ShapResult(Result):
    """`table`: mean |SHAP| per feature over out-of-fold repertoires, largest first. The arrays hold the
    per-repertoire detail for `plots.shap_beeswarm`: a row per repertoire, a column per `table` row."""
    label: str
    target: str
    auc: float
    values: np.ndarray            # out-of-fold SHAP values
    data: np.ndarray              # the features' values
    others: np.ndarray            # per repertoire, the summed SHAP values of the features not in `table`
    n_others: int                 # how many features that sum covers
    in_target: np.ndarray         # whether the repertoire belongs to `target`
    batch: Optional[np.ndarray]   # the repertoire's batch; None if the design has none


def shap_importance(design: StudyDesign, base: BaselineResult, *, target: Optional[str], top: int) -> ShapResult:
    """Which features carry a label class: by default the most recoverable one."""
    import shap

    candidates = base.table.filter(pl.col("variable") == design.label)
    if target:
        candidates = candidates.filter(pl.col("class") == target)
    if not candidates.height:
        raise ValueError("no label class matches; see the baseline table")
    best = candidates.sort("auc", descending=True).row(0, named=True)
    fit = base.fits[(design.label, best["class"])]
    X, y = design.repertoires.X[base.rows], fit.y
    values = np.full(X.shape, np.nan)  # out-of-fold SHAP value per repertoire and feature
    for train, test in fit.splits:
        if 0 < y[train].sum() < len(train):
            model = clone(base.rf).fit(X[train], y[train])
            sv = shap.TreeExplainer(model, model_output="raw").shap_values(X[test])
            values[test] = sv[1] if isinstance(sv, list) else (sv[:, :, 1] if sv.ndim == 3 else sv)
    ok = ~np.isnan(values[:, 0])
    xc, sc = X[ok] - X[ok].mean(axis=0), values[ok] - values[ok].mean(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):  # a constant feature has no direction
        r = (xc * sc).sum(axis=0) / np.sqrt((xc ** 2).sum(axis=0) * (sc ** 2).sum(axis=0))
    # direction: does a higher value push towards the class (positive) or away from it (negative)?
    table = (pl.DataFrame({"feature": design.repertoires.feature_names, "mean_abs_shap": np.abs(values[ok]).mean(axis=0),
                           "value-SHAP correlation": r})
             .with_row_index("column")
             .with_columns(pl.when(pl.col("value-SHAP correlation") > 0).then(pl.lit("positive"))
                           .when(pl.col("value-SHAP correlation") < 0).then(pl.lit("negative")).otherwise(pl.lit("none"))
                           .alias("association"))
             .sort("mean_abs_shap", descending=True).head(top))
    cols = table["column"].to_numpy()
    rest = np.setdiff1d(np.arange(X.shape[1]), cols)
    batch = (design.repertoires.metadata[design.batch].cast(pl.String).fill_null("unknown").to_numpy()[base.rows][ok]
             if design.batch else None)
    return ShapResult(table.drop("column"), design.label, best["class"], best["auc"], values[ok][:, cols],
                      X[ok][:, cols], values[ok][:, rest].sum(axis=1), len(rest), y[ok], batch)
