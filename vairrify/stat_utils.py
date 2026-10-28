"""How strongly sample metadata shape a feature matrix, and how entangled they are.

Status: prototype. NumPy ports of tools that are usually reached through R:

- PERMANOVA -- ``vegan::adonis2``. One fit, read four ways (the ``test``
  column): each variable ``alone``; ``sequential`` (added in the given order,
  vegan's ``by = "terms"``); ``marginal`` (each adjusted for all others,
  ``by = "margin"``); and the whole ``model`` (``by = NULL``).
- Dispersion -- ``vegan::betadisper`` + ``permutest``: do groups differ in
  *spread* rather than in where their centre sits?
- Variation partitioning -- ``vegan::varpart``: unique and shared adjusted R²
  of 2-4 sets of variables.
- PVCA -- Bioconductor ``pvca::pvcaBatchAssess``, minus its ``vsn2``
  microarray normalisation.
- :func:`metadata_associations` -- how strongly the metadata variables are tied
  to *each other*, before looking at ``X`` at all.

Everything starts from :func:`prepare_design`, which applies one row mask to
``X`` and the metadata and records what it dropped. Rows of ``X`` and the
metadata must already match positionally; nothing here reorders them.

Preprocessing is the caller's job. With the default Euclidean distance every
feature counts in proportion to its variance, so mixed-unit feature blocks
usually want :func:`standardize` first -- after which a block with many
features outweighs a block with few.

Effect sizes from different tools answer different questions and do not add up
across tools. PERMANOVA R² are fractions of the total sum of squares of the
chosen distance; varpart fractions are *adjusted* R² of the same; PVCA
proportions are eigenvalue-weighted variance components of a mixed model fitted
to a few PCs of a sample-correlation matrix.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass
from typing import Callable, Mapping, Optional, Sequence, Union

import numpy as np
import polars as pl
from scipy import linalg, optimize, stats
from scipy.spatial.distance import pdist, squareform

__all__ = [
    "Design",
    "DispersionResult",
    "PvcaResult",
    "VarpartResult",
    "cramers_v",
    "dispersion",
    "metadata_associations",
    "permanova",
    "prepare_design",
    "pvca",
    "standardize",
    "variation_partitioning",
]

Distance = Union[str, Callable[[np.ndarray, np.ndarray], float]]
Permutations = Union[int, np.ndarray]

# vegan's slack when counting permuted statistics >= the observed one
_EPS = np.sqrt(np.finfo(float).eps)
# compared after stripping whitespace and lower-casing
_NA_STRINGS = ("", "na", "n/a", "nan", "none", "null", "-", "--", "?")


# --------------------------------------------------------------------------
# Design: one row mask for X and metadata
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class Design:
    """Rows of ``X`` and metadata that survive missing-value handling.

    Attributes:
        X: Kept feature rows, in original order.
        metadata: The analysed columns for the kept rows. Categorical columns
            are ``pl.String``, continuous ones ``pl.Float64``.
        continuous: Continuous variable names.
        categorical: Categorical variable names.
        kept: Boolean mask over the *original* rows.
        exclusions: Samples lost per reason. A sample missing several things is
            counted under each, so the reasons need not sum to the total.
        strata: Integer block codes for restricted permutations (e.g. patient
            when a patient contributes several samples), or ``None`` for free
            permutation.
    """

    X: np.ndarray
    metadata: pl.DataFrame
    continuous: tuple[str, ...]
    categorical: tuple[str, ...]
    kept: np.ndarray
    exclusions: pl.DataFrame
    strata: Optional[np.ndarray] = None

    @property
    def variables(self) -> tuple[str, ...]:
        return self.continuous + self.categorical

    @property
    def n(self) -> int:
        return self.X.shape[0]

    def with_X(self, X: np.ndarray) -> "Design":
        """Same samples, different features -- e.g. after :func:`standardize`."""
        X = np.asarray(X, dtype=float)
        if X.shape[0] != self.n:
            raise ValueError(f"X has {X.shape[0]} rows, design has {self.n}")
        return dataclasses.replace(self, X=X)


def prepare_design(
    X: np.ndarray,
    metadata: pl.DataFrame,
    *,
    continuous: Sequence[str] = (),
    categorical: Sequence[str] = (),
    strata: Optional[str] = None,
    missing: str = "drop",
    missing_label: str = "unknown",
    na_strings: Sequence[str] = _NA_STRINGS,
) -> Design:
    """Align ``X`` with the analysed metadata columns under one row mask.

    Args:
        X: Samples x features. Rows with any non-finite value are dropped.
        metadata: One row per row of ``X``, in the same order.
        continuous: Numeric variables. Null, NaN and +-inf count as missing.
        categorical: Grouping variables. Any dtype; values are compared as
            strings, so integer-coded folds or batches belong here. Null and
            the strings in ``na_strings`` count as missing.
        strata: Optional column whose levels are permutation blocks. Samples
            are only shuffled within a block. Missing values are dropped.
        missing: ``"drop"`` removes samples with a missing categorical value;
            ``"level"`` keeps them as a level named ``missing_label``.
            Missing continuous values are always dropped.
        missing_label: Level name used by ``missing="level"``.
        na_strings: Strings treated as missing in categorical (and string-typed
            continuous) columns, matched after stripping and lower-casing.

    Returns:
        The kept rows, plus an account of what was dropped and why.
    """
    X = np.asarray(X, dtype=float)
    if X.ndim != 2:
        raise ValueError(f"X must be 2-D, got shape {X.shape}")
    if X.shape[0] != metadata.height:
        raise ValueError(f"X has {X.shape[0]} rows but metadata has {metadata.height}")
    if missing not in ("drop", "level"):
        raise ValueError(f"missing must be 'drop' or 'level', got {missing!r}")
    continuous, categorical = tuple(continuous), tuple(categorical)
    if not continuous + categorical:
        raise ValueError("no variables given")
    if both := set(continuous) & set(categorical):
        raise ValueError(
            f"variables listed as both continuous and categorical: {sorted(both)}"
        )
    wanted = continuous + categorical + ((strata,) if strata else ())
    if absent := [c for c in wanted if c not in metadata.columns]:
        raise ValueError(f"columns not in metadata: {absent}")

    na_strings = tuple(s.strip().lower() for s in na_strings)

    def textual_na(s: pl.Series) -> np.ndarray:
        text = s.cast(pl.String)
        return (
            text.is_null() | text.str.strip_chars().str.to_lowercase().is_in(na_strings)
        ).to_numpy()

    columns: dict[str, pl.Series] = {}
    reasons: dict[str, np.ndarray] = {}
    for name in continuous:
        s = metadata[name]
        if s.dtype == pl.String:
            s = pl.Series(
                name, np.where(textual_na(s), None, s.to_numpy()), dtype=pl.String
            )
        values = s.cast(pl.Float64).fill_null(np.nan).to_numpy()
        reasons[f"{name} missing"] = ~np.isfinite(values)
        columns[name] = pl.Series(name, values, dtype=pl.Float64)
    for name in categorical:
        s = metadata[name].cast(pl.String)
        na = textual_na(metadata[name])
        if missing == "level":
            s = pl.Series(
                name, np.where(na, missing_label, s.to_numpy()), dtype=pl.String
            )
        else:
            reasons[f"{name} missing"] = na
        columns[name] = s
    strata_values = None
    if strata:
        strata_values = metadata[strata].cast(pl.String).to_numpy()
        reasons[f"{strata} (strata) missing"] = textual_na(metadata[strata])
    reasons["X not finite"] = ~np.isfinite(X).all(axis=1)

    dropped = np.zeros(X.shape[0], dtype=bool)
    for bad in reasons.values():
        dropped |= bad
    kept = ~dropped

    exclusions = pl.DataFrame(
        {
            "reason": list(reasons) + ["dropped (any reason)", "kept"],
            "n_samples": [int(b.sum()) for b in reasons.values()]
            + [int(dropped.sum()), int(kept.sum())],
        }
    )
    frame = pl.DataFrame([columns[c] for c in continuous + categorical]).filter(
        pl.Series(kept)
    )
    for name in categorical:
        if frame[name].n_unique() < 2:
            raise ValueError(
                f"{name!r} has fewer than two levels after dropping missing samples"
            )
    codes = None
    if strata_values is not None:
        codes = np.unique(strata_values[kept].astype(str), return_inverse=True)[1]
    return Design(
        X=X[kept],
        metadata=frame,
        continuous=continuous,
        categorical=categorical,
        kept=kept,
        exclusions=exclusions,
        strata=codes,
    )


def standardize(X: np.ndarray) -> np.ndarray:
    """Centre each feature and scale it to unit (population) variance.

    Matches ``sklearn.preprocessing.StandardScaler``: constant features become 0.
    """
    X = np.asarray(X, dtype=float)
    sd = X.std(axis=0)
    return (X - X.mean(axis=0)) / np.where(sd > 0, sd, 1.0)


# --------------------------------------------------------------------------
# Shared linear algebra
# --------------------------------------------------------------------------


def _codes(values) -> tuple[np.ndarray, np.ndarray]:
    """Sorted levels and the integer code of each value."""
    return np.unique(np.asarray(values).astype(str), return_inverse=True)


def _term_matrix(design: Design, name: str) -> np.ndarray:
    """Centred model-matrix columns of one variable (treatment dummies for factors)."""
    values = design.metadata[name].to_numpy()
    if name in design.continuous:
        v = values.astype(float)
        return (v - v.mean())[:, None]
    levels, codes = _codes(values)
    dummies = (codes[:, None] == np.arange(1, len(levels))[None, :]).astype(float)
    return dummies - dummies.mean(axis=0)


def _orthonormal(
    A: np.ndarray, against: Optional[np.ndarray] = None, tol: float = 1e-7
) -> np.ndarray:
    """Orthonormal basis of the part of span(A) not already in span(``against``).

    ``against`` must itself be orthonormal. Directions shorter than ``tol``
    times A's longest column count as aliased, as in R's ``qr()``.
    """
    n = A.shape[0]
    if A.shape[1] == 0:
        return np.zeros((n, 0))
    scale = np.linalg.norm(A, axis=0).max()
    if against is not None and against.shape[1]:
        A = A - against @ (against.T @ A)
    U, s, _ = np.linalg.svd(A, full_matrices=False)
    return U[:, s > tol * scale]


def _double_centre(M: np.ndarray) -> np.ndarray:
    return M - M.mean(axis=0, keepdims=True) - M.mean(axis=1, keepdims=True) + M.mean()


def _gower_coordinates(
    X: np.ndarray, distance: Distance
) -> tuple[np.ndarray, np.ndarray]:
    """Coordinates reproducing Gower's centred matrix of ``distance``.

    Returns ``(Y_pos, Y_neg)`` with ``G = -J D² J / 2 = Y_pos Y_posᵀ - Y_neg Y_negᵀ``.
    ``Y_neg`` has columns only for non-Euclidean distances, whose G has negative
    eigenvalues. Every sum of squares below is a trace of G, which these
    coordinates turn into ordinary squared norms.
    """
    n, p = X.shape
    if distance == "euclidean":
        Xc = X - X.mean(axis=0)
        if p < n:
            return Xc, np.zeros((n, 0))
        G = Xc @ Xc.T
    else:
        G = -0.5 * _double_centre(squareform(pdist(X, metric=distance)) ** 2)
    lam, U = np.linalg.eigh(G)
    tol = max(_EPS, _EPS * np.abs(lam).max())
    pos, neg = lam > tol, lam < -tol
    return U[:, pos] * np.sqrt(lam[pos]), U[:, neg] * np.sqrt(-lam[neg])


def _total_ss(coords: tuple[np.ndarray, np.ndarray]) -> float:
    """tr(G): total sum of squares."""
    pos, neg = coords
    return float(np.sum(pos**2) - np.sum(neg**2))


def _projected_ss(Q: np.ndarray, coords: tuple[np.ndarray, np.ndarray]) -> float:
    """tr(Qᵀ G Q): sum of squares of the data captured by the basis Q."""
    pos, neg = coords
    return float(np.sum((Q.T @ pos) ** 2) - np.sum((Q.T @ neg) ** 2))


def _permutation_matrix(
    n: int, permutations: Permutations, strata: Optional[np.ndarray], seed
) -> np.ndarray:
    """``(n_perm, n)`` row indices; row k means "sample i takes the data of sample perm[k, i]".

    An explicit matrix is passed through (0-based), which is how results are
    reproduced exactly against vegan given the same permutations.
    """
    if not np.isscalar(permutations):
        perms = np.asarray(permutations, dtype=np.intp)
        if perms.ndim != 2 or perms.shape[1] != n:
            raise ValueError(f"permutation matrix must have shape (n_perm, {n})")
        return perms
    rng = np.random.default_rng(seed)
    n_perm = int(permutations)
    if strata is None:
        return rng.permuted(np.tile(np.arange(n), (n_perm, 1)), axis=1)
    perms = np.tile(np.arange(n), (n_perm, 1))
    for block in np.unique(strata):
        idx = np.flatnonzero(strata == block)
        perms[:, idx] = idx[
            rng.permuted(np.tile(np.arange(len(idx)), (n_perm, 1)), axis=1)
        ]
    return perms


def _p_value(f_perm: np.ndarray, f_obs: np.ndarray) -> np.ndarray:
    """vegan's permutation p: (#{F* >= F - eps} + 1) / (n_perm + 1)."""
    with np.errstate(invalid="ignore"):
        p = ((f_perm >= f_obs - _EPS).sum(axis=0) + 1) / (f_perm.shape[0] + 1)
    return np.where(np.isfinite(f_obs), p, np.nan)


def _check_variables(design: Design, names: Sequence[str]) -> None:
    if unknown := [v for v in names if v not in design.variables]:
        raise ValueError(
            f"not in the design: {unknown}; design has {list(design.variables)}"
        )


# --------------------------------------------------------------------------
# Metadata vs metadata
# --------------------------------------------------------------------------


def cramers_v(x, y) -> float:
    """Bias-corrected Cramér's V (Bergsma 2013), so unrelated factors sit near 0.

    ``x`` and ``y`` are equal-length sequences of labels; values are compared as
    strings.
    """
    _, xi = _codes(x)
    _, yi = _codes(y)
    table = np.zeros((xi.max() + 1, yi.max() + 1))
    np.add.at(table, (xi, yi), 1)
    n = table.sum()
    r, k = table.shape
    expected = np.outer(table.sum(axis=1), table.sum(axis=0)) / n
    chi2 = ((table - expected) ** 2 / expected).sum()
    phi2 = max(0.0, chi2 / n - (k - 1) * (r - 1) / (n - 1))
    denom = min(k - (k - 1) ** 2 / (n - 1), r - (r - 1) ** 2 / (n - 1)) - 1
    return float(np.sqrt(phi2 / denom)) if denom > 0 else np.nan


def _correlation_ratio(groups, values: np.ndarray) -> float:
    """η: share of a continuous variable's spread explained by a factor, square-rooted."""
    _, codes = _codes(groups)
    means = np.bincount(codes, values) / np.bincount(codes)
    total = ((values - values.mean()) ** 2).sum()
    between = ((means[codes] - values.mean()) ** 2).sum()
    return float(np.sqrt(between / total)) if total > 0 else np.nan


def metadata_associations(
    design: Design, variables: Optional[Sequence[str]] = None
) -> pl.DataFrame:
    """How strongly each pair of metadata variables is tied, on a 0-1 scale.

    This is about the metadata only. When two variables are strongly tied, no
    analysis of ``X`` can say which of them a shared effect belongs to.

    Measures: bias-corrected Cramér's V (two factors), correlation ratio η
    (factor vs continuous), |Spearman ρ| (two continuous). They share a 0-1
    range but are different statistics; compare within a measure.

    Returns:
        Long frame with both orders of every pair and the diagonal (``1.0``):
        ``var_a``, ``var_b``, ``measure``, ``value``.
    """
    variables = design.variables if variables is None else tuple(variables)
    _check_variables(design, variables)
    md = design.metadata
    rows = []
    for a in variables:
        for b in variables:
            a_cat, b_cat = a in design.categorical, b in design.categorical
            if a == b:
                measure, value = "self", 1.0
            elif a_cat and b_cat:
                measure, value = (
                    "cramers_v",
                    cramers_v(md[a].to_numpy(), md[b].to_numpy()),
                )
            elif a_cat or b_cat:
                group, cont = (a, b) if a_cat else (b, a)
                measure, value = (
                    "eta",
                    _correlation_ratio(md[group].to_numpy(), md[cont].to_numpy()),
                )
            else:
                measure, value = (
                    "abs_spearman",
                    abs(float(stats.spearmanr(md[a], md[b]).statistic)),
                )
            rows.append({"var_a": a, "var_b": b, "measure": measure, "value": value})
    return pl.DataFrame(
        rows,
        schema={
            "var_a": pl.String,
            "var_b": pl.String,
            "measure": pl.String,
            "value": pl.Float64,
        },
    )


# --------------------------------------------------------------------------
# PERMANOVA
# --------------------------------------------------------------------------


def permanova(
    design: Design,
    *,
    terms: Optional[Sequence[str]] = None,
    distance: Distance = "euclidean",
    permutations: Permutations = 999,
    seed=None,
) -> pl.DataFrame:
    """PERMANOVA (``vegan::adonis2``) of the design's variables, read four ways.

    ``test`` says which question a row answers:

    - ``alone``: the variable as the only predictor. Large here but small
      under ``marginal`` means its effect is shared with other variables.
    - ``sequential``: added after the variables before it in ``terms``
      (``adonis2(by = "terms")``). Order matters whenever variables overlap.
    - ``marginal``: added last, after all the others
      (``adonis2(by = "margin")``). What it explains *beyond* everything else.
    - ``model``: all variables together (``adonis2(by = NULL)``), plus the
      ``residual`` row.

    Permutations shuffle whole samples (within strata if the design has them)
    and are shared by all rows, as in vegan. Like vegan, this tests location
    only under equal spread: groups that differ in spread (:func:`dispersion`)
    can produce a significant result too, especially with unequal group sizes.

    Args:
        design: From :func:`prepare_design`.
        terms: Variables in sequential order. Defaults to continuous then
            categorical, as given to :func:`prepare_design`.
        distance: ``"euclidean"`` (on ``design.X`` as given) or any
            ``scipy.spatial.distance.pdist`` metric.
        permutations: Number of permutations, or an explicit ``(n_perm, n)``
            matrix of 0-based row indices.
        seed: For the permutations.

    Returns:
        ``test``, ``term``, ``df``, ``sum_sq``, ``r2`` (of the total sum of
        squares), ``pseudo_F``, ``p_value``, ``n``.
    """
    terms = design.variables if terms is None else tuple(terms)
    _check_variables(design, terms)
    n = design.n
    coords = _gower_coordinates(design.X, distance)
    total = _total_ss(coords)
    blocks = [_term_matrix(design, t) for t in terms]

    sequential, Q = [], np.zeros((n, 0))
    for A in blocks:
        B = _orthonormal(A, Q)
        sequential.append(B)
        Q = np.hstack([Q, B])
    m = Q.shape[1]
    df_res = n - 1 - m
    if df_res <= 0:
        raise ValueError(
            f"{m} model degrees of freedom leave no residual with {n} samples"
        )
    marginal = []
    for i, A in enumerate(blocks):
        others = [b for j, b in enumerate(blocks) if j != i]
        rest = _orthonormal(np.hstack(others)) if others else np.zeros((n, 0))
        marginal.append(_orthonormal(A, rest))
    alone = [_orthonormal(A) for A in blocks]

    # Every basis lies inside the model space span(Q). Expressing each as
    # Q @ C means one m x m matrix S = Qᵀ G Q per permutation gives all of
    # their sums of squares as tr(Cᵀ S C).
    tests = (
        [("alone", t, B) for t, B in zip(terms, alone)]
        + [("sequential", t, B) for t, B in zip(terms, sequential)]
        + [("marginal", t, B) for t, B in zip(terms, marginal)]
        + [("model", "model", Q)]
    )
    coefs = [Q.T @ B for *_, B in tests]
    df = np.array([B.shape[1] for *_, B in tests], dtype=float)
    is_alone = np.array([kind == "alone" for kind, *_ in tests])

    def sums_of_squares(Qp: np.ndarray) -> np.ndarray:
        pos, neg = (Qp.T @ c for c in coords)
        S = pos @ pos.T - neg @ neg.T
        return np.array([np.sum(C * (S @ C)) for C in coefs])

    def f_stats(ss: np.ndarray) -> np.ndarray:
        with np.errstate(divide="ignore", invalid="ignore"):
            resid = np.where(
                is_alone, (total - ss) / (n - 1 - df), (total - ss[-1]) / df_res
            )
            f = (ss / df) / resid
        return np.where((df > 0) & (n - 1 - df > 0), f, np.nan)

    ss_obs = sums_of_squares(Q)
    f_obs = f_stats(ss_obs)
    perms = _permutation_matrix(n, permutations, design.strata, seed)
    # Qᵀ Y[perm] == Q[perm⁻¹]ᵀ Y: permuting Q's rows is the cheaper copy
    f_perm = np.array([f_stats(sums_of_squares(Q[np.argsort(p)])) for p in perms])
    p = _p_value(f_perm, f_obs)

    out = pl.DataFrame(
        {
            "test": [kind for kind, *_ in tests] + ["model"],
            "term": [term for _, term, _ in tests] + ["residual"],
            "df": np.append(df, df_res).astype(int),
            "sum_sq": np.append(ss_obs, total - ss_obs[-1]),
            "r2": np.append(ss_obs, total - ss_obs[-1]) / total,
            "pseudo_F": np.append(f_obs, np.nan),
            "p_value": np.append(p, np.nan),
        }
    )
    return out.with_columns(
        pl.lit(n).alias("n"), pl.col("pseudo_F", "p_value").fill_nan(None)
    )


# --------------------------------------------------------------------------
# Dispersion (PERMDISP)
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class DispersionResult:
    """Spread of each group around its own centre.

    Attributes:
        tests: Per variable, does mean distance-to-centre differ between groups:
            ``variable``, ``centre``, ``n_groups``, ``df_groups``,
            ``df_residual``, ``F``, ``p_value``.
        groups: ``variable``, ``level``, ``n``, ``mean_distance``.
        distances: Per sample and variable: ``sample`` (row in the design),
            ``variable``, ``level``, ``distance``.
    """

    tests: pl.DataFrame
    groups: pl.DataFrame
    distances: pl.DataFrame


def _spatial_median(
    P: np.ndarray, tol: float = 1e-10, max_iter: int = 10_000
) -> np.ndarray:
    """Point minimising the summed Euclidean distance to the rows of P (Weiszfeld)."""
    centre = np.median(P, axis=0)  # vegan's starting point
    for _ in range(max_iter):
        d = np.maximum(np.linalg.norm(P - centre, axis=1), 1e-12)
        new = (P / d[:, None]).sum(axis=0) / (1 / d).sum()
        if np.linalg.norm(new - centre) <= tol * (1 + np.linalg.norm(centre)):
            return new
        centre = new
    return centre


def _distances_to_centres(
    coords: tuple[np.ndarray, np.ndarray], codes: np.ndarray, centre: str
) -> np.ndarray:
    """betadisper's z_ij. Real and imaginary axes get separate centres, as in vegan."""
    sq = np.zeros(len(codes))
    for sign, Y in zip((1.0, -1.0), coords):
        if Y.shape[1] == 0:
            continue
        for g in np.unique(codes):
            idx = codes == g
            c = Y[idx].mean(axis=0) if centre == "centroid" else _spatial_median(Y[idx])
            sq[idx] += sign * ((Y[idx] - c) ** 2).sum(axis=1)
    return np.sqrt(np.maximum(sq, 0.0))


def _one_way_f(values: np.ndarray, onehot: np.ndarray) -> np.ndarray:
    """F of a one-way ANOVA along the last axis of ``values``; ``onehot`` is n x groups."""
    n, k = onehot.shape
    fitted = ((values @ onehot) / onehot.sum(axis=0)) @ onehot.T
    mss = ((fitted - values.mean(axis=-1, keepdims=True)) ** 2).sum(axis=-1)
    rss = ((values - fitted) ** 2).sum(axis=-1)
    return (mss / (k - 1)) / (rss / (n - k))


def dispersion(
    design: Design,
    *,
    variables: Optional[Sequence[str]] = None,
    distance: Distance = "euclidean",
    centre: str = "median",
    permutations: Permutations = 999,
    seed=None,
) -> DispersionResult:
    """Do groups differ in spread? (``vegan::betadisper`` + ``permutest``).

    Each sample's distance to its group's centre, then a one-way ANOVA of
    those distances across groups, with a permutation p-value from shuffled
    ANOVA residuals. A significant result means some groups are more
    scattered (e.g. noisier samples), which is a different finding from the
    shift in location PERMANOVA is designed to detect -- and one that can
    leak into a PERMANOVA result.

    Args:
        design: From :func:`prepare_design`.
        variables: Categorical variables to test; defaults to all of them.
        distance: As in :func:`permanova`.
        centre: ``"median"`` (spatial median, vegan's default) or ``"centroid"``.
        permutations: As in :func:`permanova`.
        seed: For the permutations.
    """
    variables = design.categorical if variables is None else tuple(variables)
    _check_variables(design, variables)
    if not_cat := [v for v in variables if v not in design.categorical]:
        raise ValueError(
            f"dispersion needs categorical variables, got continuous {not_cat}"
        )
    if centre not in ("median", "centroid"):
        raise ValueError(f"centre must be 'median' or 'centroid', got {centre!r}")
    n = design.n
    coords = _gower_coordinates(design.X, distance)
    perms = _permutation_matrix(n, permutations, design.strata, seed)

    tests, groups, distances = [], [], []
    for v in variables:
        levels, codes = _codes(design.metadata[v].to_numpy())
        k = len(levels)
        z = _distances_to_centres(coords, codes, centre)
        counts = np.bincount(codes, minlength=k)
        means = np.bincount(codes, z, minlength=k) / counts
        onehot = (codes[:, None] == np.arange(k)[None, :]).astype(float)
        f_obs = _one_way_f(z, onehot)
        # vegan permutes the ANOVA residuals, not the distances themselves
        f_perm = _one_way_f((z - means[codes])[perms], onehot)
        tests.append(
            {
                "variable": v,
                "centre": centre,
                "n_groups": k,
                "df_groups": k - 1,
                "df_residual": n - k,
                "F": float(f_obs),
                "p_value": float(_p_value(f_perm[:, None], np.array([f_obs]))[0]),
            }
        )
        groups.append(
            pl.DataFrame(
                {"variable": v, "level": levels, "n": counts, "mean_distance": means}
            )
        )
        distances.append(
            pl.DataFrame(
                {
                    "sample": np.arange(n),
                    "variable": v,
                    "level": levels[codes],
                    "distance": z,
                }
            )
        )
    return DispersionResult(
        tests=pl.DataFrame(tests),
        groups=pl.concat(groups),
        distances=pl.concat(distances),
    )


# --------------------------------------------------------------------------
# Variation partitioning
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class VarpartResult:
    """Unique and shared variation of 2-4 sets of variables (``vegan::varpart``).

    Attributes:
        fractions: Venn-diagram regions as adjusted R²: ``fraction`` (e.g.
            ``"study only"``, ``"study & disease"``, ``"residual"``),
            ``sets`` (list), ``n_sets``, ``adj_r2``. Shared fractions are
            not tests, can be slightly negative, and are where confounded
            variables put their common effect.
        models: The regressions they come from: ``sets``, ``label``, ``df``,
            ``r2``, ``adj_r2``.
    """

    fractions: pl.DataFrame
    models: pl.DataFrame


def variation_partitioning(
    design: Design,
    sets: Mapping[str, Sequence[str]],
    *,
    distance: Distance = "euclidean",
) -> VarpartResult:
    """Split explained variation into parts unique to, and shared between, sets.

    Fits the (distance-based) redundancy analysis of every union of sets,
    adjusts each R² for its degrees of freedom, and solves the Venn diagram
    by inclusion-exclusion. For two sets A and B: ``A only = R²(A+B) - R²(B)``,
    ``A & B = R²(A) + R²(B) - R²(A+B)``.

    Args:
        design: From :func:`prepare_design`.
        sets: Name -> variables, e.g. ``{"study": ["study_name"],
            "biology": ["disease", "age"]}``. 2-4 sets, no shared variables.
        distance: As in :func:`permanova`.
    """
    names = list(sets)
    if not 2 <= len(names) <= 4:
        raise ValueError(f"needs 2-4 sets, got {len(names)}")
    flat = [v for vs in sets.values() for v in vs]
    if len(flat) != len(set(flat)):
        raise ValueError("a variable appears in more than one set")
    _check_variables(design, flat)
    n, k = design.n, len(names)
    coords = _gower_coordinates(design.X, distance)
    total = _total_ss(coords)
    blocks = [
        np.hstack([_term_matrix(design, v) for v in sets[name]]) for name in names
    ]

    def members(mask: int) -> list[str]:
        return [names[i] for i in range(k) if mask >> i & 1]

    adj, models = {0: 0.0}, []
    for mask in range(1, 2**k):
        Q = _orthonormal(
            np.hstack([blocks[i] for i in range(k) if mask >> i & 1]), tol=1e-6
        )
        m = Q.shape[1]
        r2 = _projected_ss(Q, coords) / total
        adj[mask] = 1 - (1 - r2) * (n - 1) / (n - m - 1) if m < n - 1 else np.nan
        models.append(
            {
                "sets": members(mask),
                "label": " + ".join(members(mask)),
                "df": m,
                "r2": r2,
                "adj_r2": adj[mask],
            }
        )

    full = 2**k - 1
    rows = []
    for region in range(1, 2**k):
        # variation inside every set of `region` and outside all the others
        outside = full & ~region
        value = sum(
            (-1) ** sub.bit_count() * (adj[full] - adj[outside | sub])
            for sub in range(2**k)
            if sub & ~region == 0
        )
        inside = members(region)
        label = f"{inside[0]} only" if len(inside) == 1 else " & ".join(inside)
        rows.append(
            {"fraction": label, "sets": inside, "n_sets": len(inside), "adj_r2": value}
        )
    rows.append(
        {"fraction": "residual", "sets": [], "n_sets": 0, "adj_r2": 1 - adj[full]}
    )
    schema = {
        "fraction": pl.String,
        "sets": pl.List(pl.String),
        "n_sets": pl.Int64,
        "adj_r2": pl.Float64,
    }
    return VarpartResult(
        fractions=pl.DataFrame(rows, schema=schema), models=pl.DataFrame(models)
    )


# --------------------------------------------------------------------------
# PVCA
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class PvcaResult:
    """Principal variance component analysis.

    Attributes:
        summary: ``effect``, ``proportion``: eigenvalue-weighted average of
            each effect's share of variance over the modelled PCs. Sums to 1.
        per_pc: ``pc``, ``eigen_proportion``, ``effect``, ``variance``,
            ``proportion`` (share within that PC).
        n_pcs: Number of PCs modelled.
    """

    summary: pl.DataFrame
    per_pc: pl.DataFrame
    n_pcs: int


def _reml_random_intercepts(
    y: np.ndarray, groupings: Sequence[np.ndarray]
) -> tuple[np.ndarray, float]:
    """REML variances of ``y ~ 1 + (1|g1) + (1|g2) + ...``, fitted the way lme4 does.

    Minimises lme4's profiled REML deviance (Bates et al. 2015, J. Stat. Softw.
    67(1), eq. 41) over θ_k = σ_k / σ. Fixed effects, random effects and
    their normal equations collapse into one (q + 1)-square Cholesky per
    evaluation, so the cost does not grow with the number of samples.

    Returns:
        Variance of each grouping's random intercept, and the residual variance.
    """
    n = len(y)
    sizes, columns = [], []
    for g in groupings:
        _, codes = _codes(g)
        sizes.append(codes.max() + 1)
        columns.append(codes[:, None] == np.arange(sizes[-1])[None, :])
    W = np.hstack(columns + [np.ones((n, 1), dtype=bool)]).astype(float)  # [Z | X]
    WtW, Wty, yty = W.T @ W, W.T @ y, float(y @ y)
    q = W.shape[1] - 1
    term_of_column = np.repeat(np.arange(len(sizes)), sizes)
    diag_z = np.arange(q)

    def solve(theta: np.ndarray) -> tuple[float, float]:
        lam = np.append(theta[term_of_column], 1.0)
        M = WtW * lam[:, None] * lam[None, :]
        M[diag_z, diag_z] += 1.0
        rhs = lam * Wty
        L = linalg.cholesky(M, lower=True, check_finite=False)
        b = linalg.cho_solve((L, True), rhs, check_finite=False)
        # log|L_θ|² + log|R_X|², and the penalised residual sum of squares
        return 2.0 * np.log(np.diag(L)).sum(), max(yty - float(b @ rhs), 1e-300)

    def deviance(theta: np.ndarray) -> float:
        logdet, r2 = solve(theta)
        return logdet + (n - 1) * (1.0 + np.log(2.0 * np.pi * r2 / (n - 1)))

    # The deviance is even in every θ_k, so θ_k = 0 is always a stationary
    # point: a gradient method bounded at 0 that touches the bound stops
    # there (it did, on half the PCs of a test set). lme4 sidesteps this with
    # derivative-free BOBYQA; optimising over all reals and folding with |θ|
    # does the same job and matched or beat lme4's deviance on every PC tried.
    fit = optimize.minimize(
        deviance,
        np.ones(len(sizes)),
        method="L-BFGS-B",
        options={"ftol": 1e-15, "gtol": 1e-10, "maxiter": 10_000},
    )
    theta = np.abs(fit.x)
    sigma2 = solve(theta)[1] / (n - 1)
    return sigma2 * theta**2, sigma2


def _partition(labels: np.ndarray) -> tuple[int, ...]:
    """Canonical form of a grouping: equal for any two labellings of the same groups."""
    _, first, codes = np.unique(labels, return_index=True, return_inverse=True)
    order = np.argsort(np.argsort(first))  # number groups by first appearance
    return tuple(order[codes])


def pvca(
    design: Design,
    *,
    factors: Optional[Sequence[str]] = None,
    threshold: float = 0.6,
    interactions: bool = True,
) -> PvcaResult:
    """Principal variance component analysis (``pvca::pvcaBatchAssess``).

    Follows the Bioconductor implementation step by step, except that it does
    not run ``vsn2()`` first -- ``design.X`` is used as given:

    1. Centre each feature across samples, then take the eigenvectors of the
       sample x sample correlation matrix.
    2. Keep the fewest PCs whose cumulative eigenvalue share exceeds
       ``threshold`` (at least 3).
    3. Per PC, fit ``eigenvector ~ (1|f1) + (1|f2) + ... + (1|f1:f2) + ...`` by
       REML: every factor and every pairwise interaction a random intercept.
       Unlike pvca, an interaction that groups the samples exactly like an
       earlier effect (one factor nested in the other) is left out.
    4. Turn each PC's variance components into shares, weight by the PC's
       eigenvalue share, and average.

    Every factor is a random effect, so only categorical variables are allowed.
    (The R code would silently treat a continuous variable as one group per
    distinct value.)

    Args:
        design: From :func:`prepare_design`.
        factors: Categorical variables; defaults to all of them.
        threshold: Cumulative eigenvalue share the modelled PCs must cover.
        interactions: Include all pairwise interactions, as pvca always does.
    """
    factors = design.categorical if factors is None else tuple(factors)
    _check_variables(design, factors)
    if not_cat := [f for f in factors if f not in design.categorical]:
        raise ValueError(
            f"PVCA models factors as random effects; continuous variables not allowed: {not_cat}"
        )
    n = design.n
    md = design.metadata
    effects: dict[str, np.ndarray] = {f: md[f].to_numpy().astype(str) for f in factors}
    if interactions:
        # An interaction that groups samples exactly like an existing effect --
        # e.g. prep:experiment when prep is nested in experiment -- would be a
        # second copy of that random intercept, and REML would split its
        # variance between the two arbitrarily. Such interactions are skipped.
        partitions = {_partition(g) for g in effects.values()}
        for i, a in enumerate(factors):
            for b in factors[i + 1 :]:
                g = np.char.add(np.char.add(effects[a], "\x1f"), effects[b])
                if (p := _partition(g)) not in partitions:
                    effects[f"{a}:{b}"] = g
                    partitions.add(p)
    for name, g in effects.items():
        n_levels = len(np.unique(g))
        if not 1 < n_levels < n:  # lme4 refuses these too
            raise ValueError(
                f"{name!r} has {n_levels} levels for {n} samples; needs between 2 and n - 1"
            )

    centred = design.X - design.X.mean(axis=0)
    eigenvalues, eigenvectors = np.linalg.eigh(np.corrcoef(centred))
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues, eigenvectors = eigenvalues[order], eigenvectors[:, order]
    shares = eigenvalues / eigenvalues.sum()

    # pvca's own loop, kept verbatim so floating-point ties resolve the same way
    count, remaining = 0, 1.0
    for share in shares[::-1]:
        remaining -= share
        if remaining <= threshold:
            count += 1
    n_pcs = min(max(count, 3), n)

    names = list(effects) + ["residual"]
    variances = np.zeros((n_pcs, len(names)))
    for i in range(n_pcs):
        components, resid = _reml_random_intercepts(
            eigenvectors[:, i], list(effects.values())
        )
        variances[i] = np.append(components, resid)
    within = variances / variances.sum(axis=1, keepdims=True)
    weighted = (within * shares[:n_pcs, None]).sum(axis=0)

    per_pc = pl.DataFrame(
        {
            "pc": np.repeat(np.arange(1, n_pcs + 1), len(names)),
            "eigen_proportion": np.repeat(shares[:n_pcs], len(names)),
            "effect": names * n_pcs,
            "variance": variances.ravel(),
            "proportion": within.ravel(),
        }
    )
    summary = pl.DataFrame({"effect": names, "proportion": weighted / weighted.sum()})
    return PvcaResult(summary=summary, per_pc=per_pc, n_pcs=n_pcs)
