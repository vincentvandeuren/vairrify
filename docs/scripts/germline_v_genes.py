"""Which V genes have non-standard usage across people, as germline variation would cause?

The evidence on the docs page "Germline variation". Reads V-gene usage from tcrio datasets
(GeneCountsSummary results; nothing is recomputed or written), keeps one repertoire per patient with
enough clonotypes, and compares each gene's usage to the dataset median:

- absent: below 0.05x the median, as a deletion on both chromosomes would give
- modes: peaks in a smoothed density of log2(usage / median) over the other patients

A gene is non-standard in a dataset when at least 2% of patients lack it, or its usage has at least two
modes that each hold 3% of patients. It is highlighted when non-standard in at least half the datasets
that assign it.

    python docs/scripts/germline_v_genes.py /data/current/tcrio_db docs/docs/img/germline_v_genes.png
"""

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import seaborn as sns
from scipy.signal import find_peaks
from scipy.stats import gaussian_kde, spearmanr
from tcr_io import operations as ops
from tcr_io.dataset import Dataset

DATASETS = {  # name: minimum clonotypes per repertoire
    "emerson_2017": 20_000,
    "delmonte_2024": 20_000,
    "gittelman_2022": 20_000,
    "rawat_2026": 20_000,
    "li_2025": 20_000,
    "karakus_2025": 20_000,
}
COLOURS = dict(zip(DATASETS, ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]))  # validated order
ABSENT, MIN_ABSENT, MIN_MODE, PROMINENCE, BANDWIDTH = 0.05, 0.02, 0.03, 0.15, 0.15
GRID = np.linspace(-5, 3, 801)  # log2(usage / median)
LINKED = ["TRBV6-2", "TRBV7-2", "TRBV10-2", "TRBV11-3", "TRBV12-5"]
NEIGHBOURS = [("TRBV10-2", "TRBV10-3"), ("TRBV11-3", "TRBV11-2"), ("TRBV12-5", "TRBV12-4"), ("TRBV28", "TRBV24-1")]


def usage(path: Path, min_clonotypes: int) -> pl.DataFrame:
    """V-gene usage, one repertoire per patient with at least `min_clonotypes` clonotypes."""
    ds = Dataset(path).select_locus("TRB")
    v = ds.get_result(ops.GeneCountsSummary.v_counts).read()
    return (ds.repertoire_meta.select("repertoire_id", "patient_id")
            .join(ds.repertoire_counts.select("repertoire_id", "n_clonotypes"), on="repertoire_id")
            .join(v, on="repertoire_id")
            .filter(pl.col("n_clonotypes") >= min_clonotypes)
            .unique("patient_id", keep="first", maintain_order=True)
            .select(v.columns[1:]))


def n_modes(r: np.ndarray) -> int:
    """Modes of log2(usage / median) among patients who use the gene, each holding >= MIN_MODE of all patients."""
    x = np.log2(r[r >= ABSENT])
    d = gaussian_kde(x, bw_method=BANDWIDTH / x.std(ddof=1))(GRID)
    peaks, _ = find_peaks(d, prominence=PROMINENCE * d.max())
    if len(peaks) < 2:
        return len(peaks)
    edges = [GRID[0]] + [GRID[a + np.argmin(d[a:b])] for a, b in zip(peaks[:-1], peaks[1:])] + [GRID[-1]]
    return int((np.histogram(x, bins=edges)[0] / len(r) >= MIN_MODE).sum())


def main(root: str, figure: str) -> None:
    ratios, rows = {}, []
    for name, min_clonotypes in DATASETS.items():
        u = usage(Path(root) / name, min_clonotypes)
        print(f"{name}: {u.height} patients")
        for g in u.columns:
            v = u[g].to_numpy()
            if np.median(v) == 0:  # not assigned, e.g. TRBV6-2/6-3 in some Adaptive exports
                continue
            r = v / np.median(v)
            ratios[(name, g)] = r
            absent, modes = float((r < ABSENT).mean()), n_modes(r)
            rows.append({"gene": g, "dataset": name, "absent": absent, "modes": modes,
                         "non-standard": absent >= MIN_ABSENT or modes >= 2})
        # genes tied to the 21-kb deletion: usage where TRBV4-3 is absent, relative to the median
        absent43 = ratios[(name, "TRBV4-3")] < ABSENT
        linked = {g: np.median(ratios[(name, g)][absent43]) for g in LINKED if (name, g) in ratios}
        # an allele assigned to a similar gene would make the two trade usage (negative correlation)
        neighbours = {f"{a}~{b}": spearmanr(ratios[(name, a)], ratios[(name, b)])[0]
                      for a, b in NEIGHBOURS if (name, a) in ratios and (name, b) in ratios}
        print("   where TRBV4-3 is absent:", ", ".join(f"{g} {x:.2f}x" for g, x in linked.items()))
        print("   neighbours (Spearman):", ", ".join(f"{k} {x:+.2f}" for k, x in neighbours.items()))
    per = pl.DataFrame(rows)
    genes = list(dict.fromkeys(per["gene"]))  # tcrio's gene order
    summary = (per.group_by("gene", maintain_order=True)
               .agg(pl.len().alias("datasets"), pl.col("non-standard").sum().alias("non-standard in"),
                    pl.col("absent").min().alias("absent min"), pl.col("absent").max().alias("absent max"),
                    pl.col("modes").median().alias("modes (median)"))
               .with_columns((pl.col("non-standard in") >= pl.col("datasets") / 2).alias("highlighted")))
    pl.Config.set_tbl_rows(60)
    print(summary.filter(pl.col("highlighted") | (pl.col("non-standard in") > 0)).with_columns(pl.col(pl.Float64).round(3)))
    highlighted = set(summary.filter(pl.col("highlighted"))["gene"])
    print("highlighted:", [g for g in genes if g in highlighted])

    sns.set_theme(style="whitegrid", context="notebook")
    bins = np.arange(0, 2.05, 0.05)
    ncols = 6
    nrows = -(-len(genes) // ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(3 * ncols, 2.4 * nrows), sharex=True)
    for ax, g in zip(axes.flat, genes):
        for d in DATASETS:
            if (d, g) in ratios:
                h, _ = np.histogram(np.clip(ratios[(d, g)], 0, 2.0), bins=bins)
                ax.stairs(h / h.sum(), bins, color=COLOURS[d], linewidth=1.3)
        if g in highlighted:
            ax.set_facecolor("#f4ede2")
        ax.set_title(g, fontsize=10, fontweight="bold" if g in highlighted else "normal")
        ax.tick_params(labelsize=8, labelbottom=True)
    for ax in axes.flat[len(genes):]:
        ax.set_visible(False)
    for col in range(ncols):  # the lowest visible panel of each column
        [ax for ax in axes[:, col] if ax.get_visible()][-1].set_xlabel("usage / dataset median", fontsize=9)
    handles = [plt.Line2D([], [], color=COLOURS[d], linewidth=2, label=d.replace("_", " ")) for d in DATASETS]
    handles.append(plt.Rectangle((0, 0), 1, 1, facecolor="#f4ede2", edgecolor="0.7", label="non-standard in at least half the datasets"))
    fig.legend(handles=handles, loc="upper center", ncol=4, bbox_to_anchor=(0.5, 1.0), fontsize=10, frameon=False)
    fig.suptitle("V-gene usage per patient, relative to the dataset median (2 = 2 or more)", y=1.025)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(figure, dpi=70, bbox_inches="tight")
    per.write_csv(Path(figure).with_suffix(".tsv"), separator="\t")
    print("saved", figure)


if __name__ == "__main__":
    main(*sys.argv[1:3])
