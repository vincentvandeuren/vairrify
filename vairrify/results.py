"""What the steps return: a main table, plus whatever the plots need."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import polars as pl


@dataclass(frozen=True)
class Result:
    """A step's result. Shows as its `table` in a notebook; the other fields feed the plots."""

    table: pl.DataFrame

    def __repr__(self) -> str:
        return repr(self.table)

    def _repr_html_(self) -> str:
        return self.table._repr_html_()


def save(out_dir: str | Path, **results: Result | pl.DataFrame) -> list[Path]:
    """Write each result's table as `<name>.tsv`; e.g. `save(dir, qc=qc, baseline=base)`."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for name, result in results.items():
        table = result.table if isinstance(result, Result) else result
        path = out_dir / f"{name}.tsv"
        table.drop(pl.selectors.by_dtype(pl.List)).write_csv(path, separator="\t")
        written.append(path)
    return written
