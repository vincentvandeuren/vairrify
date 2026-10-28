# Installation

## Requirements

- Python 3.10 or newer
- [tcrio](https://github.com/vincentvandeuren/tcr_qc), which reads and processes the repertoires. It is
  installed automatically, from PyPI.
- polars below 2.0: tcrio's compiled extension is tied to the polars version it was built against.

## pip

```bash
pip install vairrify                                              # TODO: once on PyPI
pip install git+https://github.com/vincentvandeuren/vairrify     # until then
```

## conda

The repository's `environment.yml` takes the scientific stack from conda-forge, and installs tcrio and
vairrify with pip:

```bash
git clone https://github.com/vincentvandeuren/vairrify
cd vairrify
conda env create -f environment.yml
conda activate vairrify
```

## For development

```bash
pip install -e ".[dev]"     # jupyterlab, nbconvert, pytest, ruff
pip install -e ".[docs]"    # mkdocs, to build these pages
```

## Compatibility notes

- TODO: tcrio dataset format versions. A dataset written by a newer tcrio can't be opened by an older
  one; say which tcrio version vairrify needs.
- Loading reruns a tcrio operation whose result is missing or out of date (`run_operations=True`).
  Opening a dataset with a different tcrio version than the one that processed it can therefore
  recompute results.
