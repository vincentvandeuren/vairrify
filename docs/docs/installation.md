# Installation

## Requirements

- Python 3.10 or newer
- [tcrio](https://github.com/vincentvandeuren/tcr_qc), which reads and processes the repertoires. It is
  installed automatically, from PyPI.
- polars (below 2 i guess, I think tcrio pins that, but I might want to remove that pin #TODO)

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