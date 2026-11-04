"""MkDocs hook: publish every rendered notebook as a downloadable .ipynb under download/.

mkdocs-jupyter turns the notebooks into pages. This hook also adds each .ipynb itself as a plain file,
so a docs page can link to it with an ordinary relative link, e.g. `download/0_ingest_data.ipynb`, which
mkdocs checks and rewrites like any other link.
"""

from pathlib import PurePosixPath

from mkdocs.plugins import event_priority
from mkdocs.structure.files import File

DOWNLOAD_DIR = "download"


@event_priority(-100)  # run after mkdocs-jupyter, so the notebooks are already pages
def on_files(files, config):
    notebooks = [f for f in files if f.src_uri.endswith(".ipynb") and f.is_documentation_page()]
    for nb in notebooks:
        name = PurePosixPath(nb.src_uri).name
        files.append(File.generated(config, f"{DOWNLOAD_DIR}/{name}", abs_src_path=nb.abs_src_path))
    return files
