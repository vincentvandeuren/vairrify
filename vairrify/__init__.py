"""vAIRRify: audit a repertoire dataset for technical variation a model could use as a shortcut.

Follows the first part of the manuscript. Every analysis uses the features chosen at load (one
feature block, or several side by side); each step returns a table and `plots` draws it.

    import vairrify as vf
    vf.feature_blocks()                                                                  # the blocks to choose from
    reps = vf.load(path, locus="TRB", features="V genes")                                # 1-3
    qc = reps.qc()                                                                       # 4
    design = reps.design(label=..., batch=..., covariates=[...], exclude=qc.failed)      # 5-6
    design.balance(); design.pca(); design.variance()                                    # 7-9: Fig. S1, 2
    base = design.baseline(); design.within_batch(base); design.shap(base)               # 10-11: Table 1, Fig. 3
"""

from . import plots
from .design import StudyDesign
from .features import FEATURE_BLOCKS, FeatureBlock, feature_blocks
from .repertoires import Repertoires
from .results import save

load = Repertoires.load

__all__ = ["FEATURE_BLOCKS", "FeatureBlock", "Repertoires", "StudyDesign", "feature_blocks", "load", "plots", "save"]
