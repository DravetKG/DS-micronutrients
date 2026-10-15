"""The HPA kidney filter.

A gene is kidney-expressed when its HPA consensus nTPM in kidney is above
``NTPM_THRESHOLD``, the threshold of the brain filter.
"""

from functools import lru_cache
from pathlib import Path

import pandas as pd

from brain import NTPM_THRESHOLD

DATA = Path(__file__).resolve().parent.parent / "input"
HPA_TISSUE_TSV = DATA / "hpa" / "rna_tissue_consensus.tsv"


@lru_cache(maxsize=None)
def kidney_genes(ntpm_threshold=NTPM_THRESHOLD):
    """Cached HPA symbols expressed in kidney."""
    hpa = pd.read_csv(HPA_TISSUE_TSV, sep="\t")
    kidney = hpa[(hpa["Tissue"] == "kidney") & (hpa["nTPM"] > ntpm_threshold)]
    return frozenset(kidney["Gene name"])
