"""The HPA brain filter.

A gene is brain-expressed when its nTPM is above ``NTPM_THRESHOLD`` in at least
one of the 13 HPA brain regions.
"""

from functools import lru_cache
from pathlib import Path

import pandas as pd

DATA = Path(__file__).resolve().parent.parent / "input"

NTPM_THRESHOLD = 1


def _expressed_genes(df, group_col, value_col, groups, threshold, verbose, label):
    """Pivot and threshold an HPA table. Returns ``(genes, group_x_gene_matrix)``."""
    df = df[df[group_col].isin(groups)]
    expr = (
        df.pivot_table(index="Gene name", columns=group_col, values=value_col,
                        aggfunc="max")
        .fillna(0.0)
    )
    genes = set(expr.index[(expr > threshold).any(axis=1)])

    if verbose:
        print(f"{label} ({len(groups)}): {list(groups)}")
        print(f"Genes with data in these {group_col.lower()}s : {len(expr):,}")
        print(f"{label}-expressed ({value_col} > {threshold} in any) : {len(genes):,}")

    return genes, expr


def brain_expressed_genes(hpa_tsv, ntpm_threshold=1.0, verbose=True):
    """Genes above ``ntpm_threshold`` in at least one of the 13 brain regions."""
    hpa = pd.read_csv(hpa_tsv, sep="\t")
    regions = sorted(hpa["Brain region"].unique())
    return _expressed_genes(hpa, "Brain region", "nTPM", regions,
                             ntpm_threshold, verbose, "Brain region")


@lru_cache(maxsize=None)
def brain_genes(ntpm_threshold=NTPM_THRESHOLD):
    """Cached HPA symbols expressed in at least one brain region."""
    genes, _ = brain_expressed_genes(DATA / "hpa" / "rna_brain_region_hpa.tsv",
                                     ntpm_threshold)
    return frozenset(genes)
