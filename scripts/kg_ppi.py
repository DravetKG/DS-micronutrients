"""The HPA brain gene list, cached in ``work/brain_genes.txt``."""

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
BRAIN_GENES_TXT = ROOT / "work" / "brain_genes.txt"
HPA_BRAIN_TSV = ROOT / "input" / "hpa" / "rna_brain_region_hpa.tsv"
NTPM_THRESHOLD = 1


def brain_genes(verbose=True):
    """Gene symbols with HPA brain-region nTPM above ``NTPM_THRESHOLD``.

    Read from ``work/brain_genes.txt`` when present, and recomputed from HPA
    when not.
    """
    if BRAIN_GENES_TXT.exists():
        genes = {g for g in BRAIN_GENES_TXT.read_text().split("\n") if g}
        if verbose:
            print(f"brain genes from {BRAIN_GENES_TXT.name}: {len(genes):,}")
        return genes

    hpa = pd.read_csv(HPA_BRAIN_TSV, sep="\t")
    expr = hpa.pivot_table(index="Gene name", columns="Brain region",
                           values="nTPM", aggfunc="max").fillna(0.0)
    genes = set(expr.index[(expr > NTPM_THRESHOLD).any(axis=1)])
    BRAIN_GENES_TXT.parent.mkdir(parents=True, exist_ok=True)
    BRAIN_GENES_TXT.write_text("\n".join(sorted(genes)))
    if verbose:
        print(f"brain genes recomputed from {HPA_BRAIN_TSV.name}: {len(genes):,}")
    return genes
