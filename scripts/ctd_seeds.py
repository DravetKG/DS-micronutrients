"""The seed gene sets of the micronutrient colocalization.

Each micronutrient set is the union of its CTD genes and the genes of its
primary pathways, cut to the HGNC protein-coding genes. The CTD genes are read
from the workbook that ``ctd_micronutrients.ipynb`` writes, which already drops
the co-treatment records and applies one reference threshold per set. The pathway genes come
from ``pathway_genes``. The keto diet has no CTD genes, so its set is only its
primary pathway genes.
"""

from functools import lru_cache
from pathlib import Path

import pandas as pd

import pathway_genes

ROOT = Path(__file__).resolve().parent.parent
WORKBOOK = ROOT / "results" / "ctd" / "ctd_micronutrient_sets.xlsx"
SHEET = "Genes_filtered"

# Gene sets given as a plain symbol list rather than by the CTD workbook. A CSV
# is read by its ``label`` column and a text file by its lines.
SYMBOL_FILES = {}

KEYS = {
    "vita": "Vitamin A",
    "vitd": "Vitamin D",
    "vitc": "Vitamin C",
    "vite": "Vitamin E",
    "vitb6": "Vitamin B6",
    "vitb12": "Vitamin B12",
    "iron": "Iron",
    "folate": "Folic Acid",
    "zinc": "Zinc",
    "selenium": "Selenium",
    "keto": "Keto diet",
}

LABELS = {
    "vita": "vitamin A",
    "vitd": "vitamin D",
    "vitc": "vitamin C",
    "vite": "vitamin E",
    "vitb6": "vitamin B6",
    "vitb12": "vitamin B12",
    "iron": "iron",
    "folate": "folic acid",
    "zinc": "zinc",
    "selenium": "selenium",
    "keto": "keto diet",
}


@lru_cache(maxsize=None)
def _filtered_genes(workbook=None):
    """The sheet of genes that the notebook kept, one row per set and gene."""
    path = Path(workbook) if workbook else WORKBOOK
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found; run ctd_micronutrients.ipynb first"
        )
    return pd.read_excel(path, sheet_name=SHEET)


@lru_cache(maxsize=None)
def symbol_file_genes(key, brain_only=False):
    """The gene symbols of one entry of ``SYMBOL_FILES``."""
    from brain import NTPM_THRESHOLD, brain_genes

    path = SYMBOL_FILES[key]
    if path.suffix == ".csv":
        labels = pd.read_csv(path)["label"]
    else:
        labels = pd.Series(path.read_text().split())
    labels = labels.dropna().astype(str).str.strip()
    symbols = {label for label in labels if label and " " not in label}
    if brain_only:
        symbols &= set(brain_genes(NTPM_THRESHOLD))
    return tuple(sorted(symbols))


def ctd_genes(key, workbook=None):
    """``(CTD symbols, n_references_used)`` of one micronutrient."""
    table = _filtered_genes(workbook)
    table = table[table["Micronutrient"] == KEYS[key]]
    threshold = int(table["RefThreshold"].iloc[0]) if len(table) else None
    return set(table["GeneSymbol"].astype(str)), threshold


def primary_pathway_genes(key):
    """The protein-coding genes of the primary pathways of one micronutrient."""
    table = pathway_genes.primary_pathway_genes()
    return set(table.loc[table["Micronutrient"].str.lower() == KEYS[key].lower(),
                         "GeneSymbol"])


def genes_for(key, workbook=None, brain_only=False):
    """``(symbols, n_references_used)``; a ``SYMBOL_FILES`` set reports ``None``.

    A micronutrient set is the union of its CTD and primary pathway genes, cut
    to the HGNC protein-coding genes.
    """
    if key in SYMBOL_FILES:
        return list(symbol_file_genes(key, brain_only)), None

    ctd, threshold = ctd_genes(key, workbook)
    symbols = ((ctd | primary_pathway_genes(key))
               & set(pathway_genes.hgnc_protein_coding()["symbol"]))
    if brain_only:
        from brain import NTPM_THRESHOLD, brain_genes
        symbols &= set(brain_genes(NTPM_THRESHOLD))
    return sorted(symbols), threshold


def provenance_table(workbook=None, keys=None):
    """One row per micronutrient and gene, with the source of the gene."""
    pathways = pathway_genes.primary_pathway_genes()
    rows = []
    for key in keys or KEYS:
        if key in SYMBOL_FILES:
            continue
        ctd, _ = ctd_genes(key, workbook)
        mine = pathways[pathways["Micronutrient"].str.lower() == KEYS[key].lower()]
        accessions = mine.groupby("GeneSymbol")["Accession"].agg(
            lambda v: "|".join(sorted(v)))
        for gene in sorted(ctd | set(accessions.index)):
            rows.append({"key": key, "micronutrient": LABELS[key],
                         "GeneSymbol": gene, "in_ctd": gene in ctd,
                         "pathways": accessions.get(gene, "")})
    return pd.DataFrame(rows)


def seed_table(workbook=None, brain_only=False, keys=None):
    """One row per gene set: gene count by source and the CTD threshold it used."""
    rows = []
    for key in keys or LABELS:
        symbols, threshold = genes_for(key, workbook, brain_only)
        row = {"key": key, "micronutrient": LABELS[key],
               "min_references_used": threshold}
        if key not in SYMBOL_FILES:
            ctd, _ = ctd_genes(key, workbook)
            pathway = primary_pathway_genes(key)
            row.update(n_ctd=len(ctd), n_pathway=len(pathway),
                       n_ctd_and_pathway=len(ctd & pathway))
        row["n_genes"] = len(symbols)
        rows.append(row)
    table = pd.DataFrame(rows)
    counts = [c for c in table if c.startswith("n_") or c == "min_references_used"]
    return table.astype({c: "Int64" for c in counts})
