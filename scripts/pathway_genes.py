"""The protein-coding genes of the primary micronutrient pathways.

The pathways are the rows of the ``PRIMARY pathways`` sheet of the micronutrient
workbook and of the keto diet workbook. Reactome genes come
from the g:Profiler GMT that the enrichment also uses, WikiPathways genes from
the WikiPathways GMT, and KEGG genes from the KEGG REST service. Every gene is
mapped to its HGNC symbol and kept only when HGNC gives it as a protein-coding
gene.
"""

from functools import lru_cache
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "input"
SELECTIONS = [DATA / "micronutrients" / "micronutrients_pathways.xlsx",
              DATA / "keto_diet" / "keto_diet_pathways.xlsx"]
SHEET = "PRIMARY pathways"
REACTOME_GMT = DATA / "gprofiler_gmt" / "hsapiens.REAC.name.gmt"
WIKIPATHWAYS_GMT = DATA / "wikipathways" / "wikipathways_Homo_sapiens.gmt"
KEGG_DIR = ROOT / "work" / "kegg"
HGNC = DATA / "hgnc" / "hgnc_complete_set.txt"

KEGG_URL = "https://rest.kegg.jp/link/hsa/{}"


@lru_cache(maxsize=None)
def hgnc_protein_coding():
    """The approved protein-coding HGNC genes, with symbol and Entrez ID."""
    table = pd.read_csv(HGNC, sep="\t", dtype=str,
                        usecols=["symbol", "locus_group", "status", "entrez_id"])
    return table[(table["status"] == "Approved")
                 & (table["locus_group"] == "protein-coding gene")]


def _read_gmt(path):
    """A GMT file as ``{first column: [genes]}``."""
    sets = {}
    for line in Path(path).read_text().splitlines():
        fields = line.split("\t")
        sets[fields[0]] = [g for g in fields[2:] if g]
    return sets


def reactome_symbols(accession):
    return _read_gmt(REACTOME_GMT)[f"REAC:{accession}"]


def wikipathways_entrez(accession):
    for name, genes in _read_gmt(WIKIPATHWAYS_GMT).items():
        if name.split("%")[2] == accession:
            return genes
    raise KeyError(f"{accession} not in {WIKIPATHWAYS_GMT.name}")


def kegg_entrez(accession):
    """The Entrez IDs of one KEGG pathway. The REST answer is kept in ``KEGG_DIR``."""
    path = KEGG_DIR / f"{accession}.tsv"
    if not path.exists():
        response = requests.get(KEGG_URL.format(accession), timeout=60)
        response.raise_for_status()
        KEGG_DIR.mkdir(parents=True, exist_ok=True)
        path.write_text(response.text)
    return [line.split("\t")[1].removeprefix("hsa:")
            for line in path.read_text().splitlines() if line]


def pathway_genes(database, accession):
    """``(protein-coding symbols, number of genes before the filter)``."""
    coding = hgnc_protein_coding()
    if database == "Reactome":
        genes = reactome_symbols(accession)
        kept = coding[coding["symbol"].isin(genes)]["symbol"]
    elif database == "WikiPathways":
        genes = wikipathways_entrez(accession)
        kept = coding[coding["entrez_id"].isin(genes)]["symbol"]
    elif database == "KEGG":
        genes = kegg_entrez(accession)
        kept = coding[coding["entrez_id"].isin(genes)]["symbol"]
    else:
        raise ValueError(f"unknown database {database!r}")
    return sorted(set(kept)), len(set(genes))


def _selection():
    """The primary pathways of all workbooks, one row per pathway."""
    return pd.concat([pd.read_excel(path, sheet_name=SHEET)
                      for path in SELECTIONS], ignore_index=True)


@lru_cache(maxsize=None)
def primary_pathway_genes():
    """One row per micronutrient, pathway and protein-coding gene."""
    selection = _selection()
    rows = []
    for _, pathway in selection.iterrows():
        symbols, _ = pathway_genes(pathway["Database"],
                                   pathway["Pathway accession"])
        rows += [{"Micronutrient": pathway["Micronutrient"],
                  "Database": pathway["Database"],
                  "Accession": pathway["Pathway accession"],
                  "Pathway": pathway["Primary pathway"],
                  "GeneSymbol": symbol} for symbol in symbols]
    return pd.DataFrame(rows)


def pathway_summary():
    """One row per primary pathway: gene count before and after the filter."""
    selection = _selection()
    rows = []
    for _, pathway in selection.iterrows():
        symbols, n_all = pathway_genes(pathway["Database"],
                                       pathway["Pathway accession"])
        rows.append({"Micronutrient": pathway["Micronutrient"],
                     "Database": pathway["Database"],
                     "Accession": pathway["Pathway accession"],
                     "Pathway": pathway["Primary pathway"],
                     "n_genes": n_all,
                     "n_protein_coding": len(symbols)})
    return pd.DataFrame(rows)
