"""Build the micronutrient gene sets from the CTD chemical-gene interactions.

Gets the curated human interactions of eight micronutrients by hierarchical
MeSH search and makes three nested sets: all, no_cotreatment, and brain_top.
"""

import argparse
import io
import time
from pathlib import Path

import pandas as pd
import requests

from brain import NTPM_THRESHOLD, brain_genes

BATCH_URL = "https://ctdbase.org/tools/batchQuery.go"

MICRONUTRIENTS = {
    "Vitamin A": "D014801",
    "Vitamin D": "D014807",
    "Vitamin C": "D001205",
    "Vitamin E": "D014810",
    "Vitamin B6": "D025101",
    "Vitamin B12": "D014805",
    "Iron": "D007501",
    "Folic Acid": "D005492",
}

ORGANISM = "Homo sapiens"
MIN_REFERENCES = 1
START_REFERENCES = 2
MIN_SET_SIZE = 6
MAX_SET_SIZE = 200

INTERACTION_COLUMNS = [
    "Micronutrient",
    "ChemicalName",
    "ChemicalID",
    "CasRN",
    "GeneSymbol",
    "GeneID",
    "Organism",
    "OrganismID",
    "Interaction",
    "InteractionActions",
    "PubMedIDs",
]


def fetch_interactions(mesh_id, search_type="hierarchicalAssociations",
                       timeout=60):
    """Get the curated chemical-gene interactions of one chemical as a table."""
    response = requests.get(
        BATCH_URL,
        params={
            "inputType": "chem",
            "inputTerms": mesh_id,
            "report": "cgixns",
            "format": "tsv",
            "actionTypes": "ANY",
            "inputTermSearchType": search_type,
        },
        timeout=timeout,
    )
    response.raise_for_status()
    text = response.text.lstrip()
    if text.startswith("#"):
        text = text[1:]
    table = pd.read_csv(io.StringIO(text), sep="\t", dtype=str).fillna("")
    table.columns = [column.strip() for column in table.columns]
    return table


def fetch_all(micronutrients=None, pause=2.0):
    """Get the human interaction records of every micronutrient as one table."""
    micronutrients = micronutrients or MICRONUTRIENTS
    tables = []
    for index, (name, mesh_id) in enumerate(micronutrients.items()):
        raw = fetch_interactions(mesh_id)
        human = raw[raw["Organism"] == ORGANISM].copy()
        human["Micronutrient"] = name
        tables.append(human[INTERACTION_COLUMNS])
        print(f"{name} ({mesh_id}): {len(raw):>6d} records, "
              f"{len(human):>6d} human, "
              f"{human['GeneSymbol'].nunique():>5d} genes")
        if index < len(micronutrients) - 1:
            time.sleep(pause)
    return pd.concat(tables, ignore_index=True)


def is_cotreatment(interactions):
    """Mark the records that report the micronutrient with another chemical."""
    return interactions["Interaction"].str.contains("co-treated", case=False,
                                                    na=False)


def split_references(values):
    """Get the set of distinct PubMed identifiers of a series of records."""
    return {
        item.strip()
        for value in values
        for item in str(value).split("|")
        if item.strip()
    }


def gene_table(interactions, min_references=MIN_REFERENCES):
    """Collapse the interaction records to one row per micronutrient and gene."""
    genes = interactions.groupby(
        ["Micronutrient", "GeneSymbol", "GeneID"], as_index=False
    ).agg(
        NumInteractions=("Interaction", "size"),
        NumReferences=("PubMedIDs", lambda values: len(split_references(values))),
        PubMedIDs=("PubMedIDs",
                   lambda values: "|".join(sorted(split_references(values)))),
        ChemicalNames=("ChemicalName",
                       lambda values: "|".join(sorted(set(values)))),
        InteractionActions=("InteractionActions", lambda values: "|".join(
            sorted({action for value in values
                    for action in str(value).split("|") if action}))),
    )
    genes = genes[genes["NumReferences"] >= min_references]
    genes["BrainExpressed"] = genes["GeneSymbol"].isin(brain_genes(NTPM_THRESHOLD))
    return genes.sort_values(
        ["Micronutrient", "NumReferences", "GeneSymbol"],
        ascending=[True, False, True],
    ).reset_index(drop=True)


def reference_threshold(counts, start=START_REFERENCES,
                        min_size=MIN_SET_SIZE, max_size=MAX_SET_SIZE):
    """Choose one set's reference threshold, keeping it between the size bounds."""
    counts = pd.Series(counts)
    size = lambda threshold: int((counts >= threshold).sum())

    if size(start) <= min_size:
        return 1

    threshold = start
    while size(threshold) >= max_size:
        if size(threshold + 1) <= min_size:
            break
        threshold += 1
    return threshold


def top_brain_genes(genes, start=START_REFERENCES, min_size=MIN_SET_SIZE,
                    max_size=MAX_SET_SIZE):
    """Keep the brain-expressed genes and cut each set by its own threshold."""
    brain = genes[genes["BrainExpressed"]].copy()
    kept = []
    for _, group in brain.groupby("Micronutrient", sort=False):
        threshold = reference_threshold(group["NumReferences"], start,
                                        min_size, max_size)
        group = group[group["NumReferences"] >= threshold].copy()
        group["RefThreshold"] = threshold
        group["Rank"] = group["NumReferences"].rank(
            method="first", ascending=False
        ).astype(int)
        kept.append(group[group["Rank"] <= max_size])
    return pd.concat(kept, ignore_index=True)




def overview(genes, micronutrients=None):
    """One row per micronutrient: the gene count and the reference range."""
    micronutrients = micronutrients or MICRONUTRIENTS
    rows = []
    for name in micronutrients:
        group = genes[genes["Micronutrient"] == name]
        rows.append({
            "Micronutrient": name,
            "MeshID": micronutrients[name],
            "NumGenes": len(group),
            "NumBrainGenes": int(group["BrainExpressed"].sum()),
            "MaxReferences": int(group["NumReferences"].max()) if len(group) else 0,
            "MinReferences": int(group["NumReferences"].min()) if len(group) else 0,
        })
    return pd.DataFrame(rows)






def write_workbook(path, tables):
    """Write one Excel workbook, one sheet per entry of ``tables``."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        for sheet, table in tables.items():
            table.to_excel(writer, sheet_name=sheet[:31], index=False)
    print(f"Written to {path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=Path("results/ctd/ctd_micronutrient_sets.xlsx"))
    parser.add_argument("--start-references", type=int,
                        default=START_REFERENCES)
    parser.add_argument("--min-size", type=int, default=MIN_SET_SIZE)
    parser.add_argument("--max-size", type=int, default=MAX_SET_SIZE)
    parser.add_argument("--pause", type=float, default=2.0)
    args = parser.parse_args()

    interactions = fetch_all(pause=args.pause)
    direct = interactions[~is_cotreatment(interactions)]

    genes_all = gene_table(interactions)
    genes_direct = gene_table(direct)
    genes_top = top_brain_genes(genes_direct, args.start_references,
                                args.min_size, args.max_size)

    write_workbook(args.output, {
        "Overview_all": overview(genes_all),
        "Overview_no_cotreat": overview(genes_direct),
        "Overview_brain_top": overview(genes_top),
        "Genes_all": genes_all,
        "Genes_no_cotreatment": genes_direct,
        "Genes_brain_top": genes_top,
    })


if __name__ == "__main__":
    main()
