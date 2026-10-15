"""Colocalization of two diseases against the micronutrient gene sets.

Each disease runs on the interactome cut to its tissue, so the seeds and the
propagation hold only tissue-expressed genes. Each pair runs in its own
interpreter and is cached. The pair table carries
seed count and median seed degree, which the label permutation does not fix.
"""

from pathlib import Path

import numpy as np
import pandas as pd

import ctd_seeds
import geneset_coloc

DISEASES = ["dravet", "cakut"]
# The keto diet is a dietary panel and not a micronutrient, so it carries its
# own label and stays separable downstream.
MICRONUTRIENTS = list(ctd_seeds.KEYS)

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results" / "micronutrient_coloc"
CACHE = OUT / "pairs"
(OUT / "tables").mkdir(parents=True, exist_ok=True)
(OUT / "figures").mkdir(parents=True, exist_ok=True)

Q_THRESH = 0.05

# The tissue of each disease, and the interactome cut to that tissue.
TISSUE_OF = {"dravet": "brain", "cakut": "kidney"}
NETWORK_OF = {"dravet": "pcnet22_brain", "cakut": "pcnet22_kidney"}


def use_disease_network(disease):
    """Switch ``geneset_coloc`` to the tissue network of ``disease``, if it is not already on it."""
    if geneset_coloc.NETWORK != NETWORK_OF[disease]:
        geneset_coloc.set_network(NETWORK_OF[disease])
    return geneset_coloc.NETWORK


def pair_tag(disease, micronutrient):
    return f"{disease}_{micronutrient}"


def _median_seed_degree(seeds, degrees):
    """The median network degree of a seed list, skipping seeds off the network."""
    values = [degrees[g] for g in seeds if g in degrees]
    return float(np.median(values)) if values else np.nan


def pair_stats_table(results, out_tsv=None):
    """One row per pair, ranked by ``z_mean_z``, with the power proxies.

    ``q_mean_z`` is BH over the pairs of one disease.
    """
    rows = []
    for (disease, micronutrient), res in results.items():
        stats = res["stats"].iloc[0]
        _, threshold = ctd_seeds.genes_for(micronutrient, brain_only=False)
        rows.append({
            "disease": geneset_coloc.LABELS[disease],
            "micronutrient": geneset_coloc.LABELS[micronutrient],
            "disease_key": disease,
            "micronutrient_key": micronutrient,
            "n_disease_seeds": int(stats[f"n_{disease}_seeds"]),
            "n_micronutrient_seeds": int(stats[f"n_{micronutrient}_seeds"]),
            "min_references_used": threshold,
            "median_seed_degree": res["median_seed_degree"],
            "n_seed_overlap": int(stats["n_seed_overlap"]),
            "coloc_size": int(stats["coloc_size"]),
            "mean_combined_z": float(stats["mean_combined_z"]),
            "contrib_both_negative": float(stats["contrib_both_negative"]),
            "z_mean_z": float(stats["z_mean_z"]),
            "p_mean_z": float(stats["p_mean_z"]),
            "z_size": float(stats["z_size"]),
            "p_size": float(stats["p_size"]),
        })
    table = pd.DataFrame(rows)
    table["q_mean_z"] = np.nan
    for disease in table["disease_key"].unique():
        mask = table["disease_key"] == disease
        table.loc[mask, "q_mean_z"] = geneset_coloc.bh(table.loc[mask, "p_mean_z"].values)
    table["significant"] = table["q_mean_z"] < Q_THRESH
    table = table.sort_values(["disease_key", "z_mean_z"], ascending=[True, False])

    if out_tsv:
        Path(out_tsv).parent.mkdir(parents=True, exist_ok=True)
        table.to_csv(out_tsv, sep="\t", index=False)
    return table.reset_index(drop=True)



def coloc_gene_matrix(gene_tables, disease):
    """Genes x micronutrient grid of combined z, for one disease.

    A gene enters when it is in the colocalized module of one pair. An empty cell
    means "not in that module", not "zero".
    """
    columns = {}
    for (disease_key, micronutrient), table in gene_tables.items():
        if disease_key != disease:
            continue
        series = table.set_index("symbol")["combined_z"]
        columns[micronutrient] = series[~series.index.duplicated()]

    if not columns:
        return pd.DataFrame()
    matrix = pd.DataFrame(columns)
    matrix = matrix.loc[[g for g in matrix.index if isinstance(g, str)]]
    matrix.columns = [geneset_coloc.LABELS[c] for c in matrix.columns]
    order = matrix.notna().sum(axis=1) * 1000 + matrix.max(axis=1).fillna(0)
    return matrix.loc[order.sort_values(ascending=False).index]


def coloc_edges(gene_tables, results):
    """Disease-micronutrient-gene edges of the colocalized modules, with role and z."""
    rows = []
    for (disease, micronutrient), table in gene_tables.items():
        res = results[(disease, micronutrient)]
        seeds_d = set(res[f"seeds_{disease}"])
        seeds_m = set(res[f"seeds_{micronutrient}"])
        for _, row in table.iterrows():
            gene = int(row["entrez"])
            role = ("both seeds" if gene in seeds_d and gene in seeds_m else
                    "disease seed" if gene in seeds_d else
                    "micronutrient seed" if gene in seeds_m else "bridge")
            rows.append({
                "disease": geneset_coloc.LABELS[disease],
                "micronutrient": geneset_coloc.LABELS[micronutrient],
                "gene": row["symbol"],
                "entrez": gene,
                "role": role,
                "combined_z": float(row["combined_z"]),
            })
    return pd.DataFrame(rows)


def run_pairs_out_of_process(diseases=None, micronutrients=None, cache=None,
                             python=None):
    """Run every pair in its own interpreter, one at a time, on its disease network.

    Each propagation holds the heats matrix in memory. A pair whose cache file
    exists is skipped, which makes a run restartable.
    """
    import subprocess
    import sys as _sys

    cache = Path(cache or CACHE)
    cache.mkdir(parents=True, exist_ok=True)
    runner = Path(__file__).resolve().parent / "run_one_pair.py"
    failed = []
    for disease in diseases or DISEASES:
        for micronutrient in micronutrients or MICRONUTRIENTS:
            target = cache / f"{disease}_{micronutrient}.pkl"
            if target.exists():
                print(f"cached {target.name}")
                continue
            print(f"running {disease} vs {micronutrient} ...")
            done = subprocess.run([python or _sys.executable, str(runner),
                                   disease, micronutrient, str(cache)])
            if done.returncode:
                failed.append(f"{disease}_{micronutrient}")
                print(f"  FAILED {disease} vs {micronutrient}")
    if failed:
        print(f"\n{len(failed)} pair(s) failed: {', '.join(failed)}")
    return cache


def load_cached(diseases=None, micronutrients=None, cache=None):
    """Read the per-pair caches back as ``(results, gene_tables)``.

    Shaped as ``_colocalize`` returns, so the cached and in-process paths are
    interchangeable downstream.
    """
    import pickle

    cache = Path(cache or CACHE)
    results, gene_tables = {}, {}
    missing = []
    for disease in diseases or DISEASES:
        for micronutrient in micronutrients or MICRONUTRIENTS:
            path = cache / f"{disease}_{micronutrient}.pkl"
            if not path.exists():
                missing.append(path.name)
                continue
            blob = pickle.load(open(path, "rb"))
            key = (disease, micronutrient)
            results[key] = {
                "stats": pd.DataFrame([blob["stats"]]),
                "coloc_genes": blob["coloc_genes"],
                f"seeds_{disease}": blob["seeds_disease"],
                f"seeds_{micronutrient}": blob["seeds_nutrient"],
                "tag": pair_tag(disease, micronutrient),
                "median_seed_degree": blob["median_seed_degree"],
            }
            gene_tables[key] = blob["coloc_genes"]
    if missing:
        print(f"not cached yet: {', '.join(missing)}")
    return results, gene_tables
