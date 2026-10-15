"""Colocalization of two gene sets on a chosen interactome, by netcoloc.
"""

import random as rn
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "input"
RESULTS = ROOT / "results" / "coloc"
TABLES = RESULTS / "tables"
sys.path.insert(0, str(ROOT / "scripts"))

from netcoloc.network_colocalization import (
    calculate_expected_overlap,
    calculate_mean_z_score_distribution,
    get_p_from_permutation_results,
)

import ctd_seeds
import kg_ppi
import kidney
import pcnet22
import propagate

# Interactome the propagation runs on: PCNet 2.2 cut to the HPA genes with
# nTPM > 1 in the disease tissue.
NETWORK = "pcnet22_brain"
NETWORKS = {
    "pcnet22_brain": lambda verbose=True: pcnet22.load_subset(
        kg_ppi.brain_genes(verbose=False), verbose),
    "pcnet22_kidney": lambda verbose=True: pcnet22.load_subset(
        kidney.kidney_genes(), verbose),
}

NUM_REPS = 1000     # degree-matched propagations behind each NPS vector
Z_COLOC = 3.0
Z_INDIVIDUAL = 1.5

# Label-permutation settings.
NETCOLOC_NUM_REPS = 1000
OVERLAP_CONTROL = "remove"

def nps_dir():
    return ROOT / "work" / "nps_sets" / NETWORK


# Panels given by symbol. Resolved to Entrez against the network node table, so
# a symbol the network does not carry fails loudly.
DRAVET_SYMBOLS = ("GABRA1", "GABRG2", "PCDH19", "SCN1A", "SCN2A", "SCN9A", "SCN1B")

# CAKUT genes: the ``Protein`` column of the first sheet of ``cakut.xlsx``.
CAKUT_XLSX = DATA / "cakut" / "cakut.xlsx"
CAKUT_SHEET = "Sheet1"
CAKUT_COLUMN = "Protein"


def cakut_symbols():
    """The CAKUT gene symbols of the ``cakut.xlsx`` protein column."""
    table = pd.read_excel(CAKUT_XLSX, sheet_name=CAKUT_SHEET)
    symbols = table[CAKUT_COLUMN].dropna().astype(str).str.strip()
    return tuple(sorted(set(symbols[symbols != ""])))


SYMBOL_SETS = {"dravet": ("Dravet", DRAVET_SYMBOLS)}

FILE_SYMBOL_SETS = {
    "cakut": ("CAKUT", cakut_symbols),
    **{k: (ctd_seeds.LABELS[k],
           (lambda k: lambda: ctd_seeds.genes_for(k, brain_only=False)[0])(k))
       for k in ctd_seeds.KEYS},
}

LABELS = {**{k: v[0] for k, v in SYMBOL_SETS.items()},
          **{k: v[0] for k, v in FILE_SYMBOL_SETS.items()}}


def set_network(name):
    """Switch the interactome. Clears the in-process network and seed caches."""
    global NETWORK
    if name not in NETWORKS:
        raise KeyError(f"unknown network {name!r}; have {sorted(NETWORKS)}")
    NETWORK = name
    _CONTEXT.clear()
    _SEEDS_CACHE.clear()
    return NETWORK


def set_pair(tag):
    """Point the output directories at one pair, ``results/{tag}``."""
    global RESULTS, TABLES
    RESULTS = ROOT / "results" / tag
    TABLES = RESULTS / "tables"
    return RESULTS


def annotation_path(tag, sources=None):
    """Annotation CSV of one comparison. ``sources`` names a non-default source
    list, which gets its own file so the default result is never overwritten.
    Tagged by network, like every other output of a comparison."""
    if sources is None:
        return TABLES / f"coloc_{tag}_annotation_{NETWORK}.csv"
    key = "_".join(s.replace(":", "") for s in sources)
    return TABLES / f"coloc_{tag}_annotation_{key}_{NETWORK}.csv"


_CONTEXT = {}


def network():
    """``(A, nodes, symbol_of, degrees)`` for the ``NETWORK`` largest component."""
    if not _CONTEXT:
        A_full, nodes_all, symbol_of = NETWORKS[NETWORK](verbose=False)
        A, nodes = propagate.largest_component(A_full, nodes_all, verbose=False)
        degrees = dict(zip(nodes.tolist(),
                           np.asarray(A.sum(axis=1)).ravel().tolist()))
        _CONTEXT.update(A=A, nodes=nodes, symbol_of=symbol_of, degrees=degrees)
    c = _CONTEXT
    return c["A"], c["nodes"], c["symbol_of"], c["degrees"]


def context():
    """``(A, nodes, symbol_of, degrees, H)`` for the ``NETWORK`` largest component.

    Cached in-process: ``H`` alone is 1.4 GB.
    """
    A, nodes, symbol_of, degrees = network()
    if "H" not in _CONTEXT:
        _CONTEXT["H"] = propagate.heats_matrix(A, network=NETWORK, verbose=False)
    return A, nodes, symbol_of, degrees, _CONTEXT["H"]


def nps_for(key, seeds, nodes, degrees, H, verbose=True):
    """NPS vector for one seed set, cached per seed set in ``nps_dir()``.

    Keyed by a hash of the seed list, so editing an input sheet invalidates it.
    """
    import hashlib

    digest = hashlib.sha1(
        ",".join(map(str, sorted(seeds))).encode()).hexdigest()[:8]
    d = nps_dir()
    d.mkdir(parents=True, exist_ok=True)
    path = d / f"{key}_{digest}.tsv"

    if path.exists():
        s = pd.read_csv(path, sep="\t", index_col=0).iloc[:, 0]
        s.index = s.index.astype(np.int64)
        return s.loc[nodes].to_numpy()

    t = time.time()
    z, _, _ = propagate.degree_matched_reference(
        H, nodes.tolist(), degrees, seeds, num_reps=NUM_REPS, alpha=propagate.ALPHA,
        minimum_bin_size=propagate.MIN_BIN_SIZE, random_seed=propagate.RANDOM_SEED,
    )
    pd.Series(z, index=nodes, name=key).to_csv(path, sep="\t")
    if verbose:
        print(f"  {key:22s} {len(seeds):4d} seeds  {time.time() - t:5.0f}s  "
              f"max NPS {z.max():6.1f}")
    return z

def symbol_seeds(key, symbol_of):
    """One symbol panel as Entrez IDs, resolved on the ``NETWORK`` nodes.

    The symbols that the network does not carry are dropped and named.
    """
    entrez_of = {s: g for g, s in symbol_of.items()}
    symbols = SYMBOL_SETS[key][1]
    missing = [s for s in symbols if s not in entrez_of]
    if missing:
        print(f"{LABELS[key]}: {len(symbols) - len(missing)}/{len(symbols)} "
              f"symbols on the {NETWORK} node table, dropped {missing}")
    return [entrez_of[s] for s in symbols if s in entrez_of]


def file_symbol_seeds(key, symbol_of):
    """One file-backed symbol list as Entrez IDs; missing symbols are dropped."""
    entrez_of = {s: g for g, s in symbol_of.items()}
    symbols = FILE_SYMBOL_SETS[key][1]()
    seeds = [entrez_of[s] for s in symbols if s in entrez_of]
    missing = len(symbols) - len(seeds)
    print(f"{LABELS[key]}: {len(seeds)}/{len(symbols)} symbols on the "
          f"{NETWORK} node table" + (f", {missing} dropped" if missing else ""))
    return seeds


_SEEDS_CACHE = {}

# The gene sets this project colocalizes. Keys double as the NPS cache key and
# the output column suffix.
SEED_SETS = {
    **{k: (lambda k: lambda s: symbol_seeds(k, s))(k) for k in SYMBOL_SETS},
    **{k: (lambda k: lambda s: file_symbol_seeds(k, s))(k)
       for k in FILE_SYMBOL_SETS},
}
NPS_KEY = {}


def seeds_of(name, symbol_of=None):
    """Seed list for a named gene set, resolved once per process."""
    if name not in _SEEDS_CACHE:
        if symbol_of is None:
            _, _, symbol_of, _ = network()
        _SEEDS_CACHE[name] = SEED_SETS[name](symbol_of)
    return _SEEDS_CACHE[name]


def _as_z(vec, nodes):
    """``vec`` as a one-column frame named ``z``, indexed by gene.

    Both package helpers select a column by that name, so a differently named
    series comes back empty.
    """
    return pd.DataFrame({"z": np.asarray(vec, dtype=float)}, index=nodes)


def _netcoloc_overlap(z_a, z_b, seeds_a, seeds_b, nodes,
                      num_reps=NETCOLOC_NUM_REPS, seed=propagate.RANDOM_SEED):
    """``calculate_expected_overlap`` at the module defaults.

    The package seeds no RNG, so both generators it can reach are seeded here.
    """
    np.random.seed(seed)
    rn.seed(seed)
    return calculate_expected_overlap(
        _as_z(z_a, nodes), _as_z(z_b, nodes),
        seed1=list(seeds_a), seed2=list(seeds_b),
        z_score_threshold=Z_COLOC,
        z1_threshold=Z_INDIVIDUAL, z2_threshold=Z_INDIVIDUAL,
        overlap_control=OVERLAP_CONTROL, num_reps=num_reps, plot=False)


def _netcoloc_mean_z(z_a, z_b, seeds_a, seeds_b, nodes,
                     num_reps=NETCOLOC_NUM_REPS, seed=propagate.RANDOM_SEED):
    """``calculate_mean_z_score_distribution`` at the same defaults.

    ``zero_double_negatives=False``: a gene both propagations avoid contributes a
    positive product, in the null as well as the observed value.
    """
    np.random.seed(seed)
    rn.seed(seed)
    return calculate_mean_z_score_distribution(
        _as_z(z_a, nodes), _as_z(z_b, nodes), num_reps=num_reps,
        zero_double_negatives=False, overlap_control=OVERLAP_CONTROL,
        seed1=list(seeds_a), seed2=list(seeds_b))


def sign_quadrants(z1, z2):
    """How each sign quadrant contributes to ``mean(z1 * z2)``.

    The raw mean is a sum of large cancelling parts, so without this a positive
    mean_combined_z cannot be read as shared proximity.
    """
    prod = z1 * z2
    n = len(prod)
    both_neg, both_pos = (z1 < 0) & (z2 < 0), (z1 > 0) & (z2 > 0)
    mixed = ~both_neg & ~both_pos
    return {"mean_combined_z": float(prod.mean()),
            "pct_both_negative": round(100 * both_neg.mean(), 1),
            "contrib_both_negative": float(prod[both_neg].sum() / n),
            "contrib_both_positive": float(prod[both_pos].sum() / n),
            "contrib_mixed_sign": float(prod[mixed].sum() / n)}


def bh(p):
    """Benjamini-Hochberg adjusted p-values, over the genes of one comparison."""
    p = np.asarray(p, dtype=float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    q = np.empty(n)
    running = 1.0
    for rank, i in zip(range(n, 0, -1), order[::-1]):
        running = min(running, p[i] * n / rank)
        q[i] = running
    return q


def _colocalize(name_a, key_a, seeds_a_raw, name_b, key_b, seeds_b_raw,
                tag, extra_stats=None, verbose=True):
    """Shared pipeline: propagate both seed sets, colocalize, permute.

    ``name_a``/``name_b`` name the output columns, ``key_a``/``key_b`` are the NPS
    cache keys, and ``tag`` names the result files.
    """
    A, nodes, symbol_of, degrees, H = context()
    pos = {int(g): i for i, g in enumerate(nodes.tolist())}

    seeds_a = [g for g in seeds_a_raw if g in pos]
    seeds_b = [g for g in seeds_b_raw if g in pos]
    if verbose:
        print(f"\non the {NETWORK} largest component: {name_a} {len(seeds_a)}/"
              f"{len(seeds_a_raw)}, {name_b} {len(seeds_b)}/{len(seeds_b_raw)}")

    both = sorted(set(seeds_a) & set(seeds_b))
    if both and verbose:
        print(f"seed overlap between the two sets: {len(both)} genes "
              f"({', '.join(symbol_of.get(g, str(g)) for g in both)}) - kept in "
              "the reported gene list and flagged by is_dual_seed, removed "
              "from the permutation by overlap_control='remove'")

    z_a = nps_for(key_a, seeds_a, nodes, degrees, H, verbose=verbose)
    z_b = nps_for(key_b, seeds_b, nodes, degrees, H, verbose=verbose)
    NPS = pd.DataFrame({name_a: z_a, name_b: z_b}, index=nodes)

    # The reported gene list covers every node, the doubly-seeded genes
    # included; ``is_dual_seed`` marks them. ``overlap_control="remove"`` still
    # drops them inside the permutation, so ``coloc_size`` counts fewer genes
    # than the table holds. The quadrant decomposition stays on the tested
    # universe, to match the permuted statistics.
    tested = nodes if not both else nodes[~np.isin(nodes, both)]
    keep = np.isin(nodes, tested)
    za_k, zb_k = z_a[keep], z_b[keep]

    prod = z_a * z_b
    hit = (prod >= Z_COLOC) & (z_a > Z_INDIVIDUAL) & (z_b > Z_INDIVIDUAL)
    overlap_genes = nodes[hit].tolist()
    if verbose:
        print(f"\nsignificance: netcoloc label permutation, "
              f"{NETCOLOC_NUM_REPS:,} reps, overlap_control={OVERLAP_CONTROL!r}, "
              f"z1=z2={Z_INDIVIDUAL}, product>={Z_COLOC} ...")

    obs_size, null_size = _netcoloc_overlap(z_a, z_b, seeds_a, seeds_b, nodes)
    obs_mean, null_mean = _netcoloc_mean_z(z_a, z_b, seeds_a, seeds_b, nodes)

    row = {**(extra_stats or {}), f"n_{name_a}_seeds": len(seeds_a),
           f"n_{name_b}_seeds": len(seeds_b), "n_seed_overlap": len(both),
           "coloc_size": float(obs_size), "n_coloc_reported": len(overlap_genes),
           "n_coloc_dual_seed": int(sum(g in set(both) for g in overlap_genes)),
           "mean_combined_z": float(obs_mean),
           **sign_quadrants(za_k, zb_k)}
    for stat, obs, null in [("size", obs_size, null_size),
                            ("mean_z", obs_mean, null_mean)]:
        null = np.asarray(null, dtype=float)
        row[f"null_{stat}_mean"] = float(null.mean())
        row[f"null_{stat}_sd"] = float(null.std())
        row[f"z_{stat}"] = float((obs - null.mean()) / null.std())
        row[f"p_{stat}"] = float(get_p_from_permutation_results(obs, null))
    stats = pd.DataFrame([row])

    TABLES.mkdir(parents=True, exist_ok=True)
    stats_path = TABLES / f"coloc_{tag}_stats_{NETWORK}.tsv"
    coloc_path = TABLES / f"coloc_{tag}_{NETWORK}.tsv"
    stats.to_csv(stats_path, sep="\t", index=False)
    if verbose:
        print("\n" + stats.T.to_string(header=False))

    seed_a_set, seed_b_set, both_set = set(seeds_a), set(seeds_b), set(both)
    coloc = pd.DataFrame({
        "entrez": overlap_genes,
        "symbol": [symbol_of.get(g) for g in overlap_genes],
        f"NPS_{name_a}": z_a[hit],
        f"NPS_{name_b}": z_b[hit],
        "combined_z": prod[hit],
        f"is_{name_a}_seed": [g in seed_a_set for g in overlap_genes],
        f"is_{name_b}_seed": [g in seed_b_set for g in overlap_genes],
        "is_dual_seed": [g in both_set for g in overlap_genes],
    }).sort_values([f"is_{name_a}_seed", f"is_{name_b}_seed", f"NPS_{name_a}"],
                   ascending=False)
    coloc.to_csv(coloc_path, sep="\t", index=False)
    if verbose:
        print(f"\nwrote {coloc_path.name} ({len(coloc)} genes), {stats_path.name}")

    return {"NPS": NPS, "coloc_genes": coloc, "stats": stats, "tag": tag,
            f"seeds_{name_a}": seeds_a, f"seeds_{name_b}": seeds_b,
            "dual_seed_genes": both,
            "perm_size": null_size, "perm_mean": null_mean,
            "obs_size": obs_size, "obs_mean": obs_mean}


def seed_qc_table(names=None, symbol_of=None):
    """Gene count, network coverage and seed degree for each registered seed set.

    Run before the pair matrix: seed set size changes behaviour under the null.
    """
    A, nodes, sym, degrees = network()
    if symbol_of is None:
        symbol_of = sym
    on_net = {int(g) for g in nodes.tolist()}
    names = list(names or SYMBOL_SETS)

    rows = []
    for key in names:
        seeds = seeds_of(key, symbol_of)
        kept = [g for g in seeds if g in on_net]
        rows.append({"geneset": key, "trait": LABELS.get(key, key),
                     "n_seeds": len(seeds), "n_on_lcc": len(kept),
                     "median_seed_degree":
                         float(np.median([degrees[g] for g in kept]))
                         if kept else np.nan})
    return pd.DataFrame(rows).set_index("geneset")


def main_geneset_pair(name_a, name_b, verbose=True):
    """Colocalize two registered gene sets."""
    _, _, symbol_of, _ = network()
    a, b = name_a, name_b
    return _colocalize(a, NPS_KEY.get(a, a), seeds_of(a, symbol_of),
                       b, NPS_KEY.get(b, b), seeds_of(b, symbol_of),
                       tag=f"{a}_{b}", verbose=verbose)


def nps_background(source, symbol_of):
    """Symbols for every gene the NPS scores were computed over.

    The g:Profiler background: the ``NETWORK`` largest component, not the whole
    genome. Cached per-pair results carry no ``NPS``, so pass ``nodes`` there.
    """
    ids = source["NPS"].index if isinstance(source, dict) else source
    return [s for s in (symbol_of.get(int(g)) for g in ids) if s]


def annotate_coloc(coloc, background, out_csv=None, expressed=None,
                   sources=("GO:BP", "GO:CC", "REAC"), reuse=False):
    """g:Profiler over ``sources`` on the colocalized gene set.

    ``background`` is identical across comparisons. ``expressed`` is a set of
    tissue-expressed symbols. It restricts query and background together, or
    the tissue expression itself comes out enriched.
    """
    from gprofiler import GProfiler

    if out_csv is None:
        out_csv = TABLES / "coloc_annotation.csv"
    if reuse and Path(out_csv).exists():
        res = pd.read_csv(out_csv)
        print(f"reusing {Path(out_csv).name} ({len(res)} terms)")
        return res

    genes = coloc["symbol"].dropna().tolist()
    if expressed is not None:
        n_before = len(genes)
        genes = [g for g in genes if g in expressed]
        background = [b for b in background if b in expressed]
        print(f"tissue filter: {n_before} -> {len(genes)} query genes, "
              f"background -> {len(background):,}")

    if not genes:
        print("  no tissue-expressed genes in this set")
        return None

    print(f"\nannotating {len(genes)} colocalized genes "
          f"against a background of {len(background):,}")
    res = GProfiler(return_dataframe=True).profile(
        organism="hsapiens", query=genes, sources=list(sources),
        background=list(background), domain_scope="custom",
        significance_threshold_method="fdr", user_threshold=0.05, no_evidences=True,
    )
    if not len(res):
        print("  no BH-significant term")
        return None
    TABLES.mkdir(parents=True, exist_ok=True)
    res.to_csv(out_csv, index=False)
    print(f"wrote {Path(out_csv).name} ({len(res)} terms)")
    return res
