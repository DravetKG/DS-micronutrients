"""PCNet 2.2 in NCBI Gene ID space.

The co-citation-free member of the PCNet 2.x family (NDEx
``8b4b54fa-e87d-11ee-9621-005056ae23aa``), cached as a sparse adjacency.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "input" / "pcnet2.2"
WORK = ROOT / "work" / "pcnet2.2"

def load_nodes(path=None):
    """``DataFrame(entrez, symbol)`` for every PCNet 2.2 node."""
    path = path or DATA / "pcnet22_nodes.json"
    raw = json.load(open(path))
    df = pd.DataFrame([{"cx_id": n["@id"], "symbol": n["n"],
                        "entrez": int(str(n["r"]).split(":")[1])} for n in raw])
    assert (df["cx_id"] == df["entrez"]).all(), \
        "CX node @id is not the NCBI Gene ID; the edge aspect cannot be read as Entrez"
    assert df["entrez"].is_unique, "duplicate NCBI Gene IDs among nodes"
    return df[["entrez", "symbol"]]


def load_edges(path=None):
    """``(src, tgt)`` arrays of raw CX edge endpoints, as Entrez IDs."""
    path = path or DATA / "pcnet22_edges.json"
    raw = json.load(open(path))
    src = np.fromiter((e["s"] for e in raw), dtype=np.int64, count=len(raw))
    tgt = np.fromiter((e["t"] for e in raw), dtype=np.int64, count=len(raw))
    return src, tgt


def standardize(nodes_df, src, tgt, verbose=True):
    """Symmetric, deduped, self-loop-free adjacency over the node table.

    Returns ``(A, nodes)``, ``nodes`` being ascending Entrez IDs. Every removed
    count is printed: an edge-count gap is usually a directed/undirected mix-up.
    """
    nodes = np.sort(nodes_df["entrez"].to_numpy())
    pos = {e: i for i, e in enumerate(nodes)}

    n_raw = len(src)
    known = np.fromiter((s in pos and t in pos for s, t in zip(src, tgt)),
                        dtype=bool, count=n_raw)
    i = np.fromiter((pos[s] for s in src[known]), dtype=np.int64, count=known.sum())
    j = np.fromiter((pos[t] for t in tgt[known]), dtype=np.int64, count=known.sum())

    loops = i == j
    i, j = i[~loops], j[~loops]

    # Orient every pair low->high, then dedupe: this collapses both repeated
    # edges and the same pair given in both directions in one step.
    lo, hi = np.minimum(i, j), np.maximum(i, j)
    key = lo.astype(np.int64) * len(nodes) + hi
    uniq = np.unique(key)
    lo_u, hi_u = uniq // len(nodes), uniq % len(nodes)

    A = sp.coo_matrix(
        (np.ones(2 * len(uniq), dtype=np.int8),
         (np.concatenate([lo_u, hi_u]), np.concatenate([hi_u, lo_u]))),
        shape=(len(nodes), len(nodes)),
    ).tocsr()
    A.data[:] = 1                      # binary; coo sums duplicates otherwise

    if verbose:
        print(f"raw CX edges            {n_raw:,}")
        print(f"  endpoint not a node   -{n_raw - int(known.sum()):,}")
        print(f"  self-loops            -{int(loops.sum()):,}")
        print(f"  duplicate/reciprocal  -{len(i) - len(uniq):,}")
        print(f"undirected edges kept    {len(uniq):,}")
        print(f"nodes                    {len(nodes):,}")

    assert (A != A.T).nnz == 0, "adjacency is not symmetric"
    assert A.diagonal().sum() == 0, "self-loops survived"
    return A, nodes


def load(verbose=True):
    """Standardized PCNet 2.2 as ``(A, nodes, symbol_of)``, cached as npz."""
    cache_A, cache_n = WORK / "pcnet22_adj.npz", WORK / "pcnet22_nodes.npy"
    nodes_df = load_nodes()
    symbol_of = dict(zip(nodes_df["entrez"], nodes_df["symbol"]))

    if cache_A.exists() and cache_n.exists():
        A = sp.load_npz(cache_A)
        nodes = np.load(cache_n)
        if verbose:
            print(f"cache: PCNet 2.2 {A.shape[0]:,} nodes, {A.nnz // 2:,} edges")
        return A, nodes, symbol_of

    src, tgt = load_edges()
    A, nodes = standardize(nodes_df, src, tgt, verbose=verbose)
    WORK.mkdir(parents=True, exist_ok=True)
    sp.save_npz(cache_A, A)
    np.save(cache_n, nodes)
    return A, nodes, symbol_of



def load_subset(symbols, verbose=True):
    """PCNet 2.2 cut to the nodes whose symbol is in ``symbols``, as ``load`` returns it."""
    A, nodes, symbol_of = load(verbose)
    keep = np.array([symbol_of[g] in symbols for g in nodes.tolist()])
    if verbose:
        print(f"tissue cut keeps {keep.sum():,}/{len(nodes):,} nodes")
    nodes = nodes[keep]
    return A[keep][:, keep].tocsr(), nodes, {g: symbol_of[g] for g in nodes.tolist()}
