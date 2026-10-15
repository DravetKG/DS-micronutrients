"""Memory-lean drop-in for ``netcoloc.netprop`` on dense interactomes.

Computes ``H = (1 - alpha) (I - alpha A D^-1)^-1`` from a sparse adjacency in
float32.
"""

import numpy as np
import scipy.linalg
import scipy.sparse as sp


def normalized_adjacency(A):
    """``A D^-1`` as CSR - column j divided by the degree of node j.

    This is the transpose of NetColoc's ``w_prime``, which it transposes before
    inverting, so building it in this orientation skips a copy.
    """
    deg = np.asarray(A.sum(axis=1)).ravel()
    if (deg == 0).any():
        raise ValueError("graph has degree-zero nodes; restrict to the "
                         "largest connected component first")
    return A.multiply(sp.csr_matrix(1.0 / deg)).tocsr()


def individual_heats_matrix(A, alpha=0.5, dtype=np.float32, verbose=True):
    """Pre-computed heat contribution of every gene to every other gene.

    Column ``j`` is the propagated heat from seeding gene ``j``.
    """
    assert 1 >= alpha >= 0, "Alpha must be between 0 and 1"
    n = A.shape[0]
    if verbose:
        print(f"densifying {n}x{n} to {dtype().itemsize * n * n / 2**30:.2f} GB "
              f"({np.dtype(dtype).name})")

    M = normalized_adjacency(A).toarray().astype(dtype, copy=False)
    M *= -alpha
    M[np.diag_indices(n)] += 1.0          # I - alpha * A D^-1, built in place

    if verbose:
        print("inverting...")
    H = scipy.linalg.inv(M, overwrite_a=True)
    H *= (1.0 - alpha)
    return H


def propagate(H, node_index, seed_genes):
    """Mean heat over the seed columns of ``H``.

    Takes a ``{symbol: row}`` dict instead of a linear scan per seed.
    """
    cols = [node_index[g] for g in seed_genes if g in node_index]
    if not cols:
        raise ValueError("no seed genes found in the interactome")
    return H[:, cols].mean(axis=1)


def test_matches_netcoloc(n=60, seed=0):
    """Check this module against ``netcoloc.netprop`` on a small graph."""
    import networkx as nx
    from netcoloc.netprop import (get_individual_heats_matrix,
                                  network_propagation)

    rng = np.random.default_rng(seed)
    G = nx.gnp_random_graph(n, 0.15, seed=seed)
    G.remove_nodes_from([v for v, d in G.degree() if d == 0])
    G = G.subgraph(max(nx.connected_components(G), key=len)).copy()
    nodes = list(G.nodes())
    G = nx.relabel_nodes(G, {v: f"G{v}" for v in nodes})
    names = [f"G{v}" for v in nodes]

    ref = get_individual_heats_matrix(G, alpha=0.5, conserve_heat=True)
    A = nx.to_scipy_sparse_array(G, nodelist=names, format="csr", dtype=np.float64)
    mine = individual_heats_matrix(A, alpha=0.5, dtype=np.float64, verbose=False)

    assert np.allclose(ref, mine, atol=1e-10), \
        f"heats matrix differs, max |d| = {np.abs(ref - mine).max():.3e}"

    seeds = list(rng.choice(names, 5, replace=False))
    ref_p = network_propagation(ref, names, seeds).loc[names].to_numpy()
    mine_p = propagate(mine, {g: i for i, g in enumerate(names)}, sorted(seeds))
    assert np.allclose(ref_p, mine_p, atol=1e-12), "propagation differs"
    print(f"OK: matches netcoloc on {len(names)} nodes "
          f"(max |d| = {np.abs(ref - mine).max():.2e})")


if __name__ == "__main__":
    test_matches_netcoloc()
