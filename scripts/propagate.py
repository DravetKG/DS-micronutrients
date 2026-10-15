"""Network propagation on an interactome and the degree-matched NPS reference.

The heat is ``F = (1 - alpha) (I - alpha W)^-1 Y0`` at alpha = 0.5. NPS is the
per-gene z of log heat against 1,000 degree-matched random seed sets.
"""

import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

import netprop_big

ALPHA = 0.5
MIN_BIN_SIZE = 10
RANDOM_SEED = 1


def heats_path(network):
    """Heats cache of one network. One file per network, never shared."""
    return ROOT / "work" / "pcnet2.2" / f"{network}_heats_a05.npy"


def largest_component(A, nodes, verbose=True):
    """Restrict to the largest connected component.

    A gene outside the seed's component receives no heat, so its z is undefined.
    """
    n_comp, labels = sp.csgraph.connected_components(A, directed=False)
    keep = labels == np.bincount(labels).argmax()
    if verbose:
        print(f"components {n_comp}; largest keeps {keep.sum():,}/{len(nodes):,} "
              f"nodes ({(~keep).sum():,} dropped)")
    return A[keep][:, keep].tocsr(), nodes[keep]


def heats_matrix(A, network, verbose=True):
    """Cached ``H`` for the network ``A``, float32 (1.1 GB on the brain network).

    Built from the sparse adjacency by ``netprop_big``, which fits in memory where
    the stock float64 path does not. A cache of the wrong shape is rebuilt.
    """
    n = A.shape[0]
    heats = heats_path(network)
    if heats.exists():
        H = np.load(heats, mmap_mode=None)
        if H.shape == (n, n):
            if verbose:
                print(f"cache: {heats.name}  {H.shape}  {H.nbytes / 2**30:.2f} GB")
            return H
        print(f"cache {heats.name} is {H.shape}, network is {n}x{n} - rebuilding")
        del H

    t = time.time()
    H = netprop_big.individual_heats_matrix(A, alpha=ALPHA, verbose=verbose)
    np.save(heats, H)
    if verbose:
        print(f"built H in {time.time() - t:.0f}s; mean column sum "
              f"{H.sum(axis=0).mean():.6f} (should be ~1.0)")
    return H



def degree_matched_reference(H, nodes, degrees, seeds, num_reps=1000, alpha=0.5,
                             minimum_bin_size=10, random_seed=1):
    """
    """
    from netcoloc.netprop_zscore import calculate_heat_zscores

    z, _, random_heats = calculate_heat_zscores(
        H, nodes, degrees, seeds, num_reps=num_reps, alpha=alpha,
        minimum_bin_size=minimum_bin_size, random_seed=random_seed,
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        log_random = np.log(random_heats)
    mu = np.nanmean(log_random, axis=0)
    sd = np.nanstd(log_random, axis=0)
    return z.loc[nodes].to_numpy(), mu, sd

