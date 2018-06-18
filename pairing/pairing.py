from copy import deepcopy
import itertools

import numpy as np
import mdtraj as md


def generate_direct_correlation(trj, cutoff=1.0):
    """
    Genrate direct correlation matrix from a COM-based mdtraj.Trajectory.

    Parameters
    ----------
    trj : mdtraj.Trajectory
        Trajectory for which "atom" sites are to be considered
    cutoff : float, default = 0.8
        Distance cutoff below which two sites are considered paired

    Returns
    -------
    direct_corr : np.ndarray, dtype=np.int32
        Direct correlation matrix
    """

    size = trj.top.n_residues
    direct_corr = np.zeros((size, size), dtype=np.int32)

    for row in range(size):
        for col in range(size):
            if row == col:
                direct_corr[row, col] = 1
            else:
                dist = md.compute_distances(trj, atom_pairs=[(row, col)])
                if dist < cutoff:
                    direct_corr[row, col] = 1
                    direct_corr[col, row] = 1

    return direct_corr


def generate_indirect_connectivity(direct_corr):
    """
    Parameters
    ----------
    direct_corr: np.ndarray, dtype=np.int32
    direct correlation matrix

    Returns
    -------
    indirect: np.ndarray, dtype=np.int32
    indirect connectivity matrix
    """
    c = deepcopy(direct_corr)
    for row in c:
        ones = np.where(row == 1)[0]
        if len(ones) == 1:
            continue
        else:
            intersect = np.maximum.reduce(c[:, ones].T)
            for ele in ones:
                c[:, ele] = intersect
    indirect = c

    return indirect


def generate_clusters(indirect):
    """
    Generate clusters by reducing the indirect matrix

    Parameters
    ----------
    indirect_corr : numpy.ndarray, dtype=np.int32
        Indirect corrlation matrix

    Returns
    -------
    clusters : numpy.ndarray, dtype=np.int32
        Matrix in which each column represents a cluster and corresponding
        row indices match the indices of sites in a given cluster.
    """
    clusters = np.unique(indirect, axis=1)
    return clusters


def analyze_clusters(clusters):
    """
    Find the average and standard deviation of cluster sizes.

    Parameters
    ----------
    clusters : numpy.ndarray, dtype=np.int32
        Matrix in which each column represents a cluster and corresponding

    Returns
    -------
    avg : float
        Average cluster size
    stdev : float
        Standard deviation of cluster sizes
    """
    cluster_sizes = np.sum(clusters, axis=0)
    avg = np.mean(cluster_sizes)
    stdev = np.std(cluster_sizes)
    return avg, stdev
