"""
pairing.py
analyze pairing and clustering of molecular systems

Handles the primary functions
"""

from copy import deepcopy
import itertools

import numpy as np
import mdtraj as md


def calc_cluster(trj, cutoff, chunk_size=500):
    """
    Collect cluster information over a trajectory
    """
    indirect_list = []
    for i, frame in enumerate(trj):
        if i % chunk_size == 0:
            frame_trj = trj[i]
            matrix = generate_direct_correlation(frame_trj, cutoff=cutoff)
        else:
            for row in range(len(matrix[0])):
                for col in range(len(matrix[0])):
                    if matrix[row][col] == 1:
                        dist = md.compute_distances(frame,
                                atom_pairs=[(row, col)])
                        if dist < cutoff:
                            continue
                        elif dist > cutoff:
                            matrix[row][col] = 0
                            matrix[col][row] = 0

        indirect = new_generate_indirect(matrix)
        reduction = generate_clusters(indirect)
        indirect_list.append(reduction)

    return indirect_list


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


def _generate_indirect_connectivity(direct_corr):
    """
    Genrate indirect correlation matrix from a direct correlation matrix

    Parameters
    ----------
    direct_corr : numpy.ndarray, dtype=np.int32
        Direct correlation matrix from which an indirect correlation matrix
        will be generated.

    Returns
    -------
    indirect_corr : numpy.ndarray, dtype=np.int32
        Indirect corrlation matrix
    """

    c = deepcopy(direct_corr)
    size = np.shape(direct_corr)
    if size[0] != size[1]:
        raise ValueError('Direct correlation matrix must be square')
    length = size[0]

    for combo in itertools.combinations([_ for _ in range(length)], 2):
        for i in range(length):
            if c[i, combo[0]] == c[i, combo[1]]:
                if c[i, combo[0]] == 0:
                    continue
                intersect = _find_intersection(c[:, combo[0]], c[:, combo[1]])
                c[:, combo[0]] = intersect
                c[:, combo[1]] = intersect

    indirect_corr = c
    return indirect_corr


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


def _find_intersection(a, b):
    """
    Find set intersection of two arrays

    Parameters
    ----------
    a : array-like
        First array to compare
    b : array-like
        Second array to compare

    Returns
    -------
    intersection : array-like
        Set intersection of a and b
    """

    intersection = np.maximum(a, b)
    return intersection


def _check_validity(c_I):
    """
    Check validity of indirect connectivity matrix

    Parameters
    ----------
    c_I : np.ndarray
    indirect connectivity matrix to test

    Returns
    -------
    Boolean 'True' or 'False'
    """

    test_indirect = _generate_indirect_connectivity(c_I)
    return (test_indirect == c_I).all()


def new_generate_indirect(direct_corr):
    """
    Iteratively call '_generate_indirect_connectivity' and
    '_check_validity' to generate valid indirect correlation
    matrices

    Parameters
    ----------
    direct_corr : np.ndarray
    direct correlation matrix

    Returns
    _______
    new_indirect : np.ndarray
    indirect connectivity matrix
    """

    new_indirect = _generate_indirect_connectivity(
            direct_corr)
    while _check_validity(new_indirect) == False:
        new_indirect = _generate_indirect_connectivity(
                new_indirect)

    return new_indirect


if __name__ == "__main__":
    # Do something if this file is invoked on its own
    print(canvas())
