"""
pairing.py
analyze pairing and clustering of molecular systems

Handles the primary functions
"""

from copy import deepcopy
import itertools

import numpy as np


def generate_indirect_connectivity(direct_corr):
    """
    Genrate indirect correlation matrix from a direct correlation matrix

    Parameters
    ----------
    direct_corr : numpy.ndarray
        Direct correlation matrix from which an indirect correlation matrix
        will be generated.

    Returns
    -------
    indirect_corr : numpy.ndarray
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
                intersect = find_intersection(c[:, combo[0]], c[:, combo[1]])
                c[:, combo[0]] = intersect
                c[:, combo[1]] = intersect

    indirect_corr = c
    return indirect_corr


def find_intersection(a, b):
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

    intersection = np.zeros(len(a))
    for i in range(len(intersection)):
        if a[i] == b[i]:
            intersection[i] = a[i]
        else:
            intersection[i] = np.max([a[i], b[i]])
    return intersection


if __name__ == "__main__":
    # Do something if this file is invoked on its own
    print(canvas())
