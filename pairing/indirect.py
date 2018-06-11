"""
indirect.py
Calculates indirect correlation matrix of molecular systems

"""

from copy import deepcopy
import itertools
import pairing.pairing

import numpy as np

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
                intersect = pairing.pairing._find_intersection(c[:, combo[0]], c[:, combo[1]])
                c[:, combo[0]] = intersect
                c[:, combo[1]] = intersect

    indirect_corr = c
    return indirect_corr
