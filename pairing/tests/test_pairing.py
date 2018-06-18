import pytest
import numpy as np
import mdtraj as md

from pairing.utils.io import get_fn
import pairing


def test_generate_direct_correlation():
    trj = md.load(get_fn('sevick1988.gro'))

    ref = np.asarray([[1, 0, 0, 0, 1],
                      [0, 1, 1, 0, 0],
                      [0, 1, 1, 0, 1],
                      [0, 0, 0, 1, 0],
                      [1, 0, 1, 0, 1]], dtype=np.int32)

    gen = pairing.generate_direct_correlation(trj, cutoff=0.8)

    assert (ref == gen).all()


def test_sevick1988():
    """Test the system desribed in the appendix of Sevick 1988,
    doi 10.1063/1.454720"""
    c_D = np.asarray([[1, 0, 0, 0, 1],
                      [0, 1, 1, 0, 0],
                      [0, 1, 1, 0, 1],
                      [0, 0, 0, 1, 0],
                      [1, 0, 1, 0, 1]], dtype=np.int32)

    c_I = np.asarray([[1, 1, 1, 0, 1],
                      [1, 1, 1, 0, 1],
                      [1, 1, 1, 0, 1],
                      [0, 0, 0, 1, 0],
                      [1, 1, 1, 0, 1]], dtype=np.int32)

    assert (c_I == pairing.generate_indirect_connectivity(c_D)).all()


def test_40_atoms():
    trj = md.load(get_fn('sevick1988.gro'))
    direct = pairing.generate_direct_correlation(trj, cutoff=0.8)
    indirect = pairing.generate_indirect_connectivity(direct)

    assert indirect.dtype == np.int32


def test_indirect_matrix_reduction():
    trj = md.load(get_fn('sevick1988.gro'))
    direct = pairing.generate_direct_correlation(trj, cutoff=0.8)
    indirect = pairing.generate_indirect_connectivity(direct)

    c_R = np.asarray([[0, 1],
                      [0, 1],
                      [0, 1],
                      [1, 0],
                      [0, 1]])

    assert (c_R == pairing.generate_clusters(indirect)).all()


def test_cluster_analysis():
    trj = md.load(get_fn('sevick1988.gro'))
    direct = pairing.generate_direct_correlation(trj, cutoff=0.8)
    indirect = pairing.generate_indirect_connectivity(direct)
    reduction = pairing.generate_clusters(indirect)

    assert pairing.analyze_clusters(reduction) == (2.5, 1.5)
