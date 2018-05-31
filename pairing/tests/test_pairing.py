"""
Unit and regression test for the pairing package.
"""

# Import package, test suite, and other packages as needed
import pytest
import numpy as np

import pairing


def test_sevick1988():
    """Test the system desribed in the appendix of Sevick 1988,
    doi 10.1063/1.454720"""
    c_D = np.asarray([[1, 0, 0, 0, 1],
                      [0, 1, 1, 0, 0],
                      [0, 1, 1, 0, 1],
                      [0, 0, 0, 1, 0],
                      [1, 0, 1, 0, 1]])

    c_I = np.asarray([[1, 1, 1, 0, 1],
                      [1, 1, 1, 0, 1],
                      [1, 1, 1, 0, 1],
                      [0, 0, 0, 1, 0],
                      [1, 1, 1, 0, 1]])

    assert (c_I == pairing.generate_indirect_connectivity(c_D)).all()
