"""
Unit and regression test for the pairing package.
"""

# Import package, test suite, and other packages as needed
import pairing
import pytest
import sys

def test_pairing_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "pairing" in sys.modules
