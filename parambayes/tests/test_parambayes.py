"""
Unit and regression test for the parambayes package.
"""

# Import package, test suite, and other packages as needed
import parambayes
import pytest
import sys

def test_parambayes_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "parambayes" in sys.modules
