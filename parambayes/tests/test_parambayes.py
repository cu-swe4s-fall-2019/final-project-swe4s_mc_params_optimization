"""
Unit and regression test for the parambayes package.
"""

# Import package, test suite, and other packages as needed
import parambayes
import unittest
import sys

def test_parambayes_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "parambayes" in sys.modules

class TestMCMCPrior(unittest.TestCase):
    def test_object_init(self):
        priors_dict = {"A":["gamma",[1,0,1]], "B":["exponential", [0,1]]}
        mcmc_priors = parambayes.MCMC_Prior(priors_dict)
        assert mcmc_priors.dnorm
        assert mcmc_priors.dgamma
        assert mcmc_priors.dgengamma
        assert mcmc_priors.duni
        assert mcmc_priors.dexp