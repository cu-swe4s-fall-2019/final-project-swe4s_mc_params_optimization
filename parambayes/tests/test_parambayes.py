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

    def test_generate_priors(self):
        priors_dict = {"A":["gamma",[1,0,1]], "B":["exponential", [0,1]]}
        mcmc_priors = parambayes.MCMC_Prior(priors_dict)
        mcmc_priors.make_priors()
        assert mcmc_priors.priors["A"]
        assert mcmc_priors.priors["B"]
        assert mcmc_priors.priors["A"]["function"].cdf(1000, *mcmc_priors.priors["A"]["values"]) == 1
        assert mcmc_priors.priors["B"]["function"].cdf(1000, *mcmc_priors.priors["B"]["values"]) == 1
    
    def test_more_distributions(self):
        priors_dict = {"paramC":["uniform",[0,1]], "paramD":["logistic", [0,1]]}
        mcmc_priors = parambayes.MCMC_Prior(priors_dict)   
        mcmc_priors.make_priors()
        assert mcmc_priors.priors["paramC"]["function"].cdf(1000) == 1
        assert mcmc_priors.priors["paramD"]["function"].cdf(1000) == 1
        
    def test_not_implemented_dist(self):
        prior_dict = {"A":["nondist",[1,2,3]]}
        self.assertRaises(KeyError, parambayes.MCMC_Prior, prior_dict)
    
    def test_bad_prior_dict(self):
        prior_dict = [1,2,3,4]
        self.assertRaises(TypeError, parambayes.MCMC_Prior, prior_dict)

        prior_dict = {1:2, 3:4}
        self.assertRaises(TypeError, parambayes.MCMC_Prior, prior_dict)
