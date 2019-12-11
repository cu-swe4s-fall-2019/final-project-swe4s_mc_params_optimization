import unittest
import parambayes.utility
from parambayes.LennardJones_2Center_correlations import LennardJones_2C
import numpy as np
class TestUtilityFunctions(unittest.TestCase):
    def test_rhol_hat_models_wrong_types(self):
        lj2c = LennardJones_2C(5)
        self.assertRaises(TypeError, parambayes.utility.rhol_hat_models,1,2,3,4,5,6)
        self.assertRaises(TypeError, parambayes.utility.rhol_hat_models, lj2c, 'nope',2,3,4,5)
        self.assertRaises(TypeError, parambayes.utility.rhol_hat_models, lj2c, 1,2,3,'nope',5)
        self.assertRaises(TypeError, parambayes.utility.rhol_hat_models, lj2c, 1,2,3,4,5)

    def test_rhol_hat_models(self):
        lj2c = LennardJones_2C(5)
        rhol_hat = parambayes.utility.rhol_hat_models(lj2c, np.array([1]), 1,1,1,1)
        assert len(rhol_hat) == 1
        rhol_hat = parambayes.utility.rhol_hat_models(lj2c, np.array([1,0.55, 0.7]), 1,1,1,1)
        assert len(rhol_hat) == 3
