import unittest
import parambayes.utility
from parambayes.utility import *
from parambayes.LennardJones_2Center_correlations import LennardJones_2C
import numpy as np


class TestUtilityFunctions(unittest.TestCase):
    def test_rhol_hat_models_wrong_types(self):
        lj2c = LennardJones_2C(5)
        self.assertRaises(TypeError,
                          parambayes.utility.rhol_hat_models,
                          1, 2, 3, 4, 5, 6
                          )
        self.assertRaises(TypeError,
                          parambayes.utility.rhol_hat_models,
                          lj2c, 'nope', 2, 3, 4, 5
                          )
        self.assertRaises(TypeError,
                          parambayes.utility.rhol_hat_models,
                          lj2c, 1, 2, 3, 'nope', 5
                          )
        self.assertRaises(TypeError,
                          parambayes.utility.rhol_hat_models,
                          lj2c, 1, 2, 3, 4, 5
                          )

    def test_rhol_hat_models(self):
        lj2c = LennardJones_2C(5)
        rhol_hat = parambayes.utility.rhol_hat_models(lj2c,
                                                      np.array([1]),
                                                      1, 1, 1, 1)
        assert len(rhol_hat) == 1
        rhol_hat = parambayes.utility.rhol_hat_models(lj2c,
                                                      np.array([1, 0.55, 0.7]),
                                                      1, 1, 1, 1)
        assert len(rhol_hat) == 3

    def test_Psat_hat_models_wrong_types(self):
        lj2c = LennardJones_2C(5)
        self.assertRaises(TypeError,
                          parambayes.utility.Psat_hat_models,
                          1, 2, 3, 4, 5, 6)
        self.assertRaises(TypeError,
                          parambayes.utility.Psat_hat_models,
                          lj2c, 'nope', 2, 3, 4, 5)
        self.assertRaises(TypeError,
                          parambayes.utility.Psat_hat_models,
                          lj2c, 1, 2, 3, 'nope', 5)
        self.assertRaises(TypeError,
                          parambayes.utility.Psat_hat_models,
                          lj2c, 1, 2, 3, 4, 5)

    def test_Psat_hat_models(self):
        lj2c = LennardJones_2C(5)
        psat_hat = parambayes.utility.Psat_hat_models(lj2c,
                                                      np.array([1]),
                                                      1, 1, 1, 1)
        assert len(psat_hat) == 1
        psat_hat = parambayes.utility.Psat_hat_models(lj2c,
                                                      np.array([1, 0.55, 0.7]),
                                                      1, 1, 1, 1)
        assert len(psat_hat) == 3

    def test_T_c_hat_models_wrong_types(self):
        lj2c = LennardJones_2C(5)
        self.assertRaises(TypeError, parambayes.utility.T_c_hat_models,
                          1, 2, 3, 4, 5)
        self.assertRaises(TypeError, parambayes.utility.T_c_hat_models,
                          lj2c, 'nope', 2, 3, 4)
        self.assertRaises(TypeError, parambayes.utility.T_c_hat_models,
                          lj2c, 1, 2, 3, 'nope')
        self.assertRaises(TypeError, parambayes.utility.T_c_hat_models,
                          lj2c, np.array([1]), 2, 3, 4)

    def test_T_c_hat_models(self):
        lj2c = LennardJones_2C(5)
        t_c_hat = parambayes.utility.T_c_hat_models(lj2c, 1, 1, 1, 1)
        self.assertAlmostEqual(t_c_hat, 2.61073904)

    def test_compute_percent_deviation(self):
        lj2c = LennardJones_2C(5)
        rhol_mrd, psat_mrd, surftens_mrd, T_c_rd = \
            parambayes.utility.computePercentDeviations(lj2c,
                                                        np.array([1]),
                                                        np.array([1]),
                                                        np.array([1]),
                                                        [1, 1, 1, 1],
                                                        np.array([1]),
                                                        np.array([1]),
                                                        np.array([1]),
                                                        np.array([1]),
                                                        rhol_hat_models,
                                                        Psat_hat_models,
                                                        SurfTens_hat_models,
                                                        T_c_hat_models
                                                        )
        assert rhol_mrd
        assert psat_mrd
        assert surftens_mrd
        assert T_c_rd

    def test_compute_percent_deviation_types(self):
        lj2c = LennardJones_2C(5)
        self.assertRaises(TypeError,
                          parambayes.utility.computePercentDeviations,
                          "nope",
                          "nope",
                          "nope",
                          "nope",
                          "nope",
                          "nope",
                          "nope",
                          "nope",
                          "nope",
                          "nope",
                          "nope",
                          "nope",
                          "nope"
                          )
        self.assertRaises(TypeError,
                          parambayes.utility.computePercentDeviations,
                          lj2c,
                          1,
                          2,
                          3,
                          "nope",
                          "nope",
                          "nope",
                          "nope",
                          "nope",
                          "nope",
                          "nope",
                          "nope",
                          "nope"
                          )
        self.assertRaises(TypeError,
                          parambayes.utility.computePercentDeviations,
                          lj2c,
                          1,
                          2,
                          3,
                          [1, 2, 3, 4],
                          1,
                          2,
                          3,
                          4,
                          "nope",
                          "nope",
                          "nope",
                          "nope"
                          )
