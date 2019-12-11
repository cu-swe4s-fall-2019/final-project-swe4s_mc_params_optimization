import unittest
import sys
import os

test_path = os.path.abspath(os.path.dirname(__file__))
pb_path = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(pb_path)
import LennardJones_2Center_correlations  # noqa: E402


class TestLennardJones2Center(unittest.TestCase):
    def setUp(self):
        M_w_N2 = 28.0134  # Example molecular weight (for N2)
        os.chdir(pb_path)
        self.LJ2C = LennardJones_2Center_correlations.LennardJones_2C(M_w_N2)
        os.chdir(test_path)
        self.assertEqual(self.LJ2C.M_w, M_w_N2)

    def test_T_c_star_hat(self):
        expected = [1.507579e+00, 2.047231e-02, -1.291671e-03,
                    3.319456e-01, 4.136462e-02, 9.755649e-03, -1.715840e-03,
                    -8.173578e-04, 2.301229e-04]
        for i in range(len(expected)):
            self.assertEqual(self.LJ2C.T_c_star_params[i], expected[i])
        expected_val = 1.669
        self.assertAlmostEqual(self.LJ2C.T_c_star_hat(2, 2),
                               expected_val, places=3)

    def test_rho_c_star_hat(self):
        expected = [3.143171e-01, 2.469999e-03, -2.422011e-04, -1.452035e-01,
                    -4.259098e-02, -2.700883e-03, 2.785485e-03, 3.007566e-04,
                    -5.084756e-04]
        for i in range(len(expected)):
            self.assertEqual(self.LJ2C.rho_c_star_params[i], expected[i])
        expected_val = 0.137
        self.assertAlmostEqual(self.LJ2C.rho_c_star_hat(2, 2),
                               expected_val, places=3)

    def test_C1_hat(self):
        expected = [4.333882, 0.1503665, -0.02085311, -1.870607, -0.7103387,
                    -0.5758677, 0.6802547, 0.1358236, -0.1916098]
        for i in range(len(expected)):
            self.assertEqual(self.LJ2C.b_c1[i], expected[i])
        expected_val = 67.978
        self.assertAlmostEqual(self.LJ2C.C1_hat(2, 2, 2),
                               expected_val, places=3)

    def test_C2_hat(self):
        expected = [-2.660590e+01, -1.144385e+00, 1.097780e-01, 1.138729e+02,
                    -1.082939e+02, 1.131664e+01, -1.732358e+01, -1.609370e+00,
                    3.026670e+00]
        for i in range(len(expected)):
            self.assertEqual(self.LJ2C.b_c2[i], expected[i])
        expected_val = 210
        self.assertAlmostEqual(self.LJ2C.C2_hat(2, 2, 2),
                               expected_val, places=3)

    def test_C3_hat(self):
        expected = [-0.1059248, -0.00355973, -0.8836935]
        for i in range(len(expected)):
            self.assertAlmostEqual(self.LJ2C.b_c3[i], expected[i])
        expected_val = 462
        self.assertAlmostEqual(self.LJ2C.C3_hat(2, 2, 2),
                               expected_val, places=3)

    def test_rho_star_hat_2CLJQ(self):
        pass

    def test_rho_hat_2CLJQ(self):
        pass

    def test_rhol_hat_2CLJQ(self):
        pass

    def test_rhov_hat_2CLJQ(self):
        pass

    def test_Psat_star_hat_2CLJQ(self):
        expected_val = 13.375682001104826
        self.assertAlmostEqual(self.LJ2C.Psat_star_hat_2CLJQ(1000, 2, 2),
                               expected_val)

    def test_Psat_hat_2CLJQ(self):
        expected_val = 61.346452923494226
        self.assertAlmostEqual(self.LJ2C.Psat_hat_2CLJQ(1000, 2, 2, 2, 2),
                               expected_val)

    def test_LJ_model(self):
        expected_val = 0.0
        self.assertAlmostEqual(self.LJ2C.LJ_model(2, 2, 2), expected_val)

    def test_Astar_hat(self):
        expected_val = -5.225219108870148
        self.assertAlmostEqual(self.LJ2C.Astar_hat(2, 2), expected_val)

    def test_ST_star_hat_2CLJQ(self):
        expected_val = [-7.15440867e+10]
        self.assertAlmostEqual(self.LJ2C.ST_star_hat_2CLJQ(1000, 20, 30)[0],
                               expected_val[0],
                               delta=25)

    def test_ST_hat_2CLJQ(self):
        expected_val = 1.8114945309985577e-07
        self.assertEqual(expected_val,
                         self.LJ2C.ST_hat_2CLJQ(1000, 400, 400, 20, 300))

    def test_T_c_hat_2CLJQ(self):
        expected_val = 738.7807737911776
        self.assertEqual(self.LJ2C.T_c_hat_2CLJQ(400, 20, 20, 300),
                         expected_val)


if __name__ == '__main__':
    unittest.main()

