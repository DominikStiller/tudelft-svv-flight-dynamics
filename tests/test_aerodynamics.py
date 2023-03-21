from unittest import TestCase

from numpy.testing import assert_allclose

from fd.analysis.aerodynamics import *


class TestAerodynamics(TestCase):
    def test_calc_true_V(self):
        assert_allclose(calc_true_V(600, 0.9), 441.937574777)
        assert_allclose(calc_true_V(200, 0.7), 198.452160482)
        assert_allclose(calc_true_V(555, 0.21), 99.1764547914)

    def test_calc_equivalent_V(self):
        assert_allclose(calc_equivalent_V(100, 1.225), 100)
        assert_allclose(calc_equivalent_V(250, 0.82), 204.540300904)
        assert_allclose(calc_equivalent_V(75, 0.105), 21.9577516413)

    def test_calc_Tc(self):
        assert_allclose(calc_Tc(2000, 300, 1.225, 15), 2.41874527589e-3)


    def test_calc_CL(self):
        assert_allclose(calc_CL(1000, 10, 1.225), 0.54421769)
        assert_allclose(
            calc_CL(np.array([1000, 15000]), np.array([10, 30]), 1.225),
            np.array([0.54421769, 0.90702948]),
            rtol=1e-01,
        )
        # assert_allclose(calc_CL([1000, 15000], [10, 30]), np.array[0.54421769, 0.90702948], rtol=1e-01)

    def test_estimate_CL_alpha(self):
        assert_allclose(
            estimate_CL_alpha(np.array([0.1, 0.2, 0.3]), np.array([0, 5, 10])),
            [0.02, 0.1, -5.0],
            rtol=1e-01,
        )
        assert_allclose(
            estimate_CL_alpha(np.array([0.11, 0.19, 0.31, 0.39]), np.array([0, 5, 10, 15])),
            [0.02, 0.1, -5.0],
            rtol=1e-01,
        )

    def test_calc_CD(self):
        assert_allclose(calc_CD(1000, 10, 1.225), 0.54421769, rtol=1e-01)
        assert_allclose(
            calc_CD(np.array([1000, 15000]), np.array([10, 30]), 1.225),
            np.array([0.54421769, 0.90702948]),
            rtol=1e-01,
        )

    def test_calc_CD0_e(self):
        assert_allclose(
            calc_CD0_e(
                np.array([0.0318, 0.0532, 0.024, 0.063]), np.array([0.5, 0.84, 0.29, 0.955])
            ),
            [0.02, 0.8],
            rtol=1e-01,
        )
        assert_allclose(
            calc_CD0_e(
                np.array([0.032, 0.053, 0.025, 0.065]), np.array([0.51, 0.83, 0.28, 0.95])
            ),
            [0.02, 0.8],
            rtol=1e-01,
        )

    def test_calc_Cmdelta(self):
        assert_allclose(calc_Cmdelta(20, 19, 2, 1, 10000, 10000, 120, 0.6), -0.037513002)
        assert_allclose(calc_Cmdelta(20, 19, 2, 1, 8000, 12000, 120, 0.6), -0.037513002)
        assert_allclose(
            calc_Cmdelta(20.01, 19.99, 1.6, 1, 8000, 12000, 110, 0.2), -0.00446435726
        )

    def test_estimate_Cmalpha(self):
        assert_allclose(estimate_Cmalpha([1, 2, 3], [0.5, 1, 1.5], -0.01), 0.005)
        assert_allclose(
            estimate_Cmalpha([1.01, 2, 3], [0.5, 1.01, 1.5], -0.01), 0.005, rtol=1e-1
        )
