from unittest import TestCase
import numpy as np
from numpy.testing import assert_allclose

from fd.analysis import aerodynamic_analysis as ana


class TestUnitconversion(TestCase):
    def test_unitconversion(self):
        assert_allclose(ana.calc_CL(1000, 10), 0.54421769)
        assert_allclose(ana.calc_CL(np.array([1000, 15000]), np.array([10, 30])), np.array([0.54421769, 0.90702948]), rtol=1e-01)
        #assert_allclose(ana.calc_CL([1000, 15000], [10, 30]), np.array[0.54421769, 0.90702948], rtol=1e-01)
        assert_allclose(ana.calc_CL_alpha(np.array([0.1, 0.2, 0.3]), np.array([0, 5, 10])), [0.02, 0.1, -5.], rtol=1e-01)
        assert_allclose(ana.calc_CL_alpha(np.array([0.11, 0.19, 0.31, 0.39]), np.array([0, 5, 10, 15])), [0.02, 0.1, -5.], rtol=1e-01)
        assert_allclose(ana.calc_CD(1000, 10), 0.54421769, rtol=1e-01)
        assert_allclose(ana.calc_CD(np.array([1000, 15000]), np.array([10, 30])), np.array([0.54421769, 0.90702948]), rtol=1e-01)
        assert_allclose(ana.calc_CD0_e(np.array([0.0318, 0.0532, 0.024, 0.063]), np.array([0.5, 0.84, 0.29, 0.955])), [0.02, 0.8], rtol=1e-01)
        assert_allclose(ana.calc_CD0_e(np.array([0.032, 0.053, 0.025, 0.065]), np.array([0.51, 0.83, 0.28, 0.95])), [0.02, 0.8], rtol=1e-01)