from unittest import TestCase
import numpy as np
from numpy.testing import assert_allclose

from fd.analysis import aerodynamic_analysis as ana


class TestUnitconversion(TestCase):
    def test_unitconversion(self):
        # Calculated from Excel sheet
        assert_allclose(ana.calc_CL(1000, 10), 0.54421769)
        assert_allclose(ana.calc_CL(np.array([1000, 15000]), np.array([10, 30])), np.array([0.54421769, 0.90702948]))
        #assert_allclose(ana.calc_CL([1000, 15000], [10, 30]), np.array[0.54421769, 0.90702948])
        assert_allclose(ana.calc_CL_alpha(np.array([0.1, 0.2, 0.3]), np.array([0, 5, 10])), [0.02, 0.1, -5.])  # randomly generated
