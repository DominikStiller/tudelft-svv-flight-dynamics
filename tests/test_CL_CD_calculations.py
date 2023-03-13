from unittest import TestCase
import numpy as np
from numpy.testing import assert_allclose

from fd.analysis import aerodynamic_analysis as ana


class TestUnitconversion(TestCase):
    def test_unitconversion(self):
        # Calculated from Excel sheet
        assert_allclose(ana.calc_CL(1000, 10), 0.54421769)
        assert_allclose(ana.calc_CL(np.array([1000, 15000]), np.array([10, 30])), np.array([0.54421769, 0.90702948]))
        assert_allclose(ana.CL([1000, 15000], [10, 30]), np.array[0.54421769, 0.90702948])
        assert_allclose(ana.calc_CL_alpha(np.array([0.1, 0.2, 0.3]), np.array([0, 5, 10])), [0.02, -5., 0.1])  # randomly generated
        assert_allclose(conversion.lbshr_to_kgs(125), 0.015749735069444)
        assert_allclose(conversion.lbshr_to_kgs(0), 0)
        assert_allclose(conversion.psi_to_Pa(25), 172368.925)
        assert_allclose(conversion.psi_to_Pa(0), 0)
        assert_allclose(conversion.ftmin_to_ms(0), 0)
        assert_allclose(conversion.ftmin_to_ms(16.1), 0.081788)
        assert_allclose(conversion.lbs_to_kg(2.2), 0.997903214)
        assert_allclose(conversion.lbs_to_kg(0), 0)
        assert_allclose(conversion.kts_to_ms(1.2), 0.617333333333333)
        assert_allclose(conversion.kts_to_ms(0), 0)
        assert_allclose(conversion.ft_to_m(6.5), 1.9812)
        assert_allclose(conversion.ft_to_m(0), 0)
        assert_allclose(conversion.C_to_K(26.2), 299.35)
        assert_allclose(conversion.C_to_K(-67.9), 205.25)
        assert_allclose(conversion.C_to_K(0), 273.15)
