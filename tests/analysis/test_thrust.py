from unittest import TestCase

from numpy.testing import assert_allclose

from fd.analysis.thrust import calculate_thrust, calc_Tc


class TestThrust(TestCase):
    def test_thrust(self):
        # Calculated from Excel sheet
        # Static temperatures in this test are calculated as temperature from ISA + dT
        assert_allclose(calculate_thrust(3000, 0.4, 268.65 + 0.5, 0.1), 4096.2853587604200)
        assert_allclose(calculate_thrust(3000, 0.4, 268.65 + 0.5, 0.09), 3510.0944255666300)
        assert_allclose(calculate_thrust(5000, 0.4, 255.65 + 0.5, 0.09), 4004.8876358752600)
        assert_allclose(calculate_thrust(5000, 0.8, 255.65 + 0.5, 0.09), 2732.5546243401900)
        assert_allclose(calculate_thrust(5000, 0.1, 255.65 + 0.5, 0.09), 5369.0542444565900)
        assert_allclose(calculate_thrust(100, 0.1, 287.50 + 0.7, 0.1), 4920.5995394974200)

    def test_calc_Tc(self):
        assert_allclose(calc_Tc(2000, 300, 1.225, 15), 2.41874527589e-3)
