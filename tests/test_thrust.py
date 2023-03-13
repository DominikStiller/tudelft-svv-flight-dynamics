from unittest import TestCase

from numpy.testing import assert_allclose

from fd.analysis.thrust import calculate_thrust


class TestThrust(TestCase):
    def test_thrust(self):
        # Calculated from Excel sheet
        assert_allclose(calculate_thrust(3000, 0.4, 269.15, 0.1), 4096.2853587604200)
        assert_allclose(calculate_thrust(3000, 0.4, 269.15, 0.09), 3510.0944255666300)
        assert_allclose(calculate_thrust(5000, 0.4, 256.15, 0.09), 4004.8876358752600)
        assert_allclose(calculate_thrust(5000, 0.8, 256.15, 0.09), 2732.5546243401900)
        assert_allclose(calculate_thrust(5000, 0.1, 256.15, 0.09), 5369.0542444565900)
        assert_allclose(calculate_thrust(100, 0.1, 288.2, 0.1), 4920.5995394974200)
