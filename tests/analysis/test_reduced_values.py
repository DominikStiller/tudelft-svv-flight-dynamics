from unittest import TestCase

from numpy.testing import assert_allclose

from fd.analysis.reduced_values import *


class TestReducedValues(TestCase):
    def test_calc_reduced_equivalent_V(self):
        assert_allclose(calc_reduced_equivalent_V(300, 60500), 300)
        assert_allclose(calc_reduced_equivalent_V(95, 100000), 73.892658634)
        assert_allclose(calc_reduced_equivalent_V(245, 20000), 426.116914708)
        
    def test_calc_reduced_elevator_deflection(self):
        assert_allclose(calc_reduced_elevator_deflection(3.0, -0.04, 0.01, 0.011), 3.00016)
        assert_allclose(calc_reduced_elevator_deflection(3, -0.04, 0.01, 0.01), 3.0)
        assert_allclose(calc_reduced_elevator_deflection(3, -0.04, 0.5, 0.1), 2.936)
        
    def test_calc_reduced_stick_force(self):
        assert_allclose(calc_reduced_stick_force(20, 60500), 20)
        assert_allclose(calc_reduced_stick_force(100, 5000), 1210)
        assert_allclose(calc_reduced_stick_force(1, 100000), 0.605)
