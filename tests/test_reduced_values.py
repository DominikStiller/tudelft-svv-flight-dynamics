from unittest import TestCase
import numpy as np
from numpy.testing import assert_allclose

from fd.analysis import reduced_values as red


class Testaeroparameters(TestCase):
    def test_aeroparameters(self):
        assert_allclose(red.calc_reduced_equivalent_V(300, 60500), 300)
        assert_allclose(red.calc_reduced_equivalent_V(95, 100000), 73.892658634)
        assert_allclose(red.calc_reduced_equivalent_V(245, 20000), 426.116914708)
        assert_allclose(red.calc_reduced_elev_deflec(3.0, -0.04, 0.01, 0.011), 3.00016)
        assert_allclose(red.calc_reduced_elev_deflec(3, -0.04, 0.01, 0.01), 3.0)
        assert_allclose(red.calc_reduced_elev_deflec(3, -0.04, 0.5, 0.1), 2.936)
        assert_allclose(red.calc_reduced_stick_force(20, 60500), 20)
        assert_allclose(red.calc_reduced_stick_force(100, 5000), 1210)
        assert_allclose(red.calc_reduced_stick_force(1, 100000), 0.605)
