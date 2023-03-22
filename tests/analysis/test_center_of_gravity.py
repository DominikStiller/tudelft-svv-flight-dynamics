from unittest import TestCase

from numpy.testing import assert_allclose

from fd.analysis.center_of_gravity import *


class TestAerodynamics(TestCase):
    def test_lin_moment_mass(self):
        assert_allclose(lin_moment_mass(), [7.238938216, 6.598248836])

    def test_get_cg(self):
        assert_allclose(calc_cg_position(1000, 80, 80, 80, 80, 80, 80, 80, 80, 80), 7.14458657)
        assert_allclose(calc_cg_position(1000, 0, 0, 0, 0, 0, 0, 0, 0, 0), 7.37828993)
        assert_allclose(calc_cg_position(1000, 80, 80, 80, 80, 80, 80, 80, 80, 80, True), 7.09932165)
