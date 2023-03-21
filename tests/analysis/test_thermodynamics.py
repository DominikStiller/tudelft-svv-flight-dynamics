from unittest import TestCase

from numpy.testing import assert_allclose

from fd.analysis.thermodynamics import *


class TestThermodynamics(TestCase):
    def test_calc_stat_pres(self):
        assert_allclose(calc_static_pressure(1000), 89870.773519)
        assert_allclose(calc_static_pressure(0), 101325)
        assert_allclose(calc_static_pressure(10672), 23815.2625371)

    def test_calc_mach(self):
        assert_allclose(calc_mach(1000, 100), 0.3116119528)
        assert_allclose(calc_mach(10672, 150), 0.85210191358)
        assert_allclose(calc_mach(0, 20), 0.05877270993)

    def test_calc_static_temp(self):
        assert_allclose(calc_static_temperature(350, 0.2), 347.222222222)
        assert_allclose(calc_static_temperature(500, 0.9), 430.292598967)
        assert_allclose(calc_static_temperature(100, 0.5), 95.2380952381)

    def test_calc_density(self):
        assert_allclose(calc_density(100000, 288), 1.20962279123)
        assert_allclose(calc_density(1005000, 600), 5.83522034489)
        assert_allclose(calc_density(80000, 200), 1.3934854555)
