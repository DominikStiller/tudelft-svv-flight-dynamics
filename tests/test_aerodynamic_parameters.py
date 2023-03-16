from unittest import TestCase
import numpy as np
from numpy.testing import assert_allclose

from fd.analysis import aerodynamic_analysis as ana


class Testaeroparameters(TestCase):
    def test_aeroparameters(self):
        assert_allclose(ana.calc_stat_pres(1000), 89870.773519)
        assert_allclose(ana.calc_stat_pres(0), 101325)
        assert_allclose(ana.calc_stat_pres(10672), 23815.2625371)
        assert_allclose(ana.calc_mach(1000, 100), 0.3116119528)
        assert_allclose(ana.calc_mach(10672, 150), 0.85210191358)
        assert_allclose(ana.calc_mach(0, 20), 0.05877270993)
        assert_allclose(ana.calc_static_temp(350, 0.2), 347.222222222)
        assert_allclose(ana.calc_static_temp(500, 0.9), 430.292598967)
        assert_allclose(ana.calc_static_temp(100, 0.5), 95.2380952381)
