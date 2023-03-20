from unittest import TestCase

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
        assert_allclose(ana.calc_true_V(600, 0.9), 441.937574777)
        assert_allclose(ana.calc_true_V(200, 0.7), 198.452160482)
        assert_allclose(ana.calc_true_V(555, 0.21), 99.1764547914)
        assert_allclose(ana.calc_rho(100000, 288), 1.20962279123)
        assert_allclose(ana.calc_rho(1005000, 600), 5.83522034489)
        assert_allclose(ana.calc_rho(80000, 200), 1.3934854555)
        assert_allclose(ana.calc_equivalent_V(100, 1.225), 100)
        assert_allclose(ana.calc_equivalent_V(250, 0.82), 204.540300904)
        assert_allclose(ana.calc_equivalent_V(75, 0.105), 21.9577516413)
        assert_allclose(ana.calc_reduced_equivalent_V(300, 60500), 300)
        assert_allclose(ana.calc_reduced_equivalent_V(95, 100000), 73.892658634)
        assert_allclose(ana.calc_reduced_equivalent_V(245, 20000), 426.116914708)
        assert_allclose(ana.calc_Tc(2000, 300, 15, 1.225), 2.41874527589e-3)
