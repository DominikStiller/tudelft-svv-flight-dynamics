import unittest
from unittest import skip

import numpy as np
from numpy.testing import assert_allclose

from fd.simulation.aircraft_model import AircraftModel
from fd.structs import AerodynamicParameters
from tests.test_simulation.constants_Cessna_Ce500 import *


class TestEigenvalues(unittest.TestCase):
    def test_type_eigenvalues_symmetric(self):
        aero_params = AerodynamicParameters
        aero_params.C_m_alpha = -0.4300
        aero_params.C_m_delta = -1.5530
        m = 4547.8
        V0 = 59.9
        rho = 0.904627056
        th0 = 0
        model = AircraftModel(aero_params)
        A, B, C, D = model.get_state_space_matrices_symmetric(m, V0, rho, th0)
        eigenvalues, eigenvectors = model.get_eigenvalues_and_eigenvectors(A)
        first = eigenvalues[0] == np.conj(eigenvalues[1])
        second = eigenvalues[2] == np.conj(eigenvalues[3])

        self.assertTupleEqual((first, second), (True, True))

    @skip
    def test_type_eigenvalues_asymmetric(self):
        aero_params = AerodynamicParameters
        aero_params.C_m_alpha = -0.4300
        aero_params.C_m_delta = -1.5530
        m = 4547.8
        V0 = 59.9
        rho = 0.904627056
        th0 = 0
        CL = 1.1360
        model = AircraftModel(aero_params)
        A, B, C, D = model.get_state_space_matrices_asymmetric(m, V0, rho, th0, CL)
        eigenvalues, eigenvectors = model.get_eigenvalues_and_eigenvectors(A)
        self.assertTrue(eigenvalues[0] < 0)
        self.assertTrue(eigenvalues[1] == np.conj(eigenvalues[2]))
        self.assertTrue(eigenvalues[3] > 0)

    @skip
    def test_shortperiod_eigenvalues(self):
        aero_params = AerodynamicParameters
        aero_params.C_m_alpha = -0.4300
        aero_params.C_m_delta = -1.5530
        m = 4547.8
        V0 = 59.9
        rho = 0.904627056
        th0 = 0
        model = AircraftModel(aero_params)
        A, B, C, D = model.get_state_space_matrices_symmetric(m, V0, rho, th0)
        print(model.get_eigenvalues_and_eigenvectors(A)[0])
        eig1, eig2 = model.get_idealized_shortperiod_eigenvalues(m, rho, V0)
        eigenvalues2 = complex(-0.039161, -0.037971) * V0 / c
        eigenvalues1 = complex(-0.039161, -0.037971) * V0 / c
        self.assertAlmostEqual(eig1, eigenvalues2)
        self.assertAlmostEqual(eig2, eigenvalues1)

    @skip
    def test_phugoid_eigenvalues(self):
        aero_params = AerodynamicParameters
        aero_params.C_m_alpha = -0.4300
        aero_params.C_m_delta = -1.5530
        m = 4547.8
        V0 = 59.9
        rho = 0.904627056
        th0 = 0
        model = AircraftModel(aero_params)
        A, B, C, D = model.get_state_space_matrices_symmetric(m, V0, rho, th0)
        # print(model.get_eigenvalues_and_eigenvectors(A)[0])
        eig1, eig2 = model.get_idealized_phugoid_eigenvalues(m, rho, V0, th0)
        eigenvalues2 = complex(-0.00029107, 0.0066006) * V0 / c
        eigenvalues1 = complex(-0.00029107, -0.0066006) * V0 / c
        self.assertAlmostEqual(eig1, eigenvalues2)
        self.assertAlmostEqual(eig2, eigenvalues1)

    @skip
    def test_aperiodicroll_eigenvalues(self):
        aero_params = AerodynamicParameters
        aero_params.C_m_alpha = -0.4300
        aero_params.C_m_delta = -1.5530
        m = 4547.8
        V0 = 59.9
        rho = 0.904627056
        th0 = 0
        CL = 1.1360
        model = AircraftModel(aero_params)
        A, B, C, D = model.get_state_space_matrices_asymmetric(m, V0, rho, th0, CL)

        A_prim = 4 * muc**2 * KY2 * (CZadot - 2 * muc)
        B_prim = (
            Cmadot * 2 * muc * (CZq + 2 * muc)
            - Cmq * 2 * muc * (CZadot - 2 * muc)
            - 2 * muc * KY2 * (CXu * (CZadot - 2 * muc) - 2 * muc * CZa)
        )
        C_prim = (
            Cma * 2 * muc * (CZq + 2 * muc)
            - Cmadot * (2 * muc * CX0 + CXu * (CZq + 2 * muc))
            + Cmq * (CXu * (CZadot - 2 * muc) - 2 * muc * CZa)
            + 2 * muc * KY2 * (CXa * CZu - CZa * CXu)
        )
        D_prim = (
            Cmu * (CXa * (CZq + 2 * muc) - CZ0 * (CZadot - 2 * muc))
            - Cma * (2 * muc * CX0 + CXu * (CZq + 2 * muc))
            + Cmadot * (CX0 * CXu - CZ0 * CZu)
            + Cmq * (CXu * CZa - CZu * CXa)
        )
        E_prim = -Cmu * (CX0 * CXa + CZ0 * CZa) + Cma * (CX0 * CXu + CZ0 * CZu)
        p = (E_prim, D_prim, C_prim, B_prim, A_prim)

        roots = np.polynomial.polynomial.polyroots(p)
        eig1 = model.get_aperiodicroll_eigenvalues(m, rho, V0, A)
        eigenvalues1 = -0.3291 * V0 / b

        self.assertAlmostEqual(roots[2], eig1)

    @skip
    def test_dutchroll_eigenvalues(self):
        aero_params = AerodynamicParameters
        aero_params.C_m_alpha = -0.4300
        aero_params.C_m_delta = -1.5530
        m = 4547.8
        V0 = 59.9
        rho = 0.904627056
        th0 = 0
        CL = 1.1360
        model = AircraftModel(aero_params)
        A, B, C, D = model.get_state_space_matrices_asymmetric(m, V0, rho, th0, CL)
        # print(model.get_eigenvalues_and_eigenvectors(A)[0])
        eig1, eig2 = model.get_dutchroll_eigenvalues(m, rho, V0, A)
        eigenvalues1 = complex(-0.0313, 0.3314) * V0 / b
        eigenvalues2 = complex(-0.0313, -0.3314) * V0 / b
        self.assertAlmostEqual(eig1, eigenvalues1)
        self.assertAlmostEqual(eig2, eigenvalues2)

    @skip
    def test_spiral_eigenvalues(self):
        aero_params = AerodynamicParameters
        aero_params.C_m_alpha = -0.4300
        aero_params.C_m_delta = -1.5530
        m = 4547.8
        V0 = 59.9
        rho = 0.904627056
        th0 = 0
        CL = 1.1360
        model = AircraftModel(aero_params)
        A, B, C, D = model.get_state_space_matrices_asymmetric(m, V0, rho, th0, CL)
        # print(model.get_eigenvalues_and_eigenvectors(A)[0])
        eig1 = model.get_spiral_eigenvalues(m, rho, V0, CL, A)
        eigenvalues1 = -0.0108 * V0 / b
        self.assertAlmostEqual(eig1, eigenvalues1)

    @skip
    def test_analytic_eigenvalues_symmetric(self):
        # In order to perform this test you need to:
        # 1. Change the imported constants file in aircraft model with the ones for cessna Ce500
        # 2. Comment any mub calculation out from the aircraft model
        aero_params = AerodynamicParameters
        aero_params.C_m_alpha = -0.4300
        aero_params.C_m_delta = -1.5530
        m = 4547.8
        V0 = 59.9
        rho = 0.904627056
        th0 = 0
        model = AircraftModel(aero_params)
        A, B, C, D = model.get_state_space_matrices_symmetric(m, V0, rho, th0)
        eigenvalues = model.get_eigenvalues_and_eigenvectors(A)[0]

        A_prim = 4 * muc**2 * KY2 * (CZadot - 2 * muc)
        B_prim = (
            Cmadot * 2 * muc * (CZq + 2 * muc)
            - Cmq * 2 * muc * (CZadot - 2 * muc)
            - 2 * muc * KY2 * (CXu * (CZadot - 2 * muc) - 2 * muc * CZa)
        )
        C_prim = (
            Cma * 2 * muc * (CZq + 2 * muc)
            - Cmadot * (2 * muc * CX0 + CXu * (CZq + 2 * muc))
            + Cmq * (CXu * (CZadot - 2 * muc) - 2 * muc * CZa)
            + 2 * muc * KY2 * (CXa * CZu - CZa * CXu)
        )
        D_prim = (
            Cmu * (CXa * (CZq + 2 * muc) - CZ0 * (CZadot - 2 * muc))
            - Cma * (2 * muc * CX0 + CXu * (CZq + 2 * muc))
            + Cmadot * (CX0 * CXu - CZ0 * CZu)
            + Cmq * (CXu * CZa - CZu * CXa)
        )
        E_prim = -Cmu * (CX0 * CXa + CZ0 * CZa) + Cma * (CX0 * CXu + CZ0 * CZu)
        p = (E_prim, D_prim, C_prim, B_prim, A_prim)
        roots = np.polynomial.polynomial.polyroots(p)

        assert_allclose(roots * V0 / c, np.sort(eigenvalues), rtol=1e-8)

    @skip
    def test_analytic_eigenvalues_asymmetric(self):
        # In order to perform this test you need to:
        # 1. Change the imported constants file in aircraft model with the ones for cessna Ce500
        # 2. Comment any mub calculation out from the aircraft model
        aero_params = AerodynamicParameters
        aero_params.C_m_alpha = -0.4300
        aero_params.C_m_delta = -1.5530
        m = 4547.8
        V0 = 59.9
        rho = 0.904627056
        th0 = 0
        CL = 1.1360
        model = AircraftModel(aero_params)
        A, B, C, D = model.get_state_space_matrices_asymmetric(m, V0, rho, th0, CL)
        eigenvalues = model.get_eigenvalues_and_eigenvectors(A)[0]
        A_prim = 16 * mub**3 * (KX2 * KZ2 - KXZ**2)
        B_prim = (
            -4
            * mub**2
            * (2 * CYb * (KX2 * KZ2 - KXZ**2) + Cnr * KX2 + Clp * KZ2 + (Clr + Cnp) * KXZ)
        )
        C_prim = (
            2
            * mub
            * (
                (CYb * Cnr - CYr * Cnb) * KX2
                + (CYb * Clp - Clb * CYp) * KZ2
                + ((CYb * Cnp - Cnb * CYp) + (CYb * Clr - Clb * CYr)) * KXZ
                + 4 * mub * Cnb * KX2
                + 4 * mub * Clb * KXZ
                + 0.5 * (Clp * Cnr - Cnp * Clr)
            )
        )
        D_prim = (
            -4 * mub * CL * (Clb * KZ2 + Cnb * KXZ)
            + 2 * mub * (Clb * Cnp - Cnb * Clp)
            + 0.5 * CYb * (Clr * Cnp - Cnr * Clp)
            + 0.5 * CYp * (Clb * Cnr - Cnb * Clr)
            + 0.5 * CYr * (Clp * Cnb - Cnp * Clb)
        )
        E_prim = CL * (Clb * Cnr - Cnb * Clr)
        p = (E_prim, D_prim, C_prim, B_prim, A_prim)
        roots = np.polynomial.polynomial.polyroots(p)

        assert_allclose(roots * V0 / b, np.sort(eigenvalues), rtol=1e-8)


if __name__ == "__main__":
    unittest.main()
