import unittest
from fd.simulation.aircraft_model import AircraftModel
from fd.structs import AerodynamicParameters
from tests.test_simulation.constants_Cessna_Ce500 import *
import numpy as np
import sympy as sp


class MyTestCase(unittest.TestCase):
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
        first = eigenvalues[0] < 0
        second = eigenvalues[1] == np.conj(eigenvalues[2])
        third = eigenvalues[3] > 0

        self.assertTupleEqual((first, second, third), (True, True, True))

    """
    def test_analytic_eigenvalues_symmetric(self):
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
        l = sp.Symbol('l')
        m = sp.Matrix([
            [CXu-2*muc*l, CXa, CZ0, 0],
            [CZu, CZa+(CZadot-2*muc)*l, -CX0, (CZq+2*muc)*c/V0],
            [0, 0, -l, -c/V0],
            [Cmu, Cma +Cmadot*l, 0, (Cmq-2*muc*KY2*l)*c/V0]])
        fuck = sp.simplify(m.det())
        roots = sp.solveset(fuck, l)
        final_roots = np.array([])
        for root in roots:
            mat = np.array([
                [CXu - 2 * muc * root, CXa, CZ0, 0],
                [CZu, CZa + (CZadot - 2 * muc) * root, -CX0, (CZq + 2 * muc) * c / V0],
                [0, 0, -root, -c / V0],
                [Cmu, Cma + Cmadot * root, 0, (Cmq - 2 * muc * KY2 * root) * c / V0]], dtype=np.complex128)
            if abs(np.linalg.det(mat)) < 1e-10:
                final_roots = np.append(final_roots, np.array([root], dtype=np.complex128))
        #coeff = np.transpose(np.flip(m.det().as_poly().coeffs()))
        #roots = np.polynomial.polynomial.polyroots(coeff) * c / V0

        self.assertTupleEqual(tuple(sorted(final_roots)), tuple(eigenvalues))

    def test_analytic_eigenvalues_asymmetric(self):
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

        self.assertTupleEqual(tuple(roots * V0 / b), tuple(eigenvalues))
    """


if __name__ == "__main__":
    unittest.main()
