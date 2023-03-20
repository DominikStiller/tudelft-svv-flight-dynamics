import unittest
from fd.simulation.aircraft_model import AircraftModel
from fd.structs import AerodynamicParameters
from tests.test_simulation.constants_Cessna_Ce500 import *
import numpy as np


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
        eigenvalues = model.get_eigenvalues_and_eigenvectors(A)[0]
        first = (eigenvalues[0] == np.conj(eigenvalues[1]))
        second = (eigenvalues[2] == np.conj(eigenvalues[3]))


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
        eigenvalues = model.get_eigenvalues_and_eigenvectors(A)[0]
        first = (eigenvalues[0]<0)
        second = (eigenvalues[1] == np.conj(eigenvalues[2]))
        third = (eigenvalues[3]>0)

        self.assertTupleEqual((first, second, third), (True, True, True))

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

        A_prim = 4 * (muc ** 2) * (KY2) * (CZadot - 2 * muc)
        B_prim = (
                Cmadot * 2 * muc * (CZq + 2 * muc)
                - Cmq * 2 * muc * (CZadot - 2 * muc)
                - 2 * muc * (KY2) * (CXu * (CZadot - 2 * muc) - 2 * muc * CZa)
        )
        C_prim = (
                Cma * 2 * muc * (CZq + 2 * muc)
                - Cmadot * (2 * muc * CX0 + CXu * (CZq + 2 * muc))
                + Cmq * (CXu * (CZadot - 2 * muc) - 2 * muc * CZa)
                + 2 * muc * (KY2) * (CXa * CZu - CZa * CXu)
        )
        D_prim = (
                Cmu * (CXa * (CZq + 2 * muc) - CZ0 * (CZadot - 2 * muc))
                - Cma * (2 * muc * CX0 + CXu * (CZq + 2 * muc))
                + Cmadot * (CX0 * CXu - CZ0 * CZu)
                + Cmq * (CXu * CZa - CZu * CXa)
        )
        E_prim = -Cmu * (CX0 * CXa + CZ0 * CZa) + Cma * (CX0 * CXu + CZ0 * CZu)
        p = (A_prim, B_prim, C_prim, D_prim, E_prim)
        roots = np.polynomial.polynomial.polyroots(p)*c/V0
        print(roots)
        print(eigenvalues)

        self.assertTupleEqual(tuple(roots), tuple(eigenvalues))

if __name__ == '__main__':
    unittest.main()
