import numpy as np
import numpy.linalg as alg
from numpy.typing import ArrayLike
from math import sin, cos, pi
import control.matlab as ml
from B24.fd.structs import AerodynamicParameters
from B24.fd.simulation.constants_Cessna_Ce500 import *


class AircraftModel:
    def __init__(self, aero_params: AerodynamicParameters):
        self.aero_params = aero_params

    def get_non_dim_masses(self, m: float, rho: float):
        """

        Args:
            m: Aircraft mass
            rho: Air density for the initial steady conditions

        Returns:
            muc: Non-dimensional aircraft mass wrt MAC
            mub: Non-dimensional aircraft mass wrt wingspan

        """
        muc = m / (rho * S * c)
        mub = m / (rho * S * b)
        return muc, mub

    def get_gravity_term_coeff(self, m: float, V0: float, rho: float, th0: float):
        """

        Args:
            m: Aircraft mass
            V0: Airspeed for the initial steady flight condition
            rho: Air density for the initial steady conditions
            th0: Pitch angle for the initial steady flight condition

        Returns:
            CX0: Gravity term coefficient in X-direction
            CZ0: Gravity term coefficient in Z-direction

        """
        W = m * g
        CX0 = W * sin(th0) / (0.5 * rho * V0**2 * S)
        CZ0 = -W * cos(th0) / (0.5 * rho * V0**2 * S)
        return CX0, CZ0

    def get_state_space_matrices_symmetric(
        self, m: float, V0: float, rho: float, th0: float
    ) -> tuple[ArrayLike, ArrayLike, ArrayLike, ArrayLike]:
        """

        Args:
            m: Aircraft mass
            V0: Airspeed for the initial steady flight condition
            rho: Air density for the initial steady conditions
            th0: Pitch angle for the initial steady flight condition

        Returns:
            A: State matrix
            B: Control matrix
            C: Output matrix
            D: Feedthrough matrix

        """
        Cma = self.aero_params.C_m_alpha
        Cmde = self.aero_params.C_m_delta
        muc = self.get_non_dim_masses(m, rho)[0]
        CX0, CZ0 = self.get_gravity_term_coeff(m, V0, rho, th0)

        # C_1*x_dot + C_2*x +C_3*u = 0
        C_1 = np.array(
            [
                [-2 * muc, 0, 0, 0],
                [0, (CZadot - 2 * muc) * c / V0, 0, 0],
                [0, 0, -c / V0, 0],
                [0, Cmadot * c / V0, 0, -2 * muc * (KY2**2) * ((c / V0) ** 2)],
            ]
        )
        C_2 = np.array(
            [
                [CXu, CXa, CZ0, CXq * c / V0],
                [CZu, CZa, -CX0, (CZq + 2 * muc) * c / V0],
                [0, 0, 0, c / V0],
                [Cmu, Cma, 0, Cmq * c / V0],
            ]
        )
        C_3 = np.array([[CXde], [CZde], [0], [Cmde]])
        A = np.matmul(-alg.inv(C_1), C_2)
        B = np.matmul(-alg.inv(C_1), C_3)
        # In order to get the state variables as output:
        C = np.eye(4)
        D = np.zeros((4, 1))
        return A, B, C, D

    def get_state_space_matrices_asymmetric(
        self, m: float, V0: float, rho: float, th0: float, CL: float
    ) -> tuple[ArrayLike, ArrayLike, ArrayLike, ArrayLike]:
        """

        Args:
            m: Aircraft mass
            V0: Airspeed for the initial steady flight condition
            rho: Air density for the initial steady conditions
            th0: Pitch angle for the initial steady flight condition
            CL: Lift coefficient for steady flight

        Returns:
            A: State matrix
            B: Control matrix
            C: Output matrix
            D: Feedthrough matrix
        """

        mub = self.get_non_dim_masses(m, rho)[-1]

        # C_1*x_dot + C_2*x +C_3*u = 0
        C_1 = np.array(
            [
                [(CYbdot - 2 * mub) * b / V0, 0, 0, 0],
                [0, -b / 2 * V0, 0, 0],
                [
                    0,
                    0,
                    -2 * mub * KX2**2 * (b / V0) ** 2,
                    2 * mub * KXZ * (b / V0) ** 2,
                ],
                [
                    Cnbdot * b / V0,
                    0,
                    2 * mub * KXZ * (b / V0) ** 2,
                    -2 * mub * KZ2**2 * (b / V0) ** 2,
                ],
            ]
        )
        C_2 = np.array(
            [
                [CYb, CL, CYp * b / (2 * V0), (CYr - 4 * mub) * b / (2 * V0)],
                [0, 0, b / (2 * V0), 0],
                [Clb, 0, Clp * b / (2 * V0), Clr],
                [Cnb, 0, Clp * b / (2 * V0), Cnr * b / (2 * V0)],
            ]
        )
        C_3 = np.array([[CYda, CYdr], [0, 0], [Clda, Cldr], [Cnda, Cndr]])
        A = np.matmul(-alg.inv(C_1), C_2)
        B = np.matmul(-alg.inv(C_1), C_3)
        # In order to get the state variables as output:
        C = np.eye(4)
        D = np.zeros((2, 2))
        return A, B, C, D

    def get_eigenvalue_and_eigenvector(
            self, A: ArrayLike):
        '''

        Args:
            A: State matrix
            B: Control matrix
            C: Output matrix
            D: Feedthrough matrix

        Returns:


        '''
        eigenvalues, eigenvectors = alg.eig(A)
        return eigenvalues, eigenvectors

