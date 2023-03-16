import numpy as np
import numpy.linalg as alg
from numpy.typing import ArrayLike
from math import sin, cos, pi
import control.matlab as ml
from fd.structs import AerodynamicParameters
from fd.simulation.constants_Cessna_Ce500 import *
import matplotlib.pyplot as plt


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
        #x = [u, alpha, theta, q]T
        C_1 = np.array(
            [
                [-2 * muc*c/V0**2, 0, 0, 0],
                [0, (CZadot - 2 * muc) * c / V0, 0, 0],
                [0, 0, -c / V0, 0],
                [0, Cmadot * c / V0, 0, -2 * muc * (KY2) * ((c / V0) ** 2)],
            ]
        )
        C_2 = np.array(
            [
                [CXu/V0, CXa, CZ0, CXq * c / V0],
                [CZu/V0, CZa, -CX0, (CZq + 2 * muc) * c / V0],
                [0, 0, 0, c / V0],
                [Cmu/V0, Cma, 0, Cmq * c / V0],
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

        #mub = self.get_non_dim_masses(m, rho)[-1]

        # C_1*x_dot + C_2*x +C_3*u = 0
        C_1 = np.array(
            [
                [(CYbdot - 2 * mub) * b / V0, 0, 0, 0],
                [0, -b / (2 * V0), 0, 0],
                [
                    0,
                    0,
                    -2 * mub * KX2 * (b**2 / V0**2),
                    2 * mub * KXZ * (b**2 / V0**2),
                ],
                [
                    Cnbdot * b / V0,
                    0,
                    2 * mub * KXZ * (b**2 / V0**2),
                    -2 * mub * KZ2 * (b**2 / V0**2),
                ],
            ]
        )
        C_2 = np.array(
            [
                [CYb, CL, CYp * b / (2 * V0), (CYr - 4 * mub) * b / (2 * V0)],
                [0, 0, b / (2 * V0), 0],
                [Clb, 0, Clp * b / (2 * V0), Clr * b / (2 * V0)],
                [Cnb, 0, Cnp * b / (2 * V0), Cnr * b / (2 * V0)],
            ]
        )
        C_3 = np.array([[CYda, CYdr], [0, 0], [Clda, Cldr], [Cnda, Cndr]])
        A = np.matmul(-alg.inv(C_1), C_2)
        B = np.matmul(-alg.inv(C_1), C_3)
        # In order to get the state variables as output:
        C = np.eye(4)
        D = np.zeros((4, 2))
        return A, B, C, D

    def get_eigenvalues_and_eigenvectors(
            self, A: ArrayLike):
        '''

        Args:
            A: State matrix
            B: Control matrix
            C: Output matrix
            D: Feedthrough matrix

        Returns:
            Eigenvalues and Eigenvectors
        '''
        eigenvalues, eigenvectors = alg.eig(A)
        return eigenvalues, eigenvectors

    def get_state_space_matrices_asymmetric_reader(
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

        #mub = self.get_non_dim_masses(m, rho)[-1]

        # C_1*x_dot + C_2*x +C_3*u = 0
        # X = [beta, phi, pb/2V, rb/2V]T
        C_1 = np.array(
            [
                [(CYbdot - 2 * mub) * b / V0, 0, 0, 0],
                [0, -b / (2 * V0), 0, 0],
                [
                    0,
                    0,
                    -4 * mub * KX2 * (b / V0),
                    4 * mub * KXZ * (b / V0),
                ],
                [
                    Cnbdot * b / V0,
                    0,
                    4 * mub * KXZ * (b / V0),
                    4 * mub * KZ2 * (b / V0),
                ],
            ]
        )
        C_2 = np.array(
            [
                [CYb, CL, CYp , (CYr - 4 * mub) ],
                [0, 0, 1, 0],
                [Clb, 0, Clp , Clr ],
                [Cnb, 0, Cnp , Cnr ],
            ]
        )
        C_3 = np.array([[CYda, CYdr], [0, 0], [Clda, Cldr], [Cnda, Cndr]])
        A = np.matmul(-alg.inv(C_1), C_2)
        B = np.matmul(-alg.inv(C_1), C_3)
        # In order to get the state variables as output:
        C = np.eye(4)
        D = np.zeros((4, 2))
        return A, B, C, D

    def get_step_input(self, maneuvre_duration, dt, input_duration, input_value, plot = False):
        t = np.arange(0, maneuvre_duration + dt, dt)
        u = np.zeros(t.shape)
        u[:int(input_duration / dt)] = input_value * np.ones(u[:int(input_duration / dt)].size)
        if plot:
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            ax.plot(t, u)
            ax.set_xlabel("Time [s]")
            ax.set_ylabel("$delta_e$")
        return t, u
    def get_response_plots_symmetric(self, sys, x0, t, u ):

        yout, t, xout = ml.lsim(sys, u, t, x0)
        fig, axs = plt.subplots(2, 2, sharex=True)

        axs[0, 0].plot(t, xout[:, 0] + V0 * np.ones(t.size))
        axs[0, 0].set_title("V [m/sec")

        axs[1, 0].plot(t, xout[:, 1])
        axs[1, 0].set_title("$alpha$ [rad]")

        axs[0, 1].plot(t, xout[:, 2])
        axs[0, 1].set_title("$theta$ [rad]")

        axs[1, 1].plot(t, xout[:, 3])
        axs[1, 1].set_title("q [rad/sec]")

        plt.show()

    def get_response_plots_asymmetric(self, sys, x0, t, u):
        yout, t, xout = ml.lsim(sys, u, t, x0)
        fig, axs = plt.subplots(2, 2, sharex=True)

        axs[0, 0].plot(t, xout[:, 0] + V0 * np.ones(t.size))
        axs[0, 0].set_title("$beta$ [rad]")

        axs[1, 0].plot(t, xout[:, 1])
        axs[1, 0].set_title("$phi$ [rad]")

        axs[0, 1].plot(t, xout[:, 2])
        axs[0, 1].set_title("p [rad/sec]")

        axs[1, 1].plot(t, xout[:, 3])
        axs[1, 1].set_title("r [rad/sec]")

        plt.show()


