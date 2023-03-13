import numpy as np
import numpy.linalg as alg
from numpy.typing import ArrayLike

from B24.fd.structs import AerodynamicParameters
from B24.fd.simulation.constants import *


class AircraftModel:
    def __init__(self, aero_params: AerodynamicParameters):
        self.aero_params = aero_params



    def get_state_space_matrices_symmetric(
        self, mass: float
    ) -> tuple[ArrayLike, ArrayLike, ArrayLike, ArrayLike]:
        # C_1*x_dot + C_2*x +C_3*u = 0
        C_1 = np.array([[-2 * muc * c / (V0 ** 2), 0, 0, 0],
                        [0, (CZadot - 2 * muc) * c / V0, 0, 0],
                        [0, 0, -c / V0, 0],
                        [0, Cmadot * c / V0, 0, -2 * muc * (KY2 ** 2) * ((c / V0) ** 2)]])
        C_2 = np.array([[CXu / V0, CXa, CZ0, CXq * c / V0],
                        [CZu / V0, CZa, -CX0, (CZq + 2 * muc) * c / V0],
                        [0, 0, 0, c / V0],
                        [Cmu / V0, Cma, 0, Cmq * c / V0]])
        C_3 = np.array([[CXde],
                        [CZde],
                        [0],
                        [Cmde]])
        A = np.matmul(-alg.inv(C_1),C_2)
        B = np.matmul(-alg.inv(C_1), C_3)
        #In order to get the state variables as output:
        C = np.eye(4)
        D = np.zeros((1, 1))
        return A, B, C, D

    def get_state_space_matrices_asymmetric(
            self, mass: float
    ) -> tuple[ArrayLike, ArrayLike, ArrayLike, ArrayLike]:
        # C_1*x_dot + C_2*x +C_3*u = 0
        C_1 = np.array([[(CYbdot-2*mub)*b/V0, 0, 0, 0],
                        [0, -b/2*V0, 0, 0],
                        [0, 0, -2*mub*KX2**2*(b/V0)**2, 2*mub*KXZ*(b/V0)**2],
                        [Cnbdot*b/V0, 0, 2*mub*KXZ*(b/V0)**2, -2*mub*KZ2**2*(b/V0)**2]])
        C_2 = np.array([[CYb, CL, CYp*b/(2*V0), (CYr-4*mub)*b/(2*V0)],
                        [0, 0, b/(2*V0), 0],
                        [Clb, 0, Clp*b/(2*V0), Clr],
                        [Cnb, 0, Clp*b/(2*V0), Cnr*b/(2*V0)]])
        C_3 = np.array([[CYda, CYdr],
                        [0, 0],
                        [Clda, Cldr],
                        [Cnda, Cndr]])
        A = np.matmul(-alg.inv(C_1), C_2)
        B = np.matmul(-alg.inv(C_1), C_3)
        #In order to get the state variables as output:
        C = np.eye(4)
        D = np.zeros((2, 2))
        return A, B, C, D

