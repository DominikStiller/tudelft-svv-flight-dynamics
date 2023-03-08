import numpy as np
from numpy.typing import ArrayLike

from fd.structs import AerodynamicParameters


class AircraftModel:
    def __init__(self, aero_params: AerodynamicParameters):
        self.aero_params = aero_params

    def get_state_space_matrices(
        self, mass: float
    ) -> tuple[ArrayLike, ArrayLike, ArrayLike, ArrayLike]:
        A = np.array()
        B = np.array()
        C = np.array()
        D = np.array()
        return A, B, C, D
