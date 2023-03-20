from fd.simulation import constants
import scipy.stats as stats
import math
import numpy as np


def calc_reduced_equivalent_V(Ve, W):
    """

    Args:
        Ve (float): Equivalent velocity[m/s]
        W (float): Weight of the aircraft[N]

    Returns (float): Reduced equivalent airspeed[m/s]

    """
    return Ve * np.sqrt(constants.Ws / W)


def calc_reduced_elev_deflec(delta_e_meas, Cmdelta, Tcs, Tc):
    """

    Args:
        delta_e_meas (float): The measured elevator deflection[deg]
        Cmdelta (float): Change in moment coefficient due to elevator defection[-]
        Tcs (float): Thrust coefficient in standard conditions[-]
        Tc (float): Thrust coefficient for conditions used[-]

    Returns (float): The reduced elevator deflection angle[deg]

    """

    return delta_e_meas - constants.CmTc / Cmdelta * (Tcs - Tc)


def calc_reduced_stick_force(Fe_aer, W):
    """

    Args:
        Fe_aer (float): The measured stick force[N]
        W (float): The actual weight of the aircraft[N]

    Returns (float): The reduced stick force[N]

    """

    return Fe_aer * constants.Ws / W
