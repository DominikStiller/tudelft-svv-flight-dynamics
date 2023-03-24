import math

import numpy as np
import scipy.stats as stats

from fd.simulation import constants


def calc_true_V(T, M):
    """

    Args:
        T (float): Static temperature[K]
        M (float): Mach number[-]

    Returns (float): True airspeed[m/s]

    """
    return M * np.sqrt(constants.gamma * constants.R * T)


def calc_dynamic_pressure(V: float, rho: float) -> float:
    """
    Calculate dynamic pressure

    Args:
        V: velocity (CAS or TAS) [m/s]
        rho: air density (sea level or actual) [kg/m^3]

    Returns:
        Dynamic pressure [Pa]
    """
    return rho * V**2 / 2


def calc_equivalent_V(Vt, rho):
    """

    Args:
        Vt (float): True airspeed[m/s]
        rho (float): Density[kg/m^3]

    Returns (float): Equivalent airspeed[m/s]

    """
    return Vt * np.sqrt(rho / constants.rho0)


def calc_CL(W: float, V: float, rho: float, S=constants.S) -> float:
    """
    Calculate CL for a given combination of W, rho, V and S.
    Args:
        W (array_like): Weight [N]
        rho (float): Air density [kg/m3]
        V (array_like): True airspeed [m/s]
        S (float): Surface area [m2]

    Returns:
        (array_like): CL [-]
    """

    return W / (calc_dynamic_pressure(V, rho) * S)


def estimate_CL_alpha(CL: float, alpha: float) -> tuple[float, float, float]:
    """
    Calculate the slope, CL-intercept and alpha intercept of the CL-alpha plot using a Theil-Sen robust linear estimator.
    Args:
        CL (array_like): CL [-]
        alpha (array_like): angle of attack [deg]

    Returns:
        CLalpha (float): slope of CL-alpha plot [1/deg]
        CL_alpha_equals0 (float): CL at alpha = 0
        alpha_0 (float): alpha at CL = 0
    """
    CLalpha, CL_alpha_equals0, _, _ = stats.theilslopes(CL, alpha, alpha=0.99)
    alpha_0 = -CL_alpha_equals0 / CLalpha

    return CLalpha, CL_alpha_equals0, alpha_0


def calc_CD(T: float, V: float, rho: float, S: float = constants.S) -> float:
    """
    This function calculates the drag coefficient CD[-] based on the thrust.

    Args:
        T (array_like): Thrust[N].
        V (array_like): True airspeed[m/s].

    Returns:
        CD (array_like): Drag coefficient CD[-].

    """
    return T / (calc_dynamic_pressure(V, rho) * S)


def estimate_CD0_e(CD: list, CL: list) -> tuple[float, float]:
    """
    This function uses the parabolic drag formula to calculate the zero lift drag, CD0[-], and the oswald
    efficiency factor, e[-].

    Args:
        CD (list): The drag coefficient CD[-].
        CL (list): The lift coefficient CL[-].

    Returns:
        CD0, e (float, float): Zero lift drag coefficient, CD0[-], oswald efficiency factor, e[-].

    """

    slope, CD0, _, _ = stats.theilslopes(CD, CL**2, alpha=0.99)
    # slope, CD0, _, _, _ = stats.linregress(CL**2, CD)
    e = 1 / (math.pi * constants.A * slope)

    return CD0, e


def estimate_Cmalpha(alpha, delta_e, Cmdelta):
    """

    Args:
        alpha (array_like): Angle of attack[deg]
        delta_e (array_like): Elevator deflection[deg]
        Cmdelta (float): Change in moment coefficient due to elevator deflection[-]

    Returns (float): Change in moment coefficient due to angle of attack[-]

    """

    slope, _, _, _ = stats.theilslopes(delta_e, alpha, alpha=0.99)
    return -slope * Cmdelta


def calc_Cmdelta(
    xcg1: float,
    xcg2: float,
    deltae1: float,
    deltae2: float,
    W: float,
    V: float,
    rho: float,
):
    """

    Args:
        xcg1 (float): X-position of the center of gravity during the first test.(aft cg)[m]
        xcg2 (float): X-position of the center of gravity during the second test.(front cg)[m]
        deltae1 (float): Deflection of the elevator during the first test.[deg]
        deltae2 (float): Deflection of the elevator during the second test.[deg]
        W (float): Weight of the aircraft during the tests.[N]
        V (float): Velocity of the aircraft during the tests.[m/s]
        rho (float): Air density.[kg/m^3]

    Returns: Cmdelta (float): The moment coefficient change due to the elevator deflection.[-]

    """
    Delta_cg = xcg2 - xcg1
    Delta_delta_e = deltae2 - deltae1
    C_N = W / (calc_dynamic_pressure(V, rho) * constants.S)
    return -1 / Delta_delta_e * C_N * Delta_cg / constants.c
