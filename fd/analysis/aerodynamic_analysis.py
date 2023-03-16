from fd.simulation import constants
import scipy.stats as stats
import math
import numpy as np


def calc_stat_pres(hp):
    """
    Calculate the static pressure from the pressure height
    Args:
        hp (float): Pressure height[m]

    Returns (float): The static pressure for the pressure height given[Pa]

    """
    return constants.p0 * (1 + constants.Tempgrad * hp / constants.Temp0) ** (
        -constants.g / (constants.Tempgrad * constants.R)
    )


def calc_mach(hp, Vc):
    """

    Args:
        hp (float): Pressure height[m]
        Vc (float): The calibrated speed[m/s]

    Returns (float): Mach number for the conditions given[-]

    """
    return np.sqrt(
        2
        / (constants.gamma - 1)
        * (
            (
                1
                + constants.p0
                / calc_stat_pres(hp)
                * (
                    (
                        1
                        + (constants.gamma - 1)
                        / (2 * constants.gamma)
                        * constants.rho0
                        / constants.p0
                        * Vc**2
                    )
                    ** (constants.gamma / (constants.gamma - 1))
                    - 1
                )
            )
            ** ((constants.gamma - 1) / constants.gamma)
            - 1
        )
    )


def calc_static_temp(Ttot, M):
    """

    Args:
        Ttot (float): Total temperature[K]
        M (float): Mach number[-]

    Returns (float): Static temperature[K]

    """
    return Ttot / (1 + (constants.gamma - 1) / 2 * M**2)


def calc_true_V(T, M):
    """

    Args:
        T (float): Static temperature[K]
        M (float): Mach number[-]

    Returns (float): True airspeed[m/s]

    """
    return M * np.sqrt(constants.gamma * constants.R * T)


def calc_rho(p, T):
    """

    Args:
        p (float): Static pressure[Pa]
        T (float): Static temperature[K]

    Returns (float): Density[kg/m^3]

    """
    return p / (constants.R * T)


def calc_equivalent_V(Vt, rho):
    """

    Args:
        Vt (float): True airspeed[m/s]
        rho (float): Density[kg/m^3]

    Returns (float): Equivalent airspeed[m/s]

    """
    return Vt * np.sqrt(rho / constants.rho0)


def calc_reduced_equivalent_V(Ve, W):
    """

    Args:
        Ve (float): Equivalent velocity[m/s]
        W (float): Weight of the aircraft[N]

    Returns (float): Reduced equivalent airspeed[m/s]

    """
    return Ve * np.sqrt(constants.Ws / W)


def calc_CL(W: float, V: float, rho: float, S=constants.S) -> float:
    """
    Calculate CL for a given combination of W, rho, V and S.
    Args:
        W (array_like): Weight [-]
        rho (float): Air density [kg/m3]
        V (array_like): True airspeed [m/s]
        S (float): Surface area [m2]

    Returns:
        (array_like): CL [-]
    """

    return 2 * W / (rho * V * V * S)


def calc_CL_alpha(CL: float, alpha: float) -> float:
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
    result = stats.theilslopes(CL, alpha, alpha=0.99)
    CLalpha = result[0]
    CL_alpha_equals0 = result[1]
    # low_slope = result.low_slope # Lower bound of the confidence interval on slope
    # high_slope = result.high_slope# Higher bound of the confidence interval on slope
    alpha_0 = -CL_alpha_equals0 / CLalpha

    return CLalpha, CL_alpha_equals0, alpha_0  # low_slope, high_slope


def calc_CD(T: float, V: float, rho) -> float:
    """
    This function calculates the drag coefficient CD[-]

    Args:
        T (array_like): Thrust[N].
        V (array_like): True airspeed[m/s].

    Returns:
        CD (array_like): Drag coefficient CD[-].

    """
    return T / (0.5 * rho * V * V * constants.S)


def calc_CD0_e(CD: list, CL: list) -> float:
    """
    This function uses the parabolic drag formula to calculate the zero lift drag, CD0[-], and the oswald
    efficiency factor, e[-].

    Args:
        CD (list): The drag coefficient CD[-].
        CL (list): The lift coefficient CL[-].

    Returns:
        CD0, e (float, float): Zero lift drag coefficient, CD0[-], oswald efficiency factor, e[-].

    """

    TheilslopesResults = stats.theilslopes(CD, CL**2, alpha=0.99)

    CD0 = TheilslopesResults[1]
    e = 1 / (math.pi * constants.A * TheilslopesResults[0])

    return CD0, e


def calc_Cmdelta(
    xcg1: float,
    xcg2: float,
    deltae1: float,
    deltae2: float,
    W1: float,
    W2: float,
    V: float,
    rho: float,
):
    """

    Args:
        xcg1 (float): X-position of the center of gravity during the first test.(aft cg)[m]
        xcg2 (float): X-position of the center of gravity during the second test.(front cg)[m]
        deltae1 (float): Deflection of the elevator during the first test.[deg]
        deltae2 (float): Deflection of the elevator during the second test.[deg]
        W1 (float): Weight of the aircraft during the first test.[N]
        W2 (float): Weight of the aircraft during the second test.[N]
        V (float): Velocity of the aircraft during the tests.[m/s]
        rho (float): Air density.[kg/m^3]

    Returns: Cmdelta (float): The moment coefficient change due to the elevator deflection.[-]

    """
    W_avg = (W1 + W2) / 2
    Delta_cg = xcg2 - xcg1
    Delta_delta_e = deltae2 - deltae1
    return -1 / Delta_delta_e * W_avg / (0.5 * rho * V**2 * constants.S) * Delta_cg / constants.c


def calc_Cmalpha(alpha, delta_e, Cmdelta):
    """

    Args:
        alpha (array_like): Angle of attack[deg]
        delta_e (array_like): Elevator deflection[deg]
        Cmdelta (float): Change in moment coefficient due to elevator deflection[-]

    Returns (float): Change in moment coefficient due to angle of attack[-]

    """

    TheilslopesResults = stats.theilslopes(delta_e, alpha, alpha=0.99)
    slope = TheilslopesResults[0]
    print(slope)
    return -slope * Cmdelta
