from fd.simulation import constants
import scipy.stats as stats
import math

def calc_CL(W: float, V: float, S = constants.S, rho = constants.rho0) -> float:
    """
    Calculate CL for a given combination of W, rho, V and S.
    Args:
        W (array_like): Weight [-]
        rho (float): Air density (set to sea level) [kg/m3]
        V (array_like): Velocity (calibrated one since we are using rho at sea level) [m/s]
        S (float): Surface area [m2]

    Returns:
        (array_like): CL [-]
    """

    return 2 * W / (rho * V**2 * S)

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
    result = stats.theilslopes(CL, alpha)
    CLalpha = result.slope
    CL_alpha_equals0 = result.intercept
    # low_slope = result.low_slope # Lower bound of the confidence interval on slope
    # high_slope = result.high_slope# Higher bound of the confidence interval on slope
    alpha_0 = -CL_alpha_equals0/CLalpha

    return CLalpha, CL_alpha_equals0, alpha_0 # low_slope, high_slope


def calc_CD(T: float, Vc: float) -> float:
    """
    This function calculates the drag coefficient CD[-]

    Args:
        T (array_like): Thrust[N].
        Vc (array_like): computed airspeed[ms-1].

    Returns:
        CD (array_like): Drag coefficient CD[-].

    """    
    return T/(0.5*constants.rho0*Vc**2*constants.S)

def calc_CD0_e(CD: list, CL: list) ->float:
    """
    This function uses the parabolic drag formula to calculate the zero lift drag, CD0[-], and the oswald
    efficiency factor, e[-].

    Args:
        CD (list): The drag coefficient CD[-].
        CL (list): The lift coefficient CL[-].

    Returns:
        CD0, e (float, float): Zero lift drag coefficient, CD0[-], oswald efficiency factor, e[-].

    """
    TheilslopesResults = stats.theilslopes(CL**2, CD)
    CD0 = TheilslopesResults.intercept
    e = 1/(math.pi*constants.A*TheilslopesResults.slope)
    
    return CD0, e