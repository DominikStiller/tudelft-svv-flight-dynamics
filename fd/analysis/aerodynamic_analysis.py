from fd.simulation import constants
import scipy.stats as stats
import math

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