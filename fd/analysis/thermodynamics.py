import numpy as np

from fd.simulation import constants


def calc_static_pressure(hp):
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
                / calc_static_pressure(hp)
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


def calc_static_temperature(Ttot, M):
    """

    Args:
        Ttot (float): Total temperature[K]
        M (float): Mach number[-]

    Returns (float): Static temperature[K]

    """
    return Ttot / (1 + (constants.gamma - 1) / 2 * M**2)


def calc_density(p, T):
    """

    Args:
        p (float): Static pressure[Pa]
        T (float): Static temperature[K]

    Returns (float): Density[kg/m^3]

    """
    return p / (constants.R * T)
