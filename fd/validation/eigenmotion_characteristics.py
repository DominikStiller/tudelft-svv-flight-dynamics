import numpy as np

from fd.simulation.constants import *


def time_constant_aperiodic_roll(eig: complex, V):
    """
    Calculating time constant for the aperiodic roll
    Args:
        eig: eigenvalue for the aperiodic roll
        V: velocity during aperiodic roll

    Returns:
        time constant
    """
    if abs(eig.imag) > 0:
        print(f"WARNING: aperiodic roll eigenvalue should be real, is {eig}")
    tau = -(1 / eig.real) * (c / V)
    return tau


def time_constant_spiral(eig: complex, V):
    """
    Calculating time constant for the spiral
    Args:
        eig: eigenvalue for the spiral
        V: velocity during spiral

    Returns:
        Time constant
    """
    if abs(eig.imag) > 0:
        print(f"WARNING: spiral eigenvalue should be real, is {eig}")
    tau = -(1 / eig.real) * (c / V)
    return tau


def characteristics_dutch_roll(eig: complex, V):
    """
    Calculating the period and time to damp to half amplitude for the Dutch roll
    Args:
        eig: eigenvalue for the Dutch roll
        V: velocity during Dutch roll

    Returns:
        Period, time to half amplitude
    """
    P = ((2 * pi) / abs(eig.imag)) * (b / V)
    T_half = (np.log(0.5) / eig.real) * (b / V)
    return P, T_half


def characteristics_phugoid(eig: complex, V):
    """
    Calculating the period and time to damp to half amplitude for the Dutch roll
    Args:
        eig: eigenvalue for the phugoid
        V: velocity during phugoid

    Returns:
        Period, time to half amplitude
    """
    P = ((2 * pi) / abs(eig.imag)) * (b / V)
    T_half = (np.log(0.5) / eig.real) * (b / V)
    return P, T_half


def characteristics_short_period(eig: complex, V):
    """
    Calculating the period and time to damp to half amplitude for the short period
    Args:
        eig: eigenvalue for the short period
        V: velocity during the short period

    Returns:
        Period, time to half amplitude
    """
    P = ((2 * pi) / abs(eig.imag)) * (b / V)
    T_half = (np.log(0.5) / eig.real) * (b / V)
    return P, T_half
