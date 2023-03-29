import numpy as np

from fd.simulation.constants import *


def time_constant_aperiodic_roll(eig: complex, V0):
    """
    Calculating time constant for the aperiodic roll
    Args:
        eig: eigenvalue for the aperiodic roll

    Returns:
        time constant
    """
    if abs(eig.imag) > 0:
        print(f"WARNING: aperiodic roll eigenvalue should be real, is {eig}")
    tau = -(1 / eig.real)
    return tau


def time_constant_spiral(eig: complex):
    """
    Calculating time constant for the spiral
    Args:
        eig: eigenvalue for the spiral

    Returns:
        Time constant
    """
    if abs(eig.imag) > 0:
        print(f"WARNING: spiral eigenvalue should be real, is {eig}")
    tau = -(1 / eig.real)
    return tau


def characteristics_dutch_roll(eig: complex):
    """
    Calculating the period and time to damp to half amplitude for the Dutch roll
    Args:
        eig: eigenvalue for the Dutch roll

    Returns:
        Period, time to half amplitude
    """
    P = (2 * pi) / abs(eig.imag)
    T_half = np.log(0.5) / eig.real
    return P, T_half


def characteristics_phugoid(eig: complex):
    """
    Calculating the period and time to damp to half amplitude for the Dutch roll
    Args:
        eig: eigenvalue for the phugoid

    Returns:
        Period, time to half amplitude
    """
    P = (2 * pi) / abs(eig.imag)
    T_half = np.log(0.5) / eig.real
    return P, T_half


def characteristics_short_period(eig: complex):
    """
    Calculating the period and time to damp to half amplitude for the short period
    Args:
        eig: eigenvalue for the short period

    Returns:
        Period, time to half amplitude
    """
    P = (2 * pi) / abs(eig.imag)
    T_half = np.log(0.5) / eig.real
    return P, T_half
