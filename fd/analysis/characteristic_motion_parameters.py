import numpy as np

from fd.simulation.constants import *


def time_constant_aperiodic_roll(eig: complex, Ve):
    """
    Calculating time constant for the aperiodic roll
    Args:
        eig: eigenvalue for the aperiodic roll
        Ve: equivalent velocity during aperiodic roll

    Returns:
        time constant
    """
    tau = -(1 / eig.real) * (c / Ve)
    return tau


def time_constant_spiral(eig: complex, Ve):
    """
    Calculating time constant for the spiral
    Args:
        eig: eigenvalue for the spiral
        Ve: equivalent velocity during spiral

    Returns:
        Time constant
    """
    tau = -(1 / eig.real) * (c / Ve)
    return tau


def characteristics_dutch_roll(eig: complex, Ve):
    """
    Calculating the period and time to damp to half amplitude for the Dutch roll
    Args:
        eig: eigenvalue for the Dutch roll
        Ve: equivalent velocity during Dutch roll

    Returns:
        Period, time to half amplitude
    """
    P = ((2 * pi) / eig.imag) * (b / Ve)
    T_half = (np.log(0.5) / eig.real) * (b / Ve)
    return P, T_half


def characteristics_phugoid(eig: complex, Ve):
    """
    Calculating the period and time to damp to half amplitude for the Dutch roll
    Args:
        eig: eigenvalue for the phugoid
        Ve: equivalent velocity during phugoid

    Returns:
        Period, time to half amplitude
    """
    P = ((2 * pi) / eig.imag) * (b / Ve)
    T_half = (np.log(0.5) / eig.real) * (b / Ve)
    return P, T_half


def characteristics_short_period(eig: complex, Ve):
    """
    Calculating the period and time to damp to half amplitude for the short period
    Args:
        eig: eigenvalue for the short period
        Ve: equivalent velocity during the short period

    Returns:
        Period, time to half amplitude
    """
    P = ((2 * pi) / eig.imag) * (b / Ve)
    T_half = (np.log(0.5) / eig.real) * (b / Ve)
    return P, T_half
