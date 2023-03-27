import numpy as np
from math import *
from fd.simulation.constants import *
from fd.simulation.aircraft_model import AircraftModel
from fd.structs import AerodynamicParameters


def time_constant_aperiodic_roll(eig, Ve):
    """ "
    Calculating time constant for the aperiodic roll
    Args:
        eig: eigenvalue for aperiodic roll
        Ve: equivalent velocity during aperiodic roll

    """
    tau = -(1 / eig) * (c / Ve)

    return tau


def time_constant_spiral(eig, Ve):
    """ "
    Calculating time constant for the spiral
    Args:
        eig: eigenvalue spiral
        Ve: equivalent velocity during spiral

    """
    tau = -(1 / eig) * (c / Ve)
    return tau


def characteristics_dutch_roll(imag_eig, real_eig, Ve):

    """
    Calculating the period and time to damp to half amplitude for the Dutch roll
    Args:
        imag_eig: imaginairy part of eigenvalue for dutch roll
        real_eig: real part of eigenvalue of dutch roll
        Ve: equivalent velocity during dutch roll

    Returns:

    """
    P = ((2 * pi) / (imag_eig)) * (b / Ve)
    T_half = (np.log(0.5) / (real_eig)) * (b / Ve)
    return P, T_half


def characteristics_phugoid(imag_eig, real_eig, Ve):

    """
    Calculating the period and time to damp to half amplitude for the Dutch roll
    Args:
        imag_eig: imaginairy part of eigenvalue for phugoid
        real_eig: real part of eigenvalue of phugoid
        Ve: equivalent velocity during phugoid

    Returns:

    """
    P = ((2 * pi) / (imag_eig)) * (b / Ve)
    T_half = (np.log(0.5) / (real_eig)) * (b / Ve)
    return P, T_half


def characteristics_short_period(imag_eig, real_eig, Ve):
    """
    Calculating the period and time to damp to half amplitude for the short period
    Args:
        imag_eig: imaginairy part of eigenvalue for the short period
        real_eig: real part of eigenvalue of the short period
        Ve: equivalent velocity during the short period

    Returns:

    """
    P = ((2 * pi) / (imag_eig)) * (b / Ve)
    T_half = (np.log(0.5) / (real_eig)) * (b / Ve)
    return P, T_half
