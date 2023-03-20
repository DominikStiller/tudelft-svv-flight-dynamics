from fd.simulation import constants
import scipy.stats as stats
import math
import numpy as np
from fd import conversion as con
import matplotlib.pyplot as plt


def lin_moment_mass():
    """

    Returns (float, float): Slope of the linear moment-vs-mass function and the intersect.

    """

    mass = np.linspace(100, 4900, 49)
    mass = np.append(mass, 5008.0)
    moment = np.array(
        [
            298.16,
            591.18,
            879.08,
            1165.42,
            1448.4,
            1732.53,
            2014.8,
            2298.84,
            2581.92,
            2866.3,
            3150.18,
            3434.52,
            3718.52,
            4003.23,
            4287.76,
            4572.24,
            4856.56,
            5141.16,
            5425.64,
            5709.9,
            5994.04,
            6278.47,
            6562.82,
            6846.96,
            7131.0,
            7415.33,
            7699.6,
            7984.34,
            8269.06,
            8554.05,
            8839.04,
            9124.8,
            9410.62,
            9696.97,
            9983.4,
            10270.08,
            10556.84,
            10843.87,
            11131.0,
            11418.2,
            11705.5,
            11993.31,
            12281.18,
            12569.04,
            12856.86,
            13144.73,
            13432.48,
            13720.56,
            14008.46,
            14320.34,
        ]
    )

    massSI = con.lbs_to_kg(mass)
    momentSI = con.inchpound_to_kgm(moment)

    result = stats.theilslopes(momentSI, massSI, alpha=0.99)

    return result[0], result[1]


slope, intersect = lin_moment_mass()


def get_cg(
    W, massP1, massP2, masscoor, mass1L, mass1R, mass2L, mass2R, mass3L, mass3R, shift=False
):
    mtot = W / constants.g
    mfuel = mtot - sum(masscoor, mass1L, mass1R, mass2L, mass2R, mass3L, mass3R, massP1, massP2)
    slope, intersect = lin_moment_mass()
    momentfuel = slope * mfuel + intersect
    xcg = (
        momentfuel
        + (massP1 + massP2) * constants.xcgP
        + masscoor * constants.xcgcoor
        + (mass1R + mass1L) * constants.xcg1
        + (mass2R + mass2L) * constants.xcg2
        + (mass3R + mass3L) * constants.xcg3
        + constants.OEW * xcgOEW
    ) / mtot
    return
