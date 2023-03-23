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
            29816,
            59118,
            87908,
            116542,
            144840,
            173253,
            201480,
            229884,
            258192,
            286630,
            315018,
            343452,
            371852,
            400323,
            428776,
            457224,
            485656,
            514116,
            542564,
            570990,
            599404,
            627847,
            656282,
            684696,
            713100,
            741533,
            769960,
            798434,
            826906,
            855405,
            883904,
            912480,
            941062,
            969697,
            998340,
            1027008,
            1055684,
            1084387,
            1113100,
            1141820,
            1170550,
            1199331,
            1228118,
            1256904,
            1285686,
            1314473,
            1343248,
            1372056,
            1400846,
            1432034,
        ]
    )

    massSI = con.lbs_to_kg(mass)
    momentSI = con.inchpound_to_kgm(moment)

    result = stats.theilslopes(momentSI, massSI, alpha=0.99)
    plt.scatter(massSI, momentSI)
    plt.plot(massSI, result[0] * massSI + result[1])
    plt.show()

    return result[0], result[1]


slope, intersect = lin_moment_mass()


def get_cg(
    mfuel, massP1, massP2, masscoor, mass1L, mass1R, mass2L, mass2R, mass3L, mass3R, shift=False
):
    mtot = (
        mfuel
        + masscoor
        + mass1L
        + massP2
        + mass2L
        + mass3L
        + mass3R
        + mass1R
        + mass2R
        + massP1
        + constants.OEW
    )
    slope, intersect = lin_moment_mass()
    momentfuel = slope * mfuel + intersect
    if shift:
        xcg = (
            momentfuel
            + (massP1 + massP2) * constants.xcgP
            + (constants.xcgP + (constants.xcgcoor - constants.xcgP) * 2 / 3) * mass3R
            + masscoor * constants.xcgcoor
            + (mass1R + mass1L) * constants.xcg1
            + (mass2R + mass2L) * constants.xcg2
            + (mass3L) * constants.xcg3
            + constants.OEW * constants.xcgOEW
        ) / mtot
    else:
        xcg = (
            momentfuel
            + (massP1 + massP2) * constants.xcgP
            + masscoor * constants.xcgcoor
            + (mass1R + mass1L) * constants.xcg1
            + (mass2R + mass2L) * constants.xcg2
            + (mass3R + mass3L) * constants.xcg3
            + constants.OEW * constants.xcgOEW
        ) / mtot
    return xcg
