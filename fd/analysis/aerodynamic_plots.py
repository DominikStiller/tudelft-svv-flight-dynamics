import math

import numpy as np

from fd.plotting import *
from fd.simulation import constants


def plot_cl_alpha(CL, alpha, Clalpha, alpha0):
    """
    Plotting og the lift slope
    Args:
        CL (array_like): Lift coefficient CL[-]
        alpha (array_like): Angle of attack alpha[deg]
        Clalpha (float): slope of the linear part of the CL-alpha curve
        alpha0 (float): crossing of the CL-alpha curve with the alpha axis

    Returns:

    """
    aa = np.linspace(min(alpha), max(alpha), 10)
    plt.plot(aa, Clalpha * (aa - alpha0), "r")
    plt.scatter(alpha, CL, marker="x", color="black", s=50)
    plt.xlabel("α")
    plt.ylabel("C_L")

    format_plot()
    plt.show()


plot_cl_alpha([0.1, 0.2, 0.3, 0.4, 0.5], [-1, 1, 3, 5, 7], 0.05, -3)


def plot_cl_cd(CL, CD, CD0, e):
    """
    Plotting of the drag polar
    Args:
        CL (array_like): The lift coefficient[-]
        CD (array_like): The drag coefficient[-]
        CD0 (float): Zero lift drag coeffcient[-]
        e (float): The oswald efficiency factor[-]

    Returns:

    """
    yy = np.linspace(min(CL), max(CL), 10)
    plt.scatter(CD, CL)
    plt.plot(CD0 + yy**2 / (math.pi * constants.A * e), yy)

    format_plot()
    plt.show()


# plot_cl_cd([0.5, 0.84, 0.29, 0.955], [0.031, 0.0532, 0.0240, 0.0630], 0.02, 0.8)
