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
    fig, ax = plt.subplots(figsize=(12, 6))

    aa = np.linspace(0, max(alpha), 20)

    ax.plot(
        aa,
        Clalpha * (aa - alpha0),
        "r",
        label="Best fit ($C_{L_\\alpha}$ = "
        + f"{Clalpha:.3} 1/°"
        + ", $\\alpha_0$ = "
        + f"{alpha0:.3} °)",
    )
    ax.scatter(alpha, CL, marker="x", color="black", s=50, label="Data")

    ax.set_xlabel(r"$\alpha$")
    ax.set_ylabel("$C_L$")

    ax.legend()

    format_plot()
    save_plot("data/", "cl_alpha")
    plt.show()


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
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

    yy = np.linspace(0, max(CL), 20)

    ax1.plot(CD0 + yy / (math.pi * constants.A * e), yy, "r")
    ax1.scatter(CD, CL**2, marker="x", color="black", s=50)
    ax1.set_xlabel("$C_D$")
    ax1.set_ylabel("$C_L^2$")

    ax2.plot(
        CD0 + yy**2 / (math.pi * constants.A * e),
        yy,
        "r",
        label="Best fit ($C_{D_0}$ = " + f"{CD0:.3}" + ", $e$ = " + f"{e:.3})",
    )
    ax2.scatter(CD, CL, marker="x", color="black", s=50, label="Data")
    ax2.set_xlabel("$C_D$")
    ax2.set_ylabel("$C_L$")
    ax2.legend()

    format_plot()
    save_plot("data/", "cl_cd")
    plt.show()
