import math

import numpy as np

from fd.plotting import *
from fd.simulation import constants
import scipy.stats as stats


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

    aa = np.linspace(alpha0, max(alpha), 20)

    ax.plot(
        aa,
        Clalpha * (aa - alpha0),
        "r",
        label="Best fit ($C_{L_\\alpha}$ = "
        + f"{Clalpha:.3} 1/rad"
        + ", $\\alpha_0$ = "
        + f"{alpha0:.3} rad)",
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



# Elevator curves
# [TODO: add fully controlled parameters to plots once these are know and its known how we want to do this]
def plot_delta_e_alpha(
    delta_e,
    alpha,
    delta_e_intercept_V_inv_squared,
    xlabel_input=r"$\alpha$ [deg]",
    ylabel_input="$\delta_e$ [deg]",
):
    """
    Plotting of the elevator trim curve: delta_e vs alpha. Should be linear.

    Args:
        delta_e (array_like): Elevator deflection [deg]
        alpha (array_like):  angle of attack [deg]

    Returns:
        delta-alpha elevator trim curve plot
    """
    # Delta_e vs alpha
    slope_alpha, y_intercept_alpha, _, _ = stats.theilslopes(delta_e, alpha, alpha=0.99)
    xx = np.linspace(0, max(alpha) * 1.05, 2)
    plt.plot(xx, slope_alpha * xx + y_intercept_alpha, "r", label=r"$\alpha$")
    plt.scatter(alpha, delta_e, marker="x", color="black", s=50)

    # Delta_e vs alpha - alpha0
    alpha0 = (delta_e_intercept_V_inv_squared - y_intercept_alpha) / slope_alpha
    slope_alpha0, y_intercept_alpha0, _, _ = stats.theilslopes(delta_e, alpha - alpha0, alpha=0.99)
    xx = np.linspace(0, max(alpha - alpha0) * 1.05, 2)
    plt.plot(xx, slope_alpha0 * xx + y_intercept_alpha0, "b", label=r"$\alpha-\alpha_0$")
    plt.scatter(alpha - alpha0, delta_e, marker="x", color="black", s=50)

    plt.legend()
    plt.xlabel(xlabel_input)
    plt.ylabel(ylabel_input)

    plt.gca().invert_yaxis()
    format_plot()
    plt.show()

    return slope_alpha, y_intercept_alpha, slope_alpha0, y_intercept_alpha0, alpha0


def plot_delta_e_V_inv_squared(
    delta_e, V_inv_squared, xlabel_input="$1/V^2$ [1/(m/s)$^2$]", ylabel_input="$\delta_e$ [deg]"
):
    """
    Plotting of the elevator trim curve for delta_e vs 1/v^2. Should be proportional to 1/v^2, so linear.

    Args:
        delta_e (array_like): elevator deflection [deg]
        V_inv_squared (array_like): reciprocal of airspeed squared [1/(m/s)^2]

    Returns:
        delta vs 1/V^2 elevator trim curve plot
    """
    slope, y_intercept, _, _ = stats.theilslopes(delta_e, V_inv_squared, alpha=0.99)
    xx = np.linspace(0, max(V_inv_squared) * 1.05, 2)
    plt.plot(xx, slope * xx + y_intercept, "r")
    plt.scatter(V_inv_squared, delta_e, marker="x", color="black", s=50)
    plt.xlabel(xlabel_input)
    plt.ylabel(ylabel_input)

    plt.gca().invert_yaxis()
    format_plot()
    plt.show()

    return slope, y_intercept


def plot_delta_e_V(
    delta_e,
    V,
    delta_e_asymptote,
    V_stall=106,
    xlabel_input="$V$ [m/s]",
    ylabel_input="$\delta_e$ [deg]",
):
    """
    Plotting of the elevator trim curve. Depending on if V or V_e is used and the deflection or reduced deflection,
    the non-reduced or reduced elevator trim curve plots can be obtained. Should be proportional to 1/v^2

    Args:
        delta_e (array_like): elevator deflection [deg]
        V (array_like): airspeed [m/s]
        V_stall (float): Stall speed in [kts]

    Returns:
        delta-V elevator trim curve plot
    """
    slope, y_intercept, _, _ = stats.theilslopes(delta_e, 1 / V**2, alpha=0.99)
    xx = np.linspace(0.95 * (V_stall * 0.514444), 200, 100)
    plt.plot(xx, y_intercept + slope * 1 / xx**2, "r")
    plt.scatter(V, delta_e, marker="x", color="black", s=50)
    plt.axvline(x=(V_stall * 0.514444), color="b", linestyle="--", label="$V_{stall}$")
    plt.axhline(y=delta_e_asymptote, linestyle=":", color="black", label="$\delta_{e,asymptote}}$")
    plt.legend()
    plt.xlabel(xlabel_input)
    plt.ylabel(ylabel_input)

    plt.gca().invert_yaxis()
    format_plot()
    plt.show()


def plot_Fe_dyn_p(
    F_e,
    dynamic_pressure,
    dynamic_pressure_stall,
    xlabel_input=r"$\frac{1}{2}\rho_0 V^2$ [Pa]",
    ylabel_input=r"$F_e$ [N]",
):
    """
    Plotting of the elevator control force curve. Plotting of the force vs dynamic pressure, should be linear.

    Args:
        F_e (array_like): Control stick force [N]
        dynamic_pressure (array_like): dynamic pressure [Pa]

    Returns:
        F_e-1/2rhoV^2 elevator control force curve
    """
    slope, y_intercept, _, _ = stats.theilslopes(F_e, dynamic_pressure, alpha=0.99)
    xx = np.linspace(dynamic_pressure_stall, max(dynamic_pressure) * 1.05, 2)
    plt.plot(xx, y_intercept + slope * xx, "r")
    plt.scatter(dynamic_pressure, F_e, marker="x", color="black", s=50)
    plt.axvline(
        x=(dynamic_pressure_stall),
        color="b",
        linestyle="--",
        label=r"$\frac{1}{2}\rho_0 V_{stall}^2$",
    )
    xx2 = np.linspace(0, dynamic_pressure_stall, 2)
    plt.plot(xx2, y_intercept + slope * xx2, "r", linestyle=":")
    plt.legend()
    plt.xlabel(xlabel_input)
    plt.ylabel(ylabel_input)

    plt.gca().invert_yaxis()
    format_plot()
    plt.show()


def plot_Fe_V(F_e, V, V_stall=106, xlabel_input="$V$ [m/s]", ylabel_input="$F_e$ [N]"):
    """
    Plotting of the elevator control force curve. Depending on if V or V_e is used and the F_e or reduced F_e,
    the non-reduced or reduced elevator control force curve plots can be obtained.

    Args:
        F_e (array_like): Control stick force [N]
        V (array_like): airspeed [m/s]
        V_stall (float): Stall speed in [kts]

    Returns:
        F_e-V elevator control force curve plot
    """
    slope, y_intercept, _, _ = stats.theilslopes(F_e, 0.5 * constants.rho0 * V**2, alpha=0.99)
    xx = np.linspace((V_stall * 0.514444), 120, 100)
    plt.plot(xx, y_intercept + slope * 0.5 * constants.rho0 * xx**2, "r")
    plt.scatter(V, F_e, marker="x", color="black", s=50)
    plt.axvline(x=(V_stall * 0.514444), color="b", linestyle="--", label="$V_{stall}$")
    xx2 = np.linspace(0, (V_stall * 0.514444), 100)
    plt.plot(xx2, y_intercept + slope * 0.5 * constants.rho0 * xx2**2, "r", linestyle=":")
    plt.legend()
    plt.xlabel(xlabel_input)
    plt.ylabel(ylabel_input)

    plt.gca().invert_yaxis()
    format_plot()
    plt.show()
