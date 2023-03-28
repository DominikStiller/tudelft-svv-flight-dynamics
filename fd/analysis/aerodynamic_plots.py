import math
from functools import partial

import numpy as np
import scipy.stats as stats

from fd.analysis.aerodynamics import calc_dynamic_pressure
from fd.plotting import *
from fd.simulation import constants
from fd.simulation.constants import rho0


def plot_cl_alpha(CL, alpha, Clalpha, alpha0):
    """
    Plotting og the lift slope
    Args:
        CL (array_like): Lift coefficient CL[-]
        alpha (array_like): Angle of attack alpha[rad]
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
def plot_delta_e_V_inv_squared(
    delta_e, V_inv_squared, xlabel_input="$1/V^2$ [1/(m/s)$^2$]", ylabel_input="$\delta_e$ [rad]"
):
    """
    Plotting of the elevator trim curve for delta_e vs 1/v^2. Should be proportional to 1/v^2, so linear.

    Args:
        delta_e (array_like): elevator deflection [rad]
        V_inv_squared (array_like): reciprocal of airspeed squared [1/(m/s)^2]

    Returns:
        delta vs 1/V^2 elevator trim curve plot
    """
    fig, ax = plt.subplots(figsize=(12, 6))

    # slope, y_intercept, _, _ = stats.theilslopes(delta_e, V_inv_squared, alpha=0.99)
    slope, y_intercept, _, _, _ = stats.linregress(V_inv_squared, delta_e)
    xx = np.linspace(0, max(V_inv_squared) * 1.05, 2)
    plt.plot(xx, slope * xx + y_intercept, "r")
    plt.scatter(V_inv_squared, delta_e, marker="x", color="black", s=50)
    plt.xlabel(xlabel_input)
    plt.ylabel(ylabel_input)

    plt.gca().invert_yaxis()
    format_plot()
    plt.show()

    return slope, y_intercept


# [TODO: add fully controlled parameters to plots once these are know and its known how we want to do this]
def plot_elevator_trim_curve(
    delta_e,
    alpha,
    cas,
    C_m_0,
    C_m_delta_e,
    cas_stall,
):
    """
    Plotting of the elevator trim curve. Depending on if V or V_e is used and the deflection or reduced deflection,
    the non-reduced or reduced elevator trim curve plots can be obtained. Should be proportional to 1/v^2

    Args:
        delta_e (array_like): elevator deflection [rad]
        alpha (array_like):  angle of attack [rad]
        cas (array_like): CAS [m/s]
        cas_stall (float): stall CAS [m/s]
    """
    # slope_V_inv_sq, y_intercept_V_inv_sq, _, _ = stats.theilslopes(
    #     delta_e, 1 / cas**2, alpha=0.99
    # )
    slope_V_inv_sq, y_intercept_V_inv_sq, _, _, _ = stats.linregress(1 / cas**2, delta_e)
    delta_e_asymptote = -C_m_0 / C_m_delta_e

    fig, (ax_alpha, ax_V) = plt.subplots(1, 2, figsize=(12, 6), sharey="all")

    # Delta_e vs alpha
    # slope_alpha, y_intercept_alpha, _, _ = stats.theilslopes(delta_e, alpha, alpha=0.99)
    slope_alpha, y_intercept_alpha, _, _, _ = stats.linregress(alpha, delta_e)
    xx_alpha = np.linspace(0, max(alpha) * 1.05, 2)
    ax_alpha.plot(xx_alpha, slope_alpha * xx_alpha + y_intercept_alpha, "r", label=r"$\alpha$")
    ax_alpha.scatter(alpha, delta_e, marker="x", color="black", s=50)

    # Delta_e vs alpha - alpha0
    alpha0 = (y_intercept_V_inv_sq - y_intercept_alpha) / slope_alpha
    # slope_alpha0, y_intercept_alpha0, _, _ = stats.theilslopes(delta_e, alpha - alpha0, alpha=0.99)
    slope_alpha0, y_intercept_alpha0, _, _, _ = stats.linregress(alpha - alpha0, delta_e)
    xx_alpha0 = np.linspace(0, max(alpha - alpha0) * 1.05, 2)
    ax_alpha.plot(
        xx_alpha0, slope_alpha0 * xx_alpha0 + y_intercept_alpha0, "b", label=r"$\alpha-\alpha_0$"
    )
    ax_alpha.scatter(alpha - alpha0, delta_e, marker="x", color="black", s=50)
    ax_alpha.set_xlabel(r"$\alpha$ [rad]")
    ax_alpha.set_ylabel(r"$\delta_e^*$ [rad]")
    ax_alpha.legend()

    # Delta_e vs V
    xx_V = np.linspace(0.95 * cas_stall, 200, 100)
    ax_V.plot(xx_V, y_intercept_V_inv_sq + slope_V_inv_sq * 1 / xx_V**2, "r")
    ax_V.scatter(cas, delta_e, marker="x", color="black", s=50)
    ax_V.axvline(x=cas_stall, color="b", linestyle="--", label="$V_{stall}$")
    ax_V.axhline(
        y=delta_e_asymptote,
        linestyle=":",
        color="black",
        label=r"$\delta_{e,V \to \infty}}$ = " + f"{delta_e_asymptote:.3} rad",
    )
    ax_V.set_xlabel(r"$\tilde V_e$ [m/s]")
    ax_V.legend()
    ax_V.invert_yaxis()

    format_plot()
    save_plot("data/", "elevator_trim_curve")
    plt.show()


def plot_elevator_control_force(F_e, cas, cas_stall):
    """
    Plotting of the elevator control force curve. Depending on if V or V_e is used and the F_e or reduced F_e,
    the non-reduced or reduced elevator control force curve plots can be obtained.

    The elevator control force is plotted against dynamic pressure and airspeed.

    Args:
        F_e (array_like): Control stick force [N]
        cas (array_like): CAS [m/s]
        cas_stall (float): stall CAS [m/s]
    """
    fig, (ax_q, ax_V) = plt.subplots(1, 2, figsize=(12, 6), sharey="all")

    dynamic_pressure_stall = calc_dynamic_pressure(cas_stall, rho0)
    dynamic_pressure = partial(calc_dynamic_pressure, rho=rho0)(cas)
    # slope, y_intercept, _, _ = stats.theilslopes(F_e, dynamic_pressure, alpha=0.99)
    slope, y_intercept, _, _, _ = stats.linregress(dynamic_pressure, F_e)

    # F_e vs q
    xx_q = np.linspace(dynamic_pressure_stall, max(dynamic_pressure) * 1.05, 2)
    ax_q.plot(xx_q, y_intercept + slope * xx_q, "r")
    ax_q.scatter(dynamic_pressure, F_e, marker="x", color="black", s=50)
    ax_q.axvline(
        x=dynamic_pressure_stall,
        color="black",
        linestyle="--",
        label=r"$\frac{1}{2}\rho_0 V_{stall}^2$",
    )

    xx2 = np.linspace(0, dynamic_pressure_stall, 2)
    ax_q.plot(xx2, y_intercept + slope * xx2, "r", linestyle=":")
    ax_q.set_xlabel(r"$\frac{1}{2}\rho_0 \tilde V_e^2$ [Pa]")
    ax_q.set_ylabel(r"$F_e^*$ [N]")
    ax_q.legend()
    ax_q.invert_yaxis()

    # F_e vs V
    xx_V = np.linspace(cas_stall, 120, 100)
    ax_V.plot(
        xx_V,
        y_intercept + slope * partial(calc_dynamic_pressure, rho=rho0)(xx_V),
        "r",
        label="Best fit",
    )
    ax_V.scatter(cas, F_e, marker="x", color="black", s=50, label="Data")
    ax_V.axvline(x=cas_stall, color="black", linestyle="--", label="$V_{stall}$")

    xx_V_stall = np.linspace(0, cas_stall, 100)
    ax_V.plot(
        xx_V_stall,
        y_intercept + slope * partial(calc_dynamic_pressure, rho=rho0)(xx_V_stall),
        "r",
        linestyle=":",
    )
    ax_V.set_xlabel(r"$\tilde V_e$ [m/s]")
    ax_V.legend()
    ax_V.invert_yaxis()

    format_plot(zeroline=True)
    save_plot("data/", "elevator_force_curve")
    plt.show()
