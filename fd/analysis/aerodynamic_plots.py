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
    aa = np.linspace(min(alpha), max(alpha), 10)
    plt.plot(aa, Clalpha * (aa - alpha0), "r")
    plt.scatter(alpha, CL, marker="x", color="black", s=50)
    plt.xlabel("α")
    plt.ylabel("C_L")

    format_plot()
    plt.show()


# plot_cl_alpha([0.1, 0.2, 0.3, 0.4, 0.5], [-1, 1, 3, 5, 7], 0.05, -3)


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


# Elevator curves
# [TODO: add fully controlled parameters to plots once these are know and its known how we want to do this]
def plot_delta_e__alpha(delta_e, alpha, xlabel_input=r"$\alpha$ [deg]", ylabel_input="$\delta_e$ [deg]"):
    """
    Plotting of the elevator trim curve: delta_e vs alpha. Should be linear.

    Args:
        delta_e (array_like): Elevator deflection [deg]
        alpha (array_like):  angle of attack [deg]

    Returns:
        delta-alpha elevator trim curve plot
    """
    delta_e_alpha_slope, delta_e_intercept, _, _ = stats.theilslopes(delta_e, alpha, alpha=0.99)
    xx = np.linspace(0, max(alpha)*1.05, 2)
    plt.plot(xx, delta_e_alpha_slope * xx + delta_e_intercept, "r")
    plt.scatter(alpha, delta_e, marker="x", color="black", s=50)
    plt.xlabel(xlabel_input)
    plt.ylabel(ylabel_input)

    format_plot()
    plt.show()

    delta_e_alpha_slope, delta_e_intercept, _, _ = stats.theilslopes(delta_e, alpha-(-1.3598111316603236), alpha=0.99)

    return delta_e_alpha_slope, delta_e_intercept

def plot_delta_e__Vinv(delta_e, Vsquaredinv, xlabel_input="$1/V^2$ [m/s]", ylabel_input="$\delta_e$ [deg]"):
    """
    Plotting of the elevator trim curve. Depending on if V or V_e is used and the deflection or reduced deflection,
    the non-reduced or reduced elevator trim curve plots can be obtained. Should be proportional to 1/v^2

    Args:
        delta_e (array_like): elevator deflection [deg]
        Vsquaredinv (array_like): airspeed squared and reciprocal [1/(m/s)^2]

    Returns:
        delta-V elevator trim curve plot
    """
    delta_e_v_slope, delta_e_intercept, _, _ = stats.theilslopes(delta_e, Vsquaredinv, alpha=0.99)
    xx = np.linspace(0, max(Vsquaredinv)*1.05, 2)
    plt.plot(xx, delta_e_v_slope * xx + delta_e_intercept, "r")
    plt.scatter(Vsquaredinv, delta_e, marker="x", color="black", s=50)
    plt.xlabel(xlabel_input)
    plt.ylabel(ylabel_input)

    format_plot()
    plt.show()

def plot_delta_e__Ve(delta_e, V, xlabel_input="$V$ [m/s]", ylabel_input="$\delta_e$ [deg]"):
    """
    Plotting of the elevator trim curve. Depending on if V or V_e is used and the deflection or reduced deflection,
    the non-reduced or reduced elevator trim curve plots can be obtained. Should be proportional to 1/v^2

    Args:
        delta_e (array_like): elevator deflection [deg]
        V (array_like): airspeed [m/s]

    Returns:
        delta-V elevator trim curve plot
    """
    delta_e_v_slope, delta_e_intercept, _, _ = stats.theilslopes(delta_e, 1/V**2, alpha=0.99)
    xx = np.linspace(106*0.514444, 200, 100)
    plt.plot(xx, delta_e_intercept + delta_e_v_slope * 1/ xx**2, "r")
    plt.scatter(V, delta_e, marker="x", color="black", s=50)
    plt.xlabel(xlabel_input)
    plt.ylabel(ylabel_input)

    format_plot()
    plt.show()

def plot_Fe_Ve(F_e, V, xlabel_input="V_e^*", ylabel_input="F_e^*"):
    """
    Plotting of the elevator control force curve. Depending on if V or V_e is used and the F_e or reduced F_e,
    the non-reduced or reduced elevator control force curve plots can be obtained.

    Args:
        F_e (array_like): Control stick force [N]
        V (array_like): airspeed [m/s]

    Returns:
        F_e-V elevator control force curve plot
    """
    plt.scatter(V, F_e, marker="x", color="black", s=50)
    plt.xlabel(xlabel_input)
    plt.ylabel(ylabel_input)

    format_plot()
    plt.show()


# Testing
import os
if os.getcwd().endswith("private"):
    os.chdir("..")

import sys
sys.path.append(".")

from fd.analysis.flight_test import FlightTest
from fd.plotting import format_plot

import matplotlib.pyplot as plt

ft = FlightTest("../../data/B24")
print(vars(ft.df_elevator_trim))
df = ft.df_elevator_trim
alpha = df["alpha"]
delta_e = df["delta_e"]
cas = df["cas"]
F_e = df["F_e"]

# a = [5.2, 6.3, 7.1, 4.3, 3.616666667]
# delta_e = [0.1, -0.35, -0.7, 0.5, 0.7]
# IAS = [154.6666667,143,135.8333333,165.6666667,173.6666667]
# d_trim = [3,3,3,3,3]
# Fe = [-1, -18, -30, 22, 38]

delta_e_alpha_slope, delta_e_intercept = plot_delta_e__alpha(delta_e,alpha)
plot_delta_e__Vinv(delta_e,1/cas **2)
plot_delta_e__Ve(delta_e, cas)
# plot_deltae_Ve(delta_e, alpha, xlabel_input="alpha")
# plot_deltae_Ve(delta_e, cas)
# plot_Fe_Ve(F_e, alpha, xlabel_input="alpha")
# plot_Fe_Ve(F_e, cas)
