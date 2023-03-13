import math


def deg_to_rad(deg):
    """Convert value from degree to radians"""
    return math.pi / 180 * deg


def lbshr_to_kgs(lbshr):
    """Convert value from lbs/hr to kg/s"""
    return 0.45359237 / 3600 * lbshr


def psi_to_Pa(x):
    """Convert value from psi to Pa"""
    return 6894.757 * psi


def ftmin_to_ms(ftmin):
    """Convert value from ft/min to m/s"""
    return 0.3048 / 60 * ftmin


def lbs_to_kg(lbs):
    """Convert mass in pounds to kilograms"""
    return 0.45359237 * lbs


def kts_to_ms(kts):
    """Convert speed in knots to meters per second"""
    return 1852 / 3600 * kts


def ft_to_m(ft):
    """Convert distance in feet to meters"""
    return 0.3048 * ft


def C_to_K(C):
    """Convert temperature in Celcius to Kelvin"""
    return 273.15 + C
