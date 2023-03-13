import math

def deg_to_rad(x):
    """Convert value from degree to radians"""
    return math.pi / 180 * x

def lbshr_to_kgs(x):
    """Convert value from lbs/hr to kg/s"""
    return 0.45359237/3600 * x
def psi_to_Pa(x):
    """Convert value from psi to Pa"""
    return 6894.757 * x

def ftmin_to_ms(x):
    """Convert value from ft/min to m/s"""
    return 0.3048/60 * x