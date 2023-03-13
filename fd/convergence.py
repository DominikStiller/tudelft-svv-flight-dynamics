import math

<<<<<<< Updated upstream
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
=======
def lbs_to_kg(lbs):
    """Convert mass in pounds to kilograms"""
    kg = 0.45359237*lbs
    return kg
    
def kts_to_ms(kts):
    """Convert speed in knots to meters per second"""
    ms = 1852/3600*kts
    return ms

def ft_to_m(ft):
    """Convert distance in feet to meters"""
    m = 0.3048*ft
    return m

def C_to_K(C):
    """Convert temperature in Celcius to Kelvin"""
    K = 273.15 + C
    return K
>>>>>>> Stashed changes
