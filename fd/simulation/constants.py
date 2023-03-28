# Citation 550 - Linear simulation

from math import pi

from fd.conversion import lbs_to_kg, in_to_m, kts_to_ms

g = 9.81  # [m/s^2] (gravity constant)

# Aircraft mass
mass_basic_empty = lbs_to_kg(9172.9)  # basic empty weight [kg]

# CG positions of components
xcgOEW = in_to_m(291.74)
xcgP = in_to_m(131)
xcgcoor = in_to_m(170)
xcg1 = in_to_m(214)
xcg2 = in_to_m(251)
xcg3 = in_to_m(288)

# Aircraft geometry
S = 30.00  # wing area [m^2]
Sh = 0.2 * S  # stabilizer area [m^2]
Sh_S = Sh / S  # [-]
lh = 0.71 * 5.968  # tail length [m]
c = 2.0569  # mean aerodynamic cord [m]
lh_c = lh / c  # [-]
b = 15.911  # wing span [m]
bh = 5.791  # stabilizer span [m]
A = b**2 / S  # wing aspect ratio [-]
Ah = bh**2 / Sh  # stabilizer aspect ratio [-]
Vh_V = 1  # [-]
ih = -2 * pi / 180  # stabilizer angle of incidence [rad]

# Constant values concerning atmosphere and gravity
rho0 = 1.2250  # air density at sea level [kg/m^3]
p0 = 101325  # air pressure at sea level [Pa]
Tempgrad = -0.0065  # temperature gradient in ISA [K/m]
Temp0 = 288.15  # temperature at sea level in ISA [K]
R = 287.05  # specific gas constant [m^2/s^2K]
gamma = 1.4  #
cas_stall = kts_to_ms(106)  # equivalent stall speed [m/s]

# Constant values concerning aircraft inertia
KX2 = 0.019
KZ2 = 0.042
KXZ = 0.002
KY2 = 1.25 * 1.114

# Aerodynamic constants
Cmac = 0  # Moment coefficient about the aerodynamic centre [-]
CNha = 2 * pi * Ah / (Ah + 2)  # Stabilizer normal force slope [-]
depsda = 4 / (A + 2)  # Downwash gradient [-]

# standard values
Ws = 60500  # standard weight from the assignment
fuel_flow_standard = 0.048  # [kg/s]

# Stability derivatives
# CX0 = W * sin(th0) / (0.5 * rho * V0**2 * S)
# CXu = -0.09500
CXu = -0.15
CXa = +0.47966  # Positive, see FD lecture notes
CXadot = +0.08330
CXq = -0.28170
CXde = -0.03728

# CZu = -0.37616
CZu = -0.45
# CZa = -5.74340
# CZa = -5.5
CZa = -6
# CZadot = -0.00350
CZadot = -0.005
CZq = -5.66290
CZde = -0.69612

Cm0 = +0.0297
# Cmu = +0.06990
Cmu = 0.1
Cmadot = +0.17800
Cmq = -8.79415
CmTc = -0.0064

CYb = -0.7500
CYbdot = 0
CYp = -0.0304
CYr = +0.8495
CYda = -0.0400
CYdr = +0.2300

# Clb = -0.10260
Clb = -0.08
# Clp = -0.71085
Clp = -0.6
# Clr = +0.23760
Clr = 0.19
# Clda = -0.23088
Clda = -0.3
# Cldr = +0.03440
Cldr = 0.035

# Cnb = +0.1348
Cnb = 0.12
Cnbdot = 0
# Cnp = -0.0602
Cnp = -0.1
# Cnr = -0.2061
Cnr = -0.28
# Cnda = -0.0120
Cnda = -0.03
Cndr = -0.0939

# Durations of the eigenmotions
# Used for data extraction and simulation
duration_phugoid = 120  # [s]
duration_short_period = 8  # [s]
duration_dutch_roll = 20  # [s]
duration_dutch_roll_yd = 10  # [s]
duration_aperiodic_roll = 12  # [s]
duration_spiral = 120  # [s]

# Lead times for eigenmotions w.r.t. timestamp
lead_phugoid = 1  # [s]
lead_short_period = 1  # [s]
lead_dutch_roll = 2  # [s]
lead_dutch_roll_yd = 3  # [s]
lead_aperiodic_roll = 1  # [s]
lead_spiral = 5  # [s]

print(Clb * Cnr - Cnb * Clr)
