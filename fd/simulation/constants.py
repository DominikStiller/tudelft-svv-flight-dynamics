# Citation 550 - Linear simulation
from fd import conversion as conv
from math import sin, cos, pi

# xcg = 0.25 * c

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
OEW = 4160.75745 # operational empty weight[kg]
xcgOEW = 7.410196 # centre of gravity for operational empty weight[m]
xcgP = 131*0.0254
xcgcoor = 170*0.0254
xcg1 = 214*0.0254
xcg2 = 251*0.0254
xcg3 = 288*0.0254

# Constant values concerning atmosphere and gravity
rho0 = 1.2250  # air density at sea level [kg/m^3]
p0 = 101325  # air pressure at sea level [Pa]
Tempgrad = -0.0065  # temperature gradient in ISA [K/m]
Temp0 = 288.15  # temperature at sea level in ISA [K]
R = 287.05  # specific gas constant [m^2/s^2K]
g = 9.81  # [m/s^2] (gravity constant)
gamma = 1.4  #

# air density [kg/m^3]
# rho = rho0 * (1 + (Tempgrad * hp0 / Temp0)) ** -((g / (Tempgrad * R)) + 1)
# W = m * g  # [N]       (aircraft weight)

# Constant values concerning aircraft inertia
# muc = m / (rho * S * c)
# mub = m / (rho * S * b)
KX2 = 0.019
KZ2 = 0.042
KXZ = 0.002
KY2 = 1.25 * 1.114

# Aerodynamic constants
Cmac = 0  # Moment coefficient about the aerodynamic centre [-]
# CNwa = CLa  # Wing normal force slope [-]
CNha = 2 * pi * Ah / (Ah + 2)  # Stabilizer normal force slope [-]
depsda = 4 / (A + 2)  # Downwash gradient [-]

# Lift and drag coefficient
# CL = 2 * W / (rho * V0**2 * S)  # Lift coefficient [-]
# CD = CD0 + (CLa * alpha0) ** 2 / (pi * A * e)  # Drag coefficient [-]

# standard values
Ws = 60500  # standard weight from the assignment

# Stability derivatives
# CX0 = W * sin(th0) / (0.5 * rho * V0**2 * S)
CXu = -0.09500
CXa = +0.47966  # Positive, see FD lecture notes
CXadot = +0.08330
CXq = -0.28170
CXde = -0.03728

# CZ0 = -W * cos(th0) / (0.5 * rho * V0**2 * S)
CZu = -0.37616
CZa = -5.74340
CZadot = -0.00350
CZq = -5.66290
CZde = -0.69612

Cm0 = +0.0297
Cmu = +0.06990
Cmadot = +0.17800
Cmq = -8.79415
CmTc = -0.0064

CYb = -0.7500
CYbdot = 0
CYp = -0.0304
CYr = +0.8495
CYda = -0.0400
CYdr = +0.2300

Clb = -0.10260
Clp = -0.71085
Clr = +0.23760
Clda = -0.23088
Cldr = +0.03440

Cnb = +0.1348
Cnbdot = 0
Cnp = -0.0602
Cnr = -0.2061
Cnda = -0.0120
Cndr = -0.0939
