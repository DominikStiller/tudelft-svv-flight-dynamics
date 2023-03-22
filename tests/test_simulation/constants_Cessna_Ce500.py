import numpy as np

# Aircraft geometry:
V = 59.9
S = 24.2  # wing area [m^2]
lh = 5.5  # tail length [m]
c = 2.022  # mean aerodynamic cord [m]
KY2 = 0.980
KX2 = 0.012
KZ2 = 0.037
KXZ = 0.002
th0 = 0
muc = 102.7
mub = 15.5
V0 = 59.9
m = 4547.8
xcg = 0.3 * c


# Stability derivatives:
CX0 = 0
CXu = -0.2199
CXa = 0.4653
CXadot = 0
CXq = 0
CXde = 0

CYb = -0.9896
CYp = -0.0870
CYr = 0.4300
CYda = 0
CYdr = 0.3037
CYbdot = 0

CZ0 = -1.1360
CZu = -2.2720
CZa = -5.1600
CZadot = -1.4300
CZq = -3.8600
CZde = -0.6238

Cmu = 0
Cma = -0.4300
Cmadot = -3.7000
Cmq = -7.0400
Cmde = -1.5530

Cnb = 0.1638
Cnp = -0.0108
Cnr = -0.1930
Cnda = 0.0286
Cndr = -0.1261
Cnbdot = 0


Clb = -0.0772
Clp = -0.3444
Clr = 0.2800
Clda = -0.2349
Cldr = 0.0286

b = 13.36
g = 9.80665
CL = 1.1360
