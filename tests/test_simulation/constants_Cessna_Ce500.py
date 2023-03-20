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


def test_Cessna_Ce500_symmetric():
    A_prim = 4 * (muc**2) * (KY2) * (CZadot - 2 * muc)
    B_prim = (
        Cmadot * 2 * muc * (CZq + 2 * muc)
        - Cmq * 2 * muc * (CZadot - 2 * muc)
        - 2 * muc * (KY2) * (CXu * (CZadot - 2 * muc) - 2 * muc * CZa)
    )
    C_prim = (
        Cma * 2 * muc * (CZq + 2 * muc)
        - Cmadot * (2 * muc * CX0 + CXu * (CZq + 2 * muc))
        + Cmq * (CXu * (CZadot - 2 * muc) - 2 * muc * CZa)
        + 2 * muc * (KY2) * (CXa * CZu - CZa * CXu)
    )
    D_prim = (
        Cmu * (CXa * (CZq + 2 * muc) - CZ0 * (CZadot - 2 * muc))
        - Cma * (2 * muc * CX0 + CXu * (CZq + 2 * muc))
        + Cmadot * (CX0 * CXu - CZ0 * CZu)
        + Cmq * (CXu * CZa - CZu * CXa)
    )
    E_prim = -Cmu * (CX0 * CXa + CZ0 * CZa) + Cma * (CX0 * CXu + CZ0 * CZu)
    p = (A_prim, B_prim, C_prim, D_prim, E_prim)
    print(np.polynomial.polynomial.polyroots(p) * V0 / c)


def test_Cessna_Ce500_asymmetric():
    A_prim = 16 * mub**3 * (KX2 * KZ2 - KXZ**2)
    B_prim = (
        -4
        * mub**2
        * (2 * CYb * (KX2 * KZ2 - KXZ**2) + Cnr * KX2 + Clp * KZ2 + (Clr + Cnp) * KXZ)
    )
    C_prim = (
        2
        * mub
        * (
            (CYb * Cnr - CYr * Cnb) * KX2
            + (CYb * Clp - Clb * CYp) * KZ2
            + ((CYb * Cnp - Cnb * CYp) + (CYb * Clr - Clb * CYr)) * KXZ
            + 4 * mub * Cnb * KX2
            + 4 * mub * Clb * KXZ
            + 0.5 * (Clp * Cnr - Cnp * Clr)
        )
    )
    D_prim = (
        -4 * mub * CL * (Clb * KZ2 + Cnb * KXZ)
        + 2 * mub * (Clb * Cnp - Cnb * Clp)
        + 0.5 * CYb * (Clr * Cnp - Cnr * Clp)
        + 0.5 * CYp * (Clb * Cnr - Cnb * Clr)
        + 0.5 * CYr * (Clp * Cnb - Cnp * Clb)
    )
    E_prim = CL * (Clb * Cnr - Cnb * Clr)
    p = (A_prim, B_prim, C_prim, D_prim, E_prim)
    print(np.polynomial.polynomial.polyroots(p))
