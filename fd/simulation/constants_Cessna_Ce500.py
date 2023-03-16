import numpy as np

#Aircraft geometry:
S = 24.2 # wing area [m^2]
lh = 5.5 # tail length [m]
c = 2.022 # mean aerodynamic cord [m]
KY2 = 0.980
th0 = 0
muc = 102.7
V0 = 59.9
m = 4547.8
xcg = 0.3*c


#Stability derivatives:
CX0 = 0
CXu = -0.2199
CXa = 0.4653
CXadot = 0
CXq = 0
CXde = 0

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

b = 15
g = 9.80665

A_prim = 4*(muc**2)*(KY2)*(CZadot-2*muc)
B_prim = Cmadot*2*muc*(CZq+2*muc)-Cmq*2*muc*(CZadot-2*muc)-2*muc*(KY2)*(CXu*(CZadot-2*muc)-2*muc*CZa)
C_prim = Cma*2*muc*(CZq+2*muc)-Cmadot*(2*muc*CX0+CXu*(CZq+2*muc))+Cmq*(CXu*(CZadot-2*muc)-2*muc*CZa)+2*muc*(KY2)*(CXa*CZu-CZa*CXu)
D_prim = Cmu*(CXa*(CZq+2*muc)-CZ0*(CZadot-2*muc))-Cma*(2*muc*CX0+CXu*(CZq+2*muc))+Cmadot*(CX0*CXu-CZ0*CZu)+Cmq*(CXu*CZa-CZu*CXa)
E_prim = -Cmu*(CX0*CXa+CZ0*CZa)+Cma*(CX0*CXu+CZ0*CZu)
p = (A_prim, B_prim, C_prim, D_prim, E_prim)
print("analytic", np.roots(p)*V0/c)

