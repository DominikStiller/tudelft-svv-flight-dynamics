from B24.fd.structs import AerodynamicParameters
from B24.fd.simulation.aircraft_model import AircraftModel
import control.matlab as ml
import numpy as np
from B24.fd.simulation.constants_Cessna_Ce500 import *

aero_params = AerodynamicParameters
aero_params.C_m_alpha = -0.4300
aero_params.C_m_delta = -1.5530
m = 4547.8
V0 = 59.9
rho = 0.904627056
th0 = 0
model = AircraftModel(aero_params)
A, B, C, D = model.get_state_space_matrices_symmetric(m, V0, rho, th0)
eigenval, eigenvec = model.get_eigenvalue_and_eigenvector(A)
sys = ml.ss(A, B, C, D)
ml.damp(sys, doprint = True)

A_prim = 4*(muc**2)*(KY2)*(CZadot-2*muc)
B_prim = Cmadot*2*muc*(CZq+2*muc)-Cmq*2*muc*(CZadot-2*muc)-2*muc*(KY2)*(CXu*(CZadot-2*muc)-2*muc*CZa)
C_prim = Cma*2*muc*(CZq+2*muc)-Cmadot*(2*muc*CX0+CXu*(CZq+2*muc))+Cmq*(CXu*(CZadot-2*muc)-2*muc*CZa)+2*muc*(KY2)*(CXa*CZu-CZa*CXu)
D_prim = Cmu*(CXa*(CZq+2*muc)-CZ0*(CZadot-2*muc))-Cma*(2*muc*CX0+CXu*(CZq+2*muc))+Cmadot*(CX0*CXu-CZ0*CZu)+Cmq*(CXu*CZa-CZu*CXa)
E_prim = -Cmu*(CX0*CXa+CZ0*CZa)+Cma*(CX0*CXu+CZ0*CZu)
p = (A_prim, B_prim, C_prim, D_prim, E_prim)
print(np.roots(p))