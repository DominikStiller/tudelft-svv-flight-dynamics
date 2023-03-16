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
CL = 1.1360
b = 13.36
model = AircraftModel(aero_params)
A, B, C, D = model.get_state_space_matrices_asymmetric(m, V0, rho, th0, CL)
eigenval, eigenvec = model.get_eigenvalue_and_eigenvector(A)
print(eigenval*b/V0)
sys = ml.ss(A, B, C, D)
ml.damp(sys, doprint = True)

