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
print(eigenval)
sys = ml.ss(A, B, C, D)
ml.damp(sys, doprint = True)

