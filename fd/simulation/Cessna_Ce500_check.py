from fd.structs import AerodynamicParameters
from fd.simulation.aircraft_model import AircraftModel
import control.matlab as ml
import numpy as np
from fd.simulation.constants_Cessna_Ce500 import *
import matplotlib.pyplot as plt

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
A, B, C, D = model.get_state_space_matrices_symmetric(m, V0, rho, th0)
#eigenval, eigenvec = model.get_eigenvalue_and_eigenvector(A)
#print(eigenval*b/V0)
sys = ml.ss(A, B, C, D)
#ml.damp(sys, doprint = True)

dt = 0.2
maneuvre_duration = 10 #seconds
t = np.arange(0, maneuvre_duration+dt, dt)
x_0 = np.array([0, 0, 0, 0])
input_duration = 2.5
input_value = -0.005 #rad
u = np.zeros(t.shape)
u[:int(input_duration/dt)] = input_value*np.ones(u[:int(input_duration/dt)].size)
plt.plot(t, u)
plt.show()
yout, t, xout = ml.lsim(sys, u, t, x_0)
fig, axs = plt.subplots(2, 2, sharex=True)

axs[0, 0].plot(t, xout[:, 0]+V0*np.ones(t.size))
axs[0, 0].set_title("Absolut deviation in airspeed")

axs[1, 0].plot(t, xout[:, 1])
axs[1, 0].set_title("Angle of attack")

axs[0, 1].plot(t, xout[:, 2])
axs[0, 1].set_title("Pitch angle")

axs[1, 1].plot(t, xout[:, 3])
axs[1, 1].set_title("Pitch rate q")

plt.show()

