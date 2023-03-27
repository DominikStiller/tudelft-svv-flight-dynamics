from fd.analysis.flight_test import FlightTest
from fd.simulation.aircraft_model import AircraftModel
from fd.simulation.simulation import Simulation
from fd.structs import AerodynamicParameters
from fd.validation.comparison import SimulatedMeasuredComparison
import control.matlab as ml
import numpy as np
from fd.plotting import *


flight_test = FlightTest("data/B24")
aero_params = AerodynamicParameters(
    C_L_alpha=4.758556374647304,
    alpha_0=-0.02312478307006348,
    C_D_0=0.023439123324849084,
    C_m_alpha=-0.5554065208385275,
    C_m_delta=-1.3380975545274032,
    e=1.0713238368125688,
)

aircraft_model = AircraftModel(aero_params)
A, B, C, D = aircraft_model.get_state_space_matrices_symmetric(5000, 150, 0.6, 0)
eig = np.linalg.eig(A)[0]

x_sym = eig.real
y_sym = eig.imag

A, B, C, D = aircraft_model.get_state_space_matrices_asymmetric(5000, 150, 0.6, 0, 0.8)
eigassym = np.linalg.eig(A)[0]


x_sym = eig.real
y_sym = eig.imag

x_assym = eigassym.real
y_assym = eigassym.imag


plt.scatter(x_sym, y_sym, marker="x")
plt.ylabel("Imaginary part")
plt.xlabel("Real part")
format_plot()
save_plot("data/", "eig_symmetric")
plt.show()

plt.scatter(x_assym, y_assym, marker="x")
plt.ylabel("Imaginary part")
plt.xlabel("Real part")
format_plot()
save_plot("data/", "eig_asymmetric")
plt.show()
