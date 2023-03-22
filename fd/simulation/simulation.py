from fd.simulation.aircraft_model import AircraftModel
from fd.structs import SimulationOutput
import control.matlab as ml
import numpy as np
import os
if os.getcwd().endswith("private"):
    os.chdir("..")
from fd.analysis.flight_test import FlightTest
from fd.plotting import format_plot
from fd.simulation.aircraft_model import AircraftModel
from fd.structs import AerodynamicParameters

import matplotlib.pyplot as plt
import sys
sys.path.append(".")


class Simulation:
    def __init__(self, model: AircraftModel):
        self.model = model

    def simulate_dutch_roll(self, df_dutch_roll) -> SimulationOutput:


        input = df_dutch_roll["delta_e"]
        V0 = df_dutch_roll['tas'][0]
        alpha0 = df_dutch_roll["alpha"][0]
        theta0 = df_dutch_roll['theta'][0]
        q0 = df_dutch_roll['q'][0]
        m = (df_dutch_roll['m'][0]+df_dutch_roll['m'][-1])/2
        state0 = np.transpose(0, alpha0, theta0, q0)
        rho0 = df_dutch_roll['rho'][0]


        model = AircraftModel(AerodynamicParameters)
        A, B, C, D = model.get_state_space_matrices_symmetric(m, V0, rho0, theta0)
        sys = ml.ss(A, B, C, D)

if __name__ == '__main__':
    sim = Simulation(AircraftModel(AerodynamicParameters(aa)))
    df_dutch_roll = FlightTest("data/B24").df_dutch_roll
    sim.simulate_dutch_roll(df_dutch_roll)




