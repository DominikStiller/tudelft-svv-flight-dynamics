from fd.simulation.aircraft_model import AircraftModel
from fd.structs import SimulationOutput
import control.matlab as ml
import pandas as pd
import numpy as np
from fd.analysis.flight_test import FlightTest
from fd.plotting import format_plot
from fd.simulation.aircraft_model import AircraftModel
from fd.structs import AerodynamicParameters
from fd.analysis.aerodynamics import calc_CL
from fd.simulation import constants


import matplotlib.pyplot as plt


class Simulation:
    def __init__(self, model: AircraftModel):
        self.model = model

    def simulate_dutch_roll(self, df_dutch_roll) -> SimulationOutput:

        input = np.transpose(np.hstack((df_dutch_roll["delta_a"], df_dutch_roll['delta_r'])))
        t = np.arange(0, len(input)*0.1, 0.1)
        V0 = df_dutch_roll['tas'][0]
        beta0 = 0
        phi0 = df_dutch_roll['phi'][0]
        theta0 = df_dutch_roll['theta'][0]
        p0 = df_dutch_roll['p'][0]
        r0 = df_dutch_roll['r'][0]
        m = (df_dutch_roll['m'][0]+df_dutch_roll['m'][-1])/2
        rho0 = df_dutch_roll['rho'][0]
        CL = calc_CL(m*constants.g*np.cos(theta0), V0, rho0)
        state0 = np.transpose(beta0, phi0, p0, r0)


        model = AircraftModel(AerodynamicParameters)
        A, B, C, D = model.get_state_space_matrices_asymmetric(m, V0, rho0, theta0, CL)
        sys = ml.ss(A, B, C, D)
        yout, t, xout = ml.lsim(sys, input, t, state0)
        result = np.hstack((t, np.transpose(yout)))
        df_result = pd.DataFrame(result, columns = ['t', 'beta', 'phi', 'p', 'r'])
        return df_result

    def simulate_dutch_roll_yd(self, df_dutch_roll_yd):
        input = np.transpose(np.hstack((df_dutch_roll_yd["delta_a"], df_dutch_roll_yd['delta_r'])))
        t = np.arange(0, len(input) * 0.1, 0.1)
        V0 = df_dutch_roll_yd['tas'][0]
        beta0 = 0
        phi0 = df_dutch_roll_yd['phi'][0]
        theta0 = df_dutch_roll_yd['theta'][0]
        p0 = df_dutch_roll_yd['p'][0]
        r0 = df_dutch_roll_yd['r'][0]
        m = (df_dutch_roll_yd['m'][0] + df_dutch_roll_yd['m'][-1]) / 2
        rho0 = df_dutch_roll_yd['rho'][0]
        CL = calc_CL(m * constants.g * np.cos(theta0), V0, rho0)
        state0 = np.transpose(beta0, phi0, p0, r0)

        model = AircraftModel(AerodynamicParameters)
        A, B, C, D = model.get_state_space_matrices_asymmetric(m, V0, rho0, theta0, CL)
        sys = ml.ss(A, B, C, D)
        yout, t, xout = ml.lsim(sys, input, t, state0)
        result = np.hstack((t, np.transpose(yout)))
        df_result = pd.DataFrame(result, columns=['t', 'beta', 'phi', 'p', 'r'])
        return df_result

    def simulate_spiral(self, df_spiral):
        input = np.transpose(np.hstack((df_spiral["delta_a"], df_spiral['delta_r'])))
        t = np.arange(0, len(input) * 0.1, 0.1)
        V0 = df_spiral['tas'][0]
        beta0 = 0
        phi0 = df_spiral['phi'][0]
        theta0 = df_spiral['theta'][0]
        p0 = df_spiral['p'][0]
        r0 = df_spiral['r'][0]
        m = (df_spiral['m'][0] + df_spiral['m'][-1]) / 2
        rho0 = df_spiral['rho'][0]
        CL = calc_CL(m * constants.g * np.cos(theta0), V0, rho0)
        state0 = np.transpose(beta0, phi0, p0, r0)

        model = AircraftModel(AerodynamicParameters)
        A, B, C, D = model.get_state_space_matrices_asymmetric(m, V0, rho0, theta0, CL)
        sys = ml.ss(A, B, C, D)
        yout, t, xout = ml.lsim(sys, input, t, state0)
        result = np.hstack((t, np.transpose(yout)))
        df_result = pd.DataFrame(result, columns=['t', 'beta', 'phi', 'p', 'r'])
        return df_result

    def simulate_aperiodic_roll(self, df_aperiodic_roll):
        input = np.transpose(np.hstack((df_aperiodic_roll["delta_a"], df_aperiodic_roll['delta_r'])))
        t = np.arange(0, len(input) * 0.1, 0.1)
        V0 = df_aperiodic_roll['tas'][0]
        beta0 = 0
        phi0 = df_aperiodic_roll['phi'][0]
        theta0 = df_aperiodic_roll['theta'][0]
        p0 = df_aperiodic_roll['p'][0]
        r0 = df_aperiodic_roll['r'][0]
        m = (df_aperiodic_roll['m'][0] + df_aperiodic_roll['m'][-1]) / 2
        rho0 = df_aperiodic_roll['rho'][0]
        CL = calc_CL(m * constants.g * np.cos(theta0), V0, rho0)
        state0 = np.transpose(beta0, phi0, p0, r0)

        model = AircraftModel(AerodynamicParameters)
        A, B, C, D = model.get_state_space_matrices_asymmetric(m, V0, rho0, theta0, CL)
        sys = ml.ss(A, B, C, D)
        yout, t, xout = ml.lsim(sys, input, t, state0)
        result = np.hstack((t, np.transpose(yout)))
        df_result = pd.DataFrame(result, columns=['t', 'beta', 'phi', 'p', 'r'])
        return df_result

    def simulate_phugoid(self, df_phugoid):
        input = np.transpose(np.hstack((df_spiral["delta_a"], df_spiral['delta_r'])))
        t = np.arange(0, len(input) * 0.1, 0.1)
        V0 = df_spiral['tas'][0]
        alpha0 = df_spiral["alpha"][0]
        theta0 = df_spiral['theta'][0]
        q0 = df_spiral['q'][0]
        m = (df_spiral['m'][0] + df_spiral['m'][-1]) / 2
        rho0 = df_spiral['rho'][0]
        CL = calc_CL(m * constants.g * np.cos(theta0), V0, rho0)
        state0 = np.transpose(0, alpha0, theta0, q0)




if __name__ == '__main__':
    sim = Simulation(AircraftModel(AerodynamicParameters()))
    df_dutch_roll = FlightTest("data/B24").df_dutch_roll
    df_out = sim.simulate_dutch_roll(df_dutch_roll)
    fig, (ax1, ax2) = plt.subplots(2, 1)
    # ax1.axvline(ft.data_sheet.timestamp_phugoid)
    ax1.plot(df_out.index, df_out["u"])
    ax2.plot(df_out.index, df_out["alpha"])
    format_plot()




