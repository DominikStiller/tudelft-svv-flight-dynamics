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
        delta_a = df_dutch_roll['delta_a']
        delta_r = df_dutch_roll['delta_r']
        input = np.column_stack((delta_a, delta_r))
        t = np.arange(0, len(input) * 0.1, 0.1)
        V0 = df_dutch_roll["tas"].iloc[0]
        beta0 = 0
        phi0 = df_dutch_roll["phi"].iloc[0]
        theta0 = df_dutch_roll["theta"].iloc[0]
        p0 = df_dutch_roll["p"].iloc[0]
        r0 = df_dutch_roll["r"].iloc[0]
        m = (df_dutch_roll["m"].iloc[0] + df_dutch_roll["m"].iloc[-1]) / 2
        rho0 = df_dutch_roll["rho"].iloc[0]
        CL = calc_CL(m * constants.g * np.cos(theta0), V0, rho0)
        state0 = np.array([beta0, phi0, p0, r0])

        A, B, C, D = self.model.get_state_space_matrices_asymmetric(m, V0, rho0, theta0, CL)
        sys = ml.ss(A, B, C, D)
        yout, t, xout = ml.lsim(sys, input, t, state0)
        result = np.hstack((np.transpose(t).reshape((len(t), 1)), yout))
        df_result = pd.DataFrame(result, columns=["t", "beta", "phi", "p", "r"])
        return df_result

    def simulate_dutch_roll_yd(self, df_dutch_roll_yd):
        input = np.transpose(np.hstack((df_dutch_roll_yd["delta_a"], df_dutch_roll_yd["delta_r"])))
        t = np.arange(0, len(input) * 0.1, 0.1)
        V0 = df_dutch_roll_yd["tas"][0]
        beta0 = 0
        phi0 = df_dutch_roll_yd["phi"][0]
        theta0 = df_dutch_roll_yd["theta"][0]
        p0 = df_dutch_roll_yd["p"][0]
        r0 = df_dutch_roll_yd["r"][0]
        m = (df_dutch_roll_yd["m"][0] + df_dutch_roll_yd["m"][-1]) / 2
        rho0 = df_dutch_roll_yd["rho"][0]
        CL = calc_CL(m * constants.g * np.cos(theta0), V0, rho0)
        state0 = np.transpose(beta0, phi0, p0, r0)

        model = AircraftModel(AerodynamicParameters)
        A, B, C, D = model.get_state_space_matrices_asymmetric(m, V0, rho0, theta0, CL)
        sys = ml.ss(A, B, C, D)
        yout, t, xout = ml.lsim(sys, input, t, state0)
        result = np.hstack((t, np.transpose(yout)))
        df_result = pd.DataFrame(result, columns=["t", "beta", "phi", "p", "r"])
        return df_result

    def simulate_spiral(self, df_spiral):
        input = np.transpose(np.hstack((df_spiral["delta_a"], df_spiral["delta_r"])))
        t = np.arange(0, len(input) * 0.1, 0.1)
        V0 = df_spiral["tas"][0]
        beta0 = 0
        phi0 = df_spiral["phi"][0]
        theta0 = df_spiral["theta"][0]
        p0 = df_spiral["p"][0]
        r0 = df_spiral["r"][0]
        m = (df_spiral["m"][0] + df_spiral["m"][-1]) / 2
        rho0 = df_spiral["rho"][0]
        CL = calc_CL(m * constants.g * np.cos(theta0), V0, rho0)
        state0 = np.transpose(beta0, phi0, p0, r0)

        model = AircraftModel(AerodynamicParameters)
        A, B, C, D = model.get_state_space_matrices_asymmetric(m, V0, rho0, theta0, CL)
        sys = ml.ss(A, B, C, D)
        yout, t, xout = ml.lsim(sys, input, t, state0)
        result = np.hstack((t, np.transpose(yout)))
        df_result = pd.DataFrame(result, columns=["t", "beta", "phi", "p", "r"])
        return df_result

    def simulate_aperiodic_roll(self, df_aperiodic_roll):
        input = np.transpose(
            np.hstack((df_aperiodic_roll["delta_a"], df_aperiodic_roll["delta_r"]))
        )
        t = np.arange(0, len(input) * 0.1, 0.1)
        V0 = df_aperiodic_roll["tas"][0]
        beta0 = 0
        phi0 = df_aperiodic_roll["phi"][0]
        theta0 = df_aperiodic_roll["theta"][0]
        p0 = df_aperiodic_roll["p"][0]
        r0 = df_aperiodic_roll["r"][0]
        m = (df_aperiodic_roll["m"][0] + df_aperiodic_roll["m"][-1]) / 2
        rho0 = df_aperiodic_roll["rho"][0]
        CL = calc_CL(m * constants.g * np.cos(theta0), V0, rho0)
        state0 = np.transpose(beta0, phi0, p0, r0)

        model = AircraftModel(AerodynamicParameters)
        A, B, C, D = model.get_state_space_matrices_asymmetric(m, V0, rho0, theta0, CL)
        sys = ml.ss(A, B, C, D)
        yout, t, xout = ml.lsim(sys, input, t, state0)
        result = np.hstack((t, np.transpose(yout)))
        df_result = pd.DataFrame(result, columns=["t", "beta", "phi", "p", "r"])
        return df_result

    def simulate_phugoid(self, df_phugoid):
        input = np.transpose(np.hstack((df_phugoid["delta_a"], df_phugoid["delta_r"])))
        t = np.arange(0, len(input) * 0.1, 0.1)
        V0 = df_phugoid["tas"][0]
        alpha0 = df_phugoid["alpha"][0]
        theta0 = df_phugoid["theta"][0]
        q0 = df_phugoid["q"][0]
        m = (df_phugoid["m"][0] + df_phugoid["m"][-1]) / 2
        rho0 = df_phugoid["rho"][0]
        CL = calc_CL(m * constants.g * np.cos(theta0), V0, rho0)
        state0 = np.transpose(0, alpha0, theta0, q0)

        model = AircraftModel(AerodynamicParameters)
        A, B, C, D = model.get_state_space_matrices_asymmetric(m, V0, rho0, theta0, CL)
        sys = ml.ss(A, B, C, D)
        yout, t, xout = ml.lsim(sys, input, t, state0)
        result = np.hstack((t, np.transpose(yout)))
        df_result = pd.DataFrame(result, columns=["t", "beta", "phi", "p", "r"])
        return df_result


if __name__ == "__main__":
    sim = Simulation(
        AircraftModel(
            AerodynamicParameters(
                C_L_alpha=0.08305247637936027,
                C_D_0=0.023439123324849084,
                C_m_alpha=-0.009693672475678806,
                C_m_delta=-0.023354208039387547,
                e=1.0713238368125688,
            )
        )
    )
    df_dutch_roll = FlightTest("data/B24").df_dutch_roll
    df_out = sim.simulate_dutch_roll(df_dutch_roll)
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1)
    # ax1.axvline(ft.data_sheet.timestamp_phugoid)
    ax1.plot(df_out["t"], df_out["beta"])
    ax1.set_ylabel("beta")
    ax2.plot(df_out["t"], df_out["phi"])
    ax2.set_ylabel("phi")
    ax3.plot(df_out["t"], df_out["p"])
    ax3.set_ylabel("p")
    ax4.plot(df_out["t"], df_out["r"])
    ax4.set_ylabel("r")
    ax4.set_xlabel("t")


    format_plot()
    plt.show()
