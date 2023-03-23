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
        data = df_dutch_roll
        delta_a = data['delta_a']
        delta_r = data['delta_r']
        input = np.column_stack((delta_a, delta_r))
        t = data.index
        V0 = data["tas"].iloc[0]
        beta0 = 0
        phi0 = data["phi"].iloc[0]
        theta0 = data["theta"].iloc[0]
        p0 = data["p"].iloc[0]
        r0 = data["r"].iloc[0]
        m = (data["m"].iloc[0] + data["m"].iloc[-1]) / 2
        rho0 = data["rho"].iloc[0]
        CL = calc_CL(m * constants.g * np.cos(theta0), V0, rho0)
        state0 = np.array([beta0, phi0, p0, r0])

        A, B, C, D = self.model.get_state_space_matrices_asymmetric(m, V0, rho0, theta0, CL)
        eigen = self.model.get_eigenvalues_and_eigenvectors(A)[0]
        sys = ml.ss(A, B, C, D)
        yout, t, xout = ml.lsim(sys, input, t, state0)
        result = np.hstack((np.transpose(t).reshape((len(t), 1)), yout))
        df_result = pd.DataFrame(result, columns=["t", "beta", "phi", "p", "r"])
        return df_result

    def simulate_dutch_roll_yd(self, df_dutch_roll_yd):
        data = df_dutch_roll_yd
        delta_a = data['delta_a']
        delta_r = data['delta_r']
        input = np.column_stack((delta_a, delta_r))
        t = data.index
        V0 = data["tas"].iloc[0]
        beta0 = 0
        phi0 = data["phi"].iloc[0]
        theta0 = data["theta"].iloc[0]
        p0 = data["p"].iloc[0]
        r0 = data["r"].iloc[0]
        m = (data["m"].iloc[0] + data["m"].iloc[-1]) / 2
        rho0 = data["rho"].iloc[0]
        CL = calc_CL(m * constants.g * np.cos(theta0), V0, rho0)
        state0 = np.array([beta0, phi0, p0, r0])

        A, B, C, D = self.model.get_state_space_matrices_asymmetric(m, V0, rho0, theta0, CL)
        eigen = self.model.get_eigenvalues_and_eigenvectors(A)[0]
        print(eigen)
        sys = ml.ss(A, B, C, D)
        yout, t, xout = ml.lsim(sys, input, t, state0)
        result = np.hstack((np.transpose(t).reshape((len(t), 1)), yout))
        df_result = pd.DataFrame(result, columns=["t", "beta", "phi", "p", "r"])
        return df_result

    def simulate_spiral(self, df_spiral):
        data = df_spiral
        delta_a = data['delta_a']
        delta_r = data['delta_r']
        input = np.column_stack((delta_a, delta_r))
        t = data.index
        V0 = data["tas"].iloc[0]
        beta0 = 0
        phi0 = data["phi"].iloc[0]
        theta0 = data["theta"].iloc[0]
        p0 = data["p"].iloc[0]
        r0 = data["r"].iloc[0]
        m = (data["m"].iloc[0] + data["m"].iloc[-1]) / 2
        rho0 = data["rho"].iloc[0]
        CL = calc_CL(m * constants.g * np.cos(theta0), V0, rho0)
        state0 = np.array([beta0, phi0, p0, r0])

        model = AircraftModel(AerodynamicParameters)
        A, B, C, D = model.get_state_space_matrices_asymmetric(m, V0, rho0, theta0, CL)
        sys = ml.ss(A, B, C, D)
        yout, t, xout = ml.lsim(sys, input, t, state0)
        result = np.hstack((np.transpose(t).reshape((len(t), 1)), yout))
        df_result = pd.DataFrame(result, columns=["t", "beta", "phi", "p", "r"])
        return df_result

    def simulate_aperiodic_roll(self, df_aperiodic_roll):
        data = df_aperiodic_roll
        delta_a = data['delta_a']
        delta_r = data['delta_r']
        input = np.column_stack((delta_a, delta_r))
        t = data.index
        V0 = data["tas"].iloc[0]
        beta0 = 0
        phi0 = data["phi"].iloc[0]
        theta0 = data["theta"].iloc[0]
        p0 = data["p"].iloc[0]
        r0 = data["r"].iloc[0]
        m = (data["m"].iloc[0] + data["m"].iloc[-1]) / 2
        rho0 = data["rho"].iloc[0]
        CL = calc_CL(m * constants.g * np.cos(theta0), V0, rho0)
        state0 = np.array([beta0, phi0, p0, r0])

        model = AircraftModel(AerodynamicParameters)
        A, B, C, D = model.get_state_space_matrices_asymmetric(m, V0, rho0, theta0, CL)
        sys = ml.ss(A, B, C, D)
        yout, t, xout = ml.lsim(sys, input, t, state0)
        result = np.hstack((np.transpose(t).reshape((len(t), 1)), yout))
        df_result = pd.DataFrame(result, columns=["t", "beta", "phi", "p", "r"])
        return df_result

    def simulate_phugoid(self, df_phugoid):
        data = df_phugoid
        delta_e = data['delta_e']
        input = delta_e
        t = data.index
        V0 = data["tas"].iloc[0]
        beta0 = 0
        u_hat = data["phi"].iloc[0]
        theta0 = data["theta"].iloc[0]
        p0 = data["p"].iloc[0]
        r0 = data["r"].iloc[0]
        m = (data["m"].iloc[0] + data["m"].iloc[-1]) / 2
        rho0 = data["rho"].iloc[0]
        CL = calc_CL(m * constants.g * np.cos(theta0), V0, rho0)
        state0 = np.array([beta0, phi0, p0, r0])

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
    df_spiral = FlightTest("data/B24").df_spiral
    df_out = sim.simulate_spiral(df_spiral)
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1)
    # ax1.axvline(ft.data_sheet.timestamp_phugoid)
    ax1.plot(df_out.index, df_out["beta"])
    ax1.set_ylabel("beta")
    ax2.plot(df_out.index, df_out["phi"])
    ax2.plot(df_out.index, df_spiral["phi"])
    ax2.set_ylabel("phi")
    ax3.plot(df_out.index, df_out["p"])
    ax3.plot(df_out.index, df_spiral["phi"])
    ax3.set_ylabel("p")
    ax4.plot(df_out.index, df_out["r"])
    ax4.plot(df_out.index, df_spiral["phi"])
    ax4.set_ylabel("r")
    ax4.set_xlabel("t")



    format_plot()
    plt.show()
