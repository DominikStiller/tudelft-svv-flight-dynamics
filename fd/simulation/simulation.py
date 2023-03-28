import control.matlab as ml
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas import DataFrame

from fd.analysis.aerodynamics import calc_CL
from fd.analysis.flight_test import FlightTest
from fd.plotting import format_plot
from fd.simulation import constants
from fd.simulation.aircraft_model import AircraftModel
from fd.structs import AerodynamicParameters


class Simulation:
    def __init__(self, model: AircraftModel):
        self.model = model

    def simulate_dutch_roll(self, df_dutch_roll) -> DataFrame:
        data = df_dutch_roll
        delta_a = -(data["delta_a"] - data["delta_a"].iloc[0])
        delta_r = -(data["delta_r"] - data["delta_a"].iloc[0])
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
        state0 = np.array([0, 0, 0, 0])

        A, B, C, D = self.model.get_state_space_matrices_asymmetric(m, V0, rho0, theta0, CL)
        eigen = self.model.get_eigenvalues_and_eigenvectors(A)[0]
        sys = ml.ss(A, B, C, D)
        yout, t, xout = ml.lsim(sys, input, t, state0)
        yout += np.array([0, phi0, p0, r0])
        result = np.hstack((np.transpose(t).reshape((len(t), 1)), yout))
        df_result = pd.DataFrame(result, columns=["t", "beta", "phi", "p", "r"])
        df_result = df_result.set_index("t", drop=True)

        return df_result

    def simulate_dutch_roll_yd(self, df_dutch_roll_yd):
        data = df_dutch_roll_yd
        delta_a = -(data["delta_a"] - data["delta_a"].iloc[0])
        delta_r = -(data["delta_r"] - data["delta_a"].iloc[0])
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
        state0 = np.array([0, 0, 0, 0])

        A, B, C, D = self.model.get_state_space_matrices_asymmetric(m, V0, rho0, theta0, CL)
        eigen = self.model.get_eigenvalues_and_eigenvectors(A)[0]
        print(eigen)
        sys = ml.ss(A, B, C, D)
        yout, t, xout = ml.lsim(sys, input, t, state0)
        yout += np.array([0, phi0, p0, r0])
        result = np.hstack((np.transpose(t).reshape((len(t), 1)), yout))
        df_result = pd.DataFrame(result, columns=["t", "beta", "phi", "p", "r"])
        df_result = df_result.set_index("t", drop=True)
        return df_result

    def simulate_spiral(self, df_spiral):
        data = df_spiral
        delta_a = data["delta_a"] - data["delta_a"].iloc[0]
        delta_r = data["delta_r"] - data["delta_a"].iloc[0]
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
        state0 = np.array([0, 0, 0, 0])

        model = AircraftModel(AerodynamicParameters)
        A, B, C, D = model.get_state_space_matrices_asymmetric(m, V0, rho0, theta0, CL)
        print(np.linalg.eig(A)[0])
        sys = ml.ss(A, B, C, D)
        yout, t, xout = ml.lsim(sys, input, t, state0)
        yout += np.array([0, phi0, p0, r0])
        result = np.hstack((np.transpose(t).reshape((len(t), 1)), yout))
        df_result = pd.DataFrame(result, columns=["t", "beta", "phi", "p", "r"])
        df_result = df_result.set_index("t", drop=True)

        return df_result

    def simulate_aperiodic_roll(self, df_aperiodic_roll):
        data = df_aperiodic_roll
        delta_a = data["delta_a"] - data["delta_a"].iloc[0]
        delta_r = data["delta_r"] - data["delta_a"].iloc[0]
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
        state0 = np.array([0, 0, 0, 0])

        model = AircraftModel(AerodynamicParameters)
        A, B, C, D = model.get_state_space_matrices_asymmetric(m, V0, rho0, theta0, CL)
        sys = ml.ss(A, B, C, D)
        yout, t, xout = ml.lsim(sys, input, t, state0)
        yout += np.array([0, phi0, p0, r0])
        result = np.hstack((np.transpose(t).reshape((len(t), 1)), yout))
        df_result = pd.DataFrame(result, columns=["t", "beta", "phi", "p", "r"])
        df_result = df_result.set_index("t", drop=True)
        return df_result

    def simulate_phugoid(self, df_phugoid):
        data = df_phugoid
        delta_e = data["delta_e"] - data["delta_e"].iloc[0]  # - data["delta_e"].iloc[0]
        t = data.index
        input = delta_e  # .reshape((len(t), 1))
        V0 = data["tas"].iloc[0]
        u_hat0 = 0
        alpha0 = data["alpha"].iloc[0]
        theta0 = data["theta"].iloc[0]
        q0 = data["q"].iloc[0]
        m = (data["m"].iloc[0] + data["m"].iloc[-1]) / 2
        rho0 = data["rho"].iloc[0]
        CL = calc_CL(m * constants.g * np.cos(theta0), V0, rho0)
        state0 = np.array([0, 0, 0, 0])  # np.array([u_hat0, alpha0, theta0, q0])

        A, B, C, D = self.model.get_state_space_matrices_symmetric(m, V0, rho0, theta0)
        print(np.linalg.eig(A)[0])
        print(np.linalg.eig(A)[1])
        sys = ml.ss(A, B, C, D)
        print(self.model.get_eigenvalues_and_eigenvectors(A)[0])
        yout, t, xout = ml.lsim(sys, input, t, state0)
        yout += np.array([u_hat0, alpha0, theta0, q0])
        result = np.hstack((np.transpose(t).reshape((len(t), 1)), yout))
        df_result = pd.DataFrame(result, columns=["t", "u_hat", "alpha", "theta", "q"])
        df_result = df_result.set_index("t", drop=True)
        return df_result

    def simulate_short_period(self, df_short_period):
        data = df_short_period
        delta_e = data["delta_e"] - data["delta_e"].iloc[0]  # - data["delta_e"].iloc[0]
        t = data.index
        input = delta_e  # .reshape((len(t), 1))
        V0 = data["tas"].iloc[0]
        u_hat0 = 0
        alpha0 = data["alpha"].iloc[0]
        theta0 = data["theta"].iloc[0]
        q0 = data["q"].iloc[0]
        m = (data["m"].iloc[0] + data["m"].iloc[-1]) / 2
        rho0 = data["rho"].iloc[0]
        CL = calc_CL(m * constants.g * np.cos(theta0), V0, rho0)
        state0 = np.array([0, 0, 0, 0])  # np.array([u_hat0, alpha0, theta0, q0])

        A, B, C, D = self.model.get_state_space_matrices_symmetric(m, V0, rho0, theta0)
        print(np.linalg.eig(A)[0])
        print(np.linalg.eig(A)[1])
        sys = ml.ss(A, B, C, D)
        print(self.model.get_eigenvalues_and_eigenvectors(A)[0])
        yout, t, xout = ml.lsim(sys, input, t, state0)
        yout += np.array([u_hat0, alpha0, theta0, q0])
        result = np.hstack((np.transpose(t).reshape((len(t), 1)), yout))
        df_result = pd.DataFrame(result, columns=["t", "u_hat", "alpha", "theta", "q"])
        df_result = df_result.set_index("t", drop=True)
        return df_result


if __name__ == "__main__":
    sim = Simulation(
        AircraftModel(
            AerodynamicParameters(
                C_L_alpha=4.758556374647304,
                alpha_0=-0.023124783070063493,
                C_D_0=0.023439123324849084,
                # C_m_alpha=-0.5554065208385275,
                C_m_alpha=-0.5,
                C_m_delta=-1.3380975545274032,
                e=1.0713238368125688,
            )
        )
    )
    df = FlightTest("data/B24").df_phugoid
    df_out = sim.simulate_phugoid(df)
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1)

    y1 = "tas"
    y2 = "alpha"
    y3 = "theta"
    y4 = "q"
    """
    y1 = "beta"
    y2 = "phi"
    y3 = "p"
    y4 = "r"
    """

    ax1.plot(df_out.index, df_out["u_hat"] * df["tas"].iloc[0] + df["tas"].iloc[0])
    # ax1.plot(df_out.index, df_out[y1])
    ax1.plot(df_out.index, df[y1], color="black")
    ax1.set_ylabel(y1)
    ax2.plot(df_out.index, df_out[y2])
    ax2.plot(df_out.index, df[y2], color="black")
    # ax2.set_ylim(-0.2, 2.7)
    ax2.set_ylabel(y2)
    ax3.plot(df_out.index, df_out[y3])
    ax3.plot(df_out.index, df[y3], color="black")
    # ax3.set_ylim(-0.3, 0.25)
    ax3.set_ylabel(y3)
    ax4.plot(df_out.index, df_out[y4])
    ax4.plot(df_out.index, df[y4], color="black")
    # ax4.set_ylim(-0.25, 0.3)
    ax4.set_ylabel(y4)
    # ax5.plot(df_out.index, df['delta_'])
    ax4.set_xlabel("t")

    format_plot()
    plt.show()
