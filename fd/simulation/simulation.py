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

    def simulate_asymmetric(self, data, flip_input=False) -> DataFrame:
        t = data.index

        delta_a = data["delta_a"] - data["delta_a"].iloc[0]
        delta_r = data["delta_r"] - data["delta_a"].iloc[0]
        input = np.column_stack((delta_a, delta_r))
        if flip_input:
            input *= -1

        phi0 = data["phi"].iloc[0]
        p0 = data["p"].iloc[0]
        r0 = data["r"].iloc[0]
        state0_absolute = np.array([0, phi0, p0, r0])

        ABCD = self.model.get_state_space_matrices_asymmetric_from_df(data)

        sys = ml.ss(*ABCD)
        yout, t, xout = ml.lsim(sys, input, t)
        yout += state0_absolute
        result = np.hstack((np.transpose(t).reshape((len(t), 1)), yout))
        df_result = pd.DataFrame(result, columns=["t", "beta", "phi", "p", "r"])
        df_result = df_result.set_index("t", drop=True)

        return df_result

    def simulate_symmetric(self, data):
        t = data.index

        delta_e = data["delta_e"] - data["delta_e"].iloc[0]
        input = delta_e

        theta0 = data["theta"].iloc[0]
        u_hat0 = 0
        alpha0 = data["alpha"].iloc[0]
        q0 = data["q"].iloc[0]
        state0_absolute = np.array([u_hat0, alpha0, theta0, q0])

        ABCD = self.model.get_state_space_matrices_symmetric_from_df(data)

        sys = ml.ss(*ABCD)
        yout, t, xout = ml.lsim(sys, input, t)
        yout += state0_absolute
        result = np.hstack((np.transpose(t).reshape((len(t), 1)), yout))
        df_result = pd.DataFrame(result, columns=["t", "u_hat", "alpha", "theta", "q"])
        df_result = df_result.set_index("t", drop=True)

        return df_result


if __name__ == "__main__":
    # sim = Simulation(
    #     AircraftModel(
    #         AerodynamicParameters(
    #             C_L_alpha=4.758556374647304,
    #             alpha_0=-0.023124783070063493,
    #             C_D_0=0.023439123324849084,
    #             # C_m_alpha=-0.5554065208385275,
    #             C_m_alpha=-0.5,
    #             C_m_delta=-1.3380975545274032,
    #             e=1.0713238368125688,
    #         )
    #     )
    # )
    ft = FlightTest("data/B24")
    df = ft.df_spiral
    aircraft_model = AircraftModel(ft.aerodynamic_parameters)
    sim = Simulation(aircraft_model)
    df_out = sim.simulate_asymmetric(df, flip_input=False)
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1)
    """
    y1 = "tas"
    y2 = "alpha"
    y3 = "theta"
    y4 = "q"
    """
    y1 = "beta"
    y2 = "phi"
    y3 = "p"
    y4 = "r"

    # ax1.plot(df_out.index, df_out["u_hat"] * df["tas"].iloc[0] + df["tas"].iloc[0])
    ax1.plot(df_out.index, df_out[y1])
    # ax1.plot(df_out.index, df[y1], color="black")
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
    ax4.set_ylim(-0.25, 0.3)
    ax4.set_ylabel(y4)
    # ax5.plot(df_out.index, df['delta_'])
    ax4.set_xlabel("t")

    format_plot()
    plt.show()
