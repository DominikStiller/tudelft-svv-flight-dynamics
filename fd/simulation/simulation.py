import control.matlab as ml
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from fd.analysis.aerodynamics import calc_CL
from fd.analysis.flight_test import FlightTest
from fd.plotting import format_plot
from fd.simulation import constants
from fd.simulation.aircraft_model import AircraftModel
from fd.structs import AerodynamicParameters


class Simulation:
    def __init__(self, model: AircraftModel):
        self.model = model

    def simulate_dutch_roll(self, df_dutch_roll) -> pd.DataFrame:
        delta_a = df_dutch_roll["delta_a"]
        delta_r = df_dutch_roll["delta_r"]
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
        df_result = df_result.set_index("t", drop=True)
        return df_result


if __name__ == "__main__":
    sim = Simulation(
        AircraftModel(
            AerodynamicParameters(
                C_L_alpha=0.08305247637936027,
                alpha_0=-1.3249524720702168,
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
