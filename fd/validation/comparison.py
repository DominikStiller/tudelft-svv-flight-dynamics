import pandas as pd
from matplotlib import pyplot as plt

from fd.analysis.flight_test import FlightTest
from fd.plotting import format_plot
from fd.simulation.simulation import Simulation


class SimulatedMeasuredComparison:
    simulated_dutch_roll: pd.DataFrame

    def __init__(self, flight_test: FlightTest, simulation: Simulation):
        self.flight_test = flight_test
        self.simulation = simulation

    def run_simulations(self):
        self.simulated_dutch_roll = self.simulation.simulate_dutch_roll(
            self.flight_test.df_dutch_roll
        )

    def plot_responses(self):
        self._plot_dutch_roll()

    def _plot_dutch_roll(self):
        fig, (ax_p, ax_r) = plt.subplots(2, 1, figsize=(12, 6))

        ax_p.plot(
            self.simulated_dutch_roll.index,
            self.simulated_dutch_roll["p"],
            label="Simulated",
        )
        ax_p.plot(
            self.flight_test.df_dutch_roll.index,
            self.flight_test.df_dutch_roll["r"],
            label="Measured",
        )
        ax_p.set_ylabel("Roll rate $p$ [°/s]")
        ax_p.legend()

        ax_r.plot(
            self.simulated_dutch_roll.index,
            self.simulated_dutch_roll["r"],
            label="Simulated",
        )
        ax_r.plot(
            self.flight_test.df_dutch_roll.index,
            self.flight_test.df_dutch_roll["p"],
            label="Measured",
        )
        ax_r.set_xlabel("Time [s]")
        ax_r.set_ylabel("Yaw rate $r$ [°/s]")
        ax_r.legend()

        format_plot()
        plt.show()
