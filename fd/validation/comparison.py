from fd.analysis.flight_test import FlightTest
from fd.simulation.simulation import Simulation


class SimulatedMeasuredComparison:
    def __init__(self, flight_test: FlightTest, simulation: Simulation):
        self.flight_test = flight_test
        self.simulation = simulation

    def run_simulations(self):
        self.simulation_output_dutch_roll = self.simulation.simulate_dutch_roll()

    def plot_responses(self):
        self._plot_dutch_roll()

    def _plot_dutch_roll(self):
        pass
