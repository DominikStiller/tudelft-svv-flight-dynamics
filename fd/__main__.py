import sys

from fd.analysis.flight_test import FlightTest
from fd.simulation.aircraft_model import AircraftModel
from fd.simulation.simulation import Simulation
from fd.validation.comparison import SimulatedMeasuredComparison
from fd.validation.comparison_eigenvalues import EigenvalueComparison

if __name__ == "__main__":
    flight_test = FlightTest(sys.argv[1])
    print(flight_test.aerodynamic_parameters)
    flight_test.make_aerodynamic_plots()

    aircraft_model = AircraftModel(flight_test.aerodynamic_parameters)
    simulation = Simulation(aircraft_model)

    comparison = SimulatedMeasuredComparison(flight_test, simulation)
    comparison.run_simulations()
    comparison.plot_responses()

    comparison_eigenvalues = EigenvalueComparison(flight_test, aircraft_model)
    comparison_eigenvalues.compare()
