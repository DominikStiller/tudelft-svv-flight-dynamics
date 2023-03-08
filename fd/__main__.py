from fd.analysis.flight_test import FlightTest
from fd.simulation.aircraft_model import AircraftModel
from fd.simulation.simulation import Simulation
from fd.validation.comparison import SimulatedMeasuredComparison

if __name__ == "__main__":
    flight_test = FlightTest("data/dataset")
    aero_params = flight_test.analyze()

    aircraft_model = AircraftModel(aero_params)
    simulation = Simulation(aircraft_model)

    comparison = SimulatedMeasuredComparison(flight_test, simulation)
    comparison.run_simulations()
    comparison.plot_responses()
