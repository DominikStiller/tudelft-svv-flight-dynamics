from fd.simulation.aircraft_model import AircraftModel
from fd.structs import SimulationOutput


class Simulation:
    def __init__(self, model: AircraftModel):
        self.model = model

    def simulate_dutch_roll(self) -> SimulationOutput:
        pass
