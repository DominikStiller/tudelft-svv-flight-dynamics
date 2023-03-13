import unittest
import numpy.testing as test
from B24.fd.simulation.constants import *
from B24.fd.structs import AerodynamicParameters
from B24.fd.simulation.aircraft_model import AircraftModel

class MyTestCase(unittest.TestCase):

    def __init__(self, aero_params: AerodynamicParameters):
        self.aero_params = aero_params
    def test_symmetric_matrices(self):
        model = AircraftModel
        A, B, C, D = model.get_state_space_matrices_symmetric()
        self.assertEqual(True, False)  # add assertion here


if __name__ == '__main__':
    unittest.main()
