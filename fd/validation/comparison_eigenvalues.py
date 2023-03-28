from fd.analysis.flight_test import FlightTest
from fd.simulation.aircraft_model import AircraftModel
from fd.validation.eigenmotion_characteristics import (
    time_constant_aperiodic_roll,
    characteristics_dutch_roll,
    time_constant_spiral,
    characteristics_phugoid,
    characteristics_short_period,
)


class EigenvalueComparison:
    def __init__(self, flight_test: FlightTest, model: AircraftModel):
        self.flight_test = flight_test
        self.model = model

    def compare(self):
        self._compare_phugoid()
        self._compare_short_period()
        self._compare_dutch_roll()
        self._compare_aperiodic_roll()
        self._compare_spiral()

    def _compare_phugoid(self):
        A, _, _, _ = self.model.get_state_space_matrices_symmetric_from_df(
            self.flight_test.df_phugoid
        )
        eigs, _ = self.model.get_eigenvalues_and_eigenvectors(A)
        eig_phugoid = eigs[2]
        V0 = self.flight_test.df_phugoid["tas"].iloc[0]
        P, T_half = characteristics_phugoid(eig_phugoid, V0)
        print(f"Phugoid: P = {P:.2} s, T_half = {T_half:.2} s")

    def _compare_short_period(self):
        A, _, _, _ = self.model.get_state_space_matrices_symmetric_from_df(
            self.flight_test.df_short_period
        )
        eigs, _ = self.model.get_eigenvalues_and_eigenvectors(A)
        eig_short_period = eigs[0]
        V0 = self.flight_test.df_short_period["tas"].iloc[0]
        P, T_half = characteristics_short_period(eig_short_period, V0)
        print(f"Short period: P = {P:.2} s, T_half = {T_half:.2} s")

    def _compare_dutch_roll(self):
        A, _, _, _ = self.model.get_state_space_matrices_asymmetric_from_df(
            self.flight_test.df_dutch_roll
        )
        eigs, _ = self.model.get_eigenvalues_and_eigenvectors(A)
        eig_dutch_roll = eigs[1]
        V0 = self.flight_test.df_dutch_roll["tas"].iloc[0]
        P, T_half = characteristics_dutch_roll(eig_dutch_roll, V0)
        print(f"Dutch roll: P = {P:.2} s, T_half = {T_half:.2} s")

    def _compare_aperiodic_roll(self):
        A, _, _, _ = self.model.get_state_space_matrices_asymmetric_from_df(
            self.flight_test.df_aperiodic_roll
        )
        eigs, _ = self.model.get_eigenvalues_and_eigenvectors(A)
        eig_aperiodic_roll = eigs[0]
        V0 = self.flight_test.df_aperiodic_roll["tas"].iloc[0]
        tau = time_constant_aperiodic_roll(eig_aperiodic_roll, V0)
        print(f"Aperiodic roll: tau = {tau:.2}")

    def _compare_spiral(self):
        A, _, _, _ = self.model.get_state_space_matrices_asymmetric_from_df(
            self.flight_test.df_spiral
        )
        eigs, _ = self.model.get_eigenvalues_and_eigenvectors(A)
        V0 = self.flight_test.df_spiral["tas"].iloc[0]
        eig_spiral = eigs[3]
        tau = time_constant_spiral(eig_spiral, V0)
        print(f"Spiral: tau = {tau:.2}")
