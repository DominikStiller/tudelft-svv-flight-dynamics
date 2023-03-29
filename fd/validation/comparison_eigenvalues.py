import numpy as np

from fd.analysis.aerodynamics import calc_CL
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
        data = self.flight_test.df_phugoid
        A, _, _, _ = self.model.get_state_space_matrices_symmetric_from_df(data)
        eigs, _ = self.model.get_eigenvalues_and_eigenvectors(A)
        eig_phugoid = eigs[2]
        P, T_half = characteristics_phugoid(eig_phugoid)
        print(f"Phugoid (simulated): P = {P:.3f} s, T_half = {T_half:.3f} s, {eig_phugoid}")

        m = (data["m"].iloc[0] + data["m"].iloc[-1]) / 2
        V0 = data["tas"].iloc[0]
        rho0 = data["rho"].iloc[0]
        theta0 = data["theta"].iloc[0]
        eig_phugoid = self.model.get_idealized_phugoid_eigenvalues(m, rho0, theta0, V0)[0]
        P, T_half = characteristics_phugoid(eig_phugoid)
        print(f"Phugoid (idealized): P = {P:.3f} s, T_half = {T_half:.3f} s")

    def _compare_short_period(self):
        data = self.flight_test.df_short_period
        A, _, _, _ = self.model.get_state_space_matrices_symmetric_from_df(data)
        eigs, _ = self.model.get_eigenvalues_and_eigenvectors(A)
        eig_short_period = eigs[0]
        P, T_half = characteristics_short_period(eig_short_period)
        print(
            f"Short period (simulated): P = {P:.3f} s, T_half = {T_half:.3f} s, {eig_short_period}"
        )

        m = (data["m"].iloc[0] + data["m"].iloc[-1]) / 2
        V0 = data["tas"].iloc[0]
        rho0 = data["rho"].iloc[0]
        eig_short_period = self.model.get_idealized_shortperiod_eigenvalues(m, rho0, V0)[0]
        P, T_half = characteristics_short_period(eig_short_period)
        print(f"Short period (idealized): P = {P:.3f} s, T_half = {T_half:.3f} s")

    def _compare_dutch_roll(self):
        data = self.flight_test.df_dutch_roll
        A, _, _, _ = self.model.get_state_space_matrices_asymmetric_from_df(data)
        eigs, _ = self.model.get_eigenvalues_and_eigenvectors(A)
        eig_dutch_roll = eigs[1]
        P, T_half = characteristics_dutch_roll(eig_dutch_roll)
        print(f"Dutch roll (simulated): P = {P:.3f} s, T_half = {T_half:.3f} s, {eig_dutch_roll}")

        m = (data["m"].iloc[0] + data["m"].iloc[-1]) / 2
        V0 = data["tas"].iloc[0]
        rho0 = data["rho"].iloc[0]
        eig_dutch_roll = self.model.get_idealized_dutchroll_eigenvalues(m, rho0, V0)[0]
        P, T_half = characteristics_dutch_roll(eig_dutch_roll)
        print(f"Dutch roll (idealized): P = {P:.3f} s, T_half = {T_half:.3f} s")

    def _compare_aperiodic_roll(self):
        data = self.flight_test.df_aperiodic_roll
        A, _, _, _ = self.model.get_state_space_matrices_asymmetric_from_df(data)
        eigs, _ = self.model.get_eigenvalues_and_eigenvectors(A)
        eig_aperiodic_roll = eigs[0]
        V0 = data["tas"].iloc[0]
        tau = time_constant_aperiodic_roll(eig_aperiodic_roll, V0)
        print(f"Aperiodic roll (simulated): tau = {tau:.3f}, {eig_aperiodic_roll}")

        m = (data["m"].iloc[0] + data["m"].iloc[-1]) / 2
        rho0 = data["rho"].iloc[0]
        eig_aperiodic_roll = self.model.get_idealized_aperiodicroll_eigenvalues(m, rho0, V0)
        tau = time_constant_aperiodic_roll(eig_aperiodic_roll, V0)
        print(f"Aperiodic roll (idealized): tau = {tau:.3f}")

    def _compare_spiral(self):
        data = self.flight_test.df_spiral
        A, _, _, _ = self.model.get_state_space_matrices_asymmetric_from_df(data)
        eigs, _ = self.model.get_eigenvalues_and_eigenvectors(A)
        eig_spiral = eigs[3]
        tau = time_constant_spiral(eig_spiral)
        print(f"Spiral (simulated): tau = {tau:.3f}, {eig_spiral}")

        m = (data["m"].iloc[0] + data["m"].iloc[-1]) / 2
        V0 = data["tas"].iloc[0]
        rho0 = data["rho"].iloc[0]
        theta0 = data["theta"].iloc[0]
        CL = calc_CL(data["W"].iloc[0] * np.cos(theta0), V0, rho0)
        eig_aperiodic_roll = self.model.get_idealized_spiral_eigenvalues(m, rho0, V0, CL)
        tau = time_constant_spiral(eig_aperiodic_roll)
        print(f"Spiral (idealized): tau = {tau:.3f}")
