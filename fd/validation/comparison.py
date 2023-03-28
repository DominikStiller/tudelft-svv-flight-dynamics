import pandas as pd
from matplotlib import pyplot as plt
from math import pi

from fd.analysis.flight_test import FlightTest
from fd.plotting import format_plot
from fd.simulation.simulation import Simulation
from fd.plotting import save_plot


class SimulatedMeasuredComparison:
    simulated_dutch_roll: pd.DataFrame

    def __init__(self, flight_test: FlightTest, simulation: Simulation):
        self.flight_test = flight_test
        self.simulation = simulation

    def run_simulations(self):
        self.simulated_dutch_roll = self.simulation.simulate_dutch_roll(
            self.flight_test.df_dutch_roll
        )
        self.simulated_phugoid = self.simulation.simulate_phugoid(self.flight_test.df_phugoid)
        self.simulated_aperiodic_roll = self.simulation.simulate_aperiodic_roll(
            self.flight_test.df_aperiodic_roll
        )
        self.simulated_dutch_roll_yd = self.simulation.simulate_dutch_roll_yd(
            self.flight_test.df_dutch_roll_yd
        )
        self.simulated_spiral = self.simulation.simulate_spiral(self.flight_test.df_spiral)
        self.simulated_short_period = self.simulation.simulate_short_period(
            self.flight_test.df_short_period
        )

    def plot_responses(self):
        self.plot_phugoid_full()
        self.plot_phugoid()
        self.plot_short_period_full()
        self.plot_short_period()
        self.plot_spiral_full()
        self.plot_spiral()
        self.plot_dutch_roll_full()
        self.plot_dutch_roll()
        self.plot_dutch_roll_yd_full()
        self.plot_dutch_roll_yd()
        self.plot_aperiodic_roll_full()
        self.plot_aperiodic_roll()
        print("Done")

    def plot_dutch_roll(self):
        fig, (ax_p, ax_r) = plt.subplots(2, 1, figsize=(12, 6))

        ax_p.plot(
            self.simulated_dutch_roll.index,
            self.simulated_dutch_roll["p"] * 180 / pi,
            label="Simulated",
        )
        ax_p.plot(
            self.flight_test.df_dutch_roll.index,
            self.flight_test.df_dutch_roll["p"] * 180 / pi,
            label="Measured",
        )
        ax_p.set_ylabel("Roll rate $p$ [°/s]")
        ax_p.legend()

        ax_r.plot(
            self.simulated_dutch_roll.index,
            self.simulated_dutch_roll["r"] * 180 / pi,
        )
        ax_r.plot(
            self.flight_test.df_dutch_roll.index,
            self.flight_test.df_dutch_roll["r"] * 180 / pi,
        )
        ax_r.set_xlabel("Time [s]")
        ax_r.set_ylabel("Yaw rate $r$ [°/s]")

        format_plot()
        save_plot("C:\SVV\Results_final", "dutch_roll")
        plt.show()

    def plot_dutch_roll_full(self):
        fig, (ax_b, ax_phi, ax_p, ax_r) = plt.subplots(4, 1, figsize=(12, 12))

        ax_b.plot(
            self.simulated_dutch_roll.index,
            self.simulated_dutch_roll["beta"] * 180 / pi,
            label="Simulated",
        )
        ax_b.set_ylabel("Sideslip angle $beta$ [°]")

        ax_phi.plot(
            self.simulated_dutch_roll.index,
            self.simulated_dutch_roll["phi"] * 180 / pi,
            label="Simulated",
        )
        ax_phi.plot(
            self.flight_test.df_dutch_roll.index,
            self.flight_test.df_dutch_roll["phi"] * 180 / pi,
            label="Measured",
        )
        ax_phi.set_ylabel("Roll angle $phi$ [°]")
        ax_phi.legend()

        ax_p.plot(
            self.simulated_dutch_roll.index,
            self.simulated_dutch_roll["p"] * 180 / pi,
            label="Simulated",
        )
        ax_p.plot(
            self.flight_test.df_dutch_roll.index,
            self.flight_test.df_dutch_roll["p"] * 180 / pi,
            label="Measured",
        )
        ax_p.set_ylabel("Roll rate $p$ [°/s]")

        ax_r.plot(
            self.simulated_dutch_roll.index,
            self.simulated_dutch_roll["r"] * 180 / pi,
        )
        ax_r.plot(
            self.flight_test.df_dutch_roll.index,
            self.flight_test.df_dutch_roll["r"] * 180 / pi,
        )
        ax_r.set_xlabel("Time [s]")
        ax_r.set_ylabel("Yaw rate $r$ [°/s]")

        format_plot()
        save_plot("C:\SVV\Results_final", "dutch_roll_full")
        plt.show()

    def plot_dutch_roll_yd(self):
        fig, (ax_p, ax_r) = plt.subplots(2, 1, figsize=(12, 6))

        ax_p.plot(
            self.simulated_dutch_roll_yd.index,
            self.simulated_dutch_roll_yd["p"] * 180 / pi,
            label="Simulated",
        )
        ax_p.plot(
            self.flight_test.df_dutch_roll_yd.index,
            self.flight_test.df_dutch_roll_yd["p"] * 180 / pi,
            label="Measured",
        )
        ax_p.set_ylabel("Roll rate $p$ [°/s]")
        ax_p.legend()

        ax_r.plot(
            self.simulated_dutch_roll_yd.index,
            self.simulated_dutch_roll_yd["r"] * 180 / pi,
        )
        ax_r.plot(
            self.flight_test.df_dutch_roll_yd.index,
            self.flight_test.df_dutch_roll_yd["r"] * 180 / pi,
        )
        ax_r.set_xlabel("Time [s]")
        ax_r.set_ylabel("Yaw rate $r$ [°/s]")

        format_plot()
        save_plot("C:\SVV\Results_final", "dutch_roll_yd")
        plt.show()

    def plot_dutch_roll_yd_full(self):
        fig, (ax_b, ax_phi, ax_p, ax_r) = plt.subplots(4, 1, figsize=(12, 12))

        ax_b.plot(
            self.simulated_dutch_roll_yd.index,
            self.simulated_dutch_roll_yd["beta"] * 180 / pi,
            label="Simulated",
        )
        ax_b.set_ylabel("Sideslip angle $beta$ [°]")

        ax_phi.plot(
            self.simulated_dutch_roll_yd.index,
            self.simulated_dutch_roll_yd["phi"] * 180 / pi,
            label="Simulated",
        )
        ax_phi.plot(
            self.flight_test.df_dutch_roll_yd.index,
            self.flight_test.df_dutch_roll_yd["phi"] * 180 / pi,
            label="Measured",
        )
        ax_phi.set_ylabel("Roll angle $phi$ [°]")
        ax_phi.legend()

        ax_p.plot(
            self.simulated_dutch_roll_yd.index,
            self.simulated_dutch_roll_yd["p"] * 180 / pi,
            label="Simulated",
        )
        ax_p.plot(
            self.flight_test.df_dutch_roll_yd.index,
            self.flight_test.df_dutch_roll_yd["p"] * 180 / pi,
            label="Measured",
        )
        ax_p.set_ylabel("Roll rate $p$ [°/s]")

        ax_r.plot(
            self.simulated_dutch_roll_yd.index,
            self.simulated_dutch_roll_yd["r"] * 180 / pi,
        )
        ax_r.plot(
            self.flight_test.df_dutch_roll_yd.index,
            self.flight_test.df_dutch_roll_yd["r"] * 180 / pi,
        )
        ax_r.set_xlabel("Time [s]")
        ax_r.set_ylabel("Yaw rate $r$ [°/s]")

        format_plot()
        save_plot("C:\SVV\Results_final", "dutch_roll_yd_full")
        plt.show()

    def plot_aperiodic_roll(self):
        fig, (ax_p, ax_r) = plt.subplots(2, 1, figsize=(12, 6))

        ax_p.plot(
            self.simulated_aperiodic_roll.index,
            self.simulated_aperiodic_roll["p"] * 180 / pi,
            label="Simulated",
        )
        ax_p.plot(
            self.flight_test.df_aperiodic_roll.index,
            self.flight_test.df_aperiodic_roll["p"] * 180 / pi,
            label="Measured",
        )
        ax_p.set_ylabel("Roll rate $p$ [°/s]")
        ax_p.legend()

        ax_r.plot(
            self.simulated_aperiodic_roll.index,
            self.simulated_aperiodic_roll["r"] * 180 / pi,
        )
        ax_r.plot(
            self.flight_test.df_aperiodic_roll.index,
            self.flight_test.df_aperiodic_roll["r"] * 180 / pi,
        )
        ax_r.set_xlabel("Time [s]")
        ax_r.set_ylabel("Yaw rate $r$ [°/s]")

        format_plot()
        save_plot("C:\SVV\Results_final", "aperiodic_roll")
        plt.show()

    def plot_aperiodic_roll_full(self):
        fig, (ax_b, ax_phi, ax_p, ax_r) = plt.subplots(4, 1, figsize=(12, 12))

        ax_b.plot(
            self.simulated_aperiodic_roll.index,
            self.simulated_aperiodic_roll["beta"] * 180 / pi,
            label="Simulated",
        )
        ax_b.set_ylabel("Sideslip angle $beta$ [°]")

        ax_phi.plot(
            self.simulated_aperiodic_roll.index,
            self.simulated_aperiodic_roll["phi"] * 180 / pi,
            label="Simulated",
        )
        ax_phi.plot(
            self.flight_test.df_aperiodic_roll.index,
            self.flight_test.df_aperiodic_roll["phi"] * 180 / pi,
            label="Measured",
        )
        ax_phi.set_ylabel("Roll angle $phi$ [°]")
        ax_phi.legend()

        ax_p.plot(
            self.simulated_aperiodic_roll.index,
            self.simulated_aperiodic_roll["p"] * 180 / pi,
            label="Simulated",
        )
        ax_p.plot(
            self.flight_test.df_aperiodic_roll.index,
            self.flight_test.df_aperiodic_roll["p"] * 180 / pi,
            label="Measured",
        )
        ax_p.set_ylabel("Roll rate $p$ [°/s]")

        ax_r.plot(
            self.simulated_aperiodic_roll.index,
            self.simulated_aperiodic_roll["r"] * 180 / pi,
        )
        ax_r.plot(
            self.flight_test.df_aperiodic_roll.index,
            self.flight_test.df_aperiodic_roll["r"] * 180 / pi,
        )
        ax_r.set_xlabel("Time [s]")
        ax_r.set_ylabel("Yaw rate $r$ [°/s]")

        format_plot()
        save_plot("C:\SVV\Results_final", "aperiodic_roll_full")
        plt.show()

    def plot_spiral(self):
        fig, (ax_p, ax_r) = plt.subplots(2, 1, figsize=(12, 6))

        ax_p.plot(
            self.simulated_spiral.index,
            self.simulated_spiral["p"] * 180 / pi,
            label="Simulated",
        )
        ax_p.plot(
            self.flight_test.df_spiral.index,
            self.flight_test.df_spiral["p"] * 180 / pi,
            label="Measured",
        )
        ax_p.set_ylabel("Roll rate $p$ [°/s]")
        ax_p.legend()

        ax_r.plot(
            self.simulated_spiral.index,
            self.simulated_spiral["r"] * 180 / pi,
        )
        ax_r.plot(
            self.flight_test.df_spiral.index,
            self.flight_test.df_spiral["r"] * 180 / pi,
        )
        ax_r.set_xlabel("Time [s]")
        ax_r.set_ylabel("Yaw rate $r$ [°/s]")

        format_plot()
        save_plot("C:\SVV\Results_final", "spiral")
        plt.show()

    def plot_spiral_full(self):
        fig, (ax_b, ax_phi, ax_p, ax_r) = plt.subplots(4, 1, figsize=(12, 12))

        ax_b.plot(
            self.simulated_spiral.index, self.simulated_spiral["beta"] * 180 / pi, label="Simulated"
        )
        ax_b.set_ylabel("Sideslip angle $beta$ [°]")

        ax_phi.plot(
            self.simulated_spiral.index, self.simulated_spiral["phi"] * 180 / pi, label="Simulated"
        )
        ax_phi.plot(
            self.flight_test.df_spiral.index,
            self.flight_test.df_spiral["phi"] * 180 / pi,
            label="Measured",
        )
        ax_phi.set_ylabel("Roll angle $phi$ [°]")
        ax_phi.legend()

        ax_p.plot(
            self.simulated_spiral.index,
            self.simulated_spiral["p"] * 180 / pi,
            label="Simulated",
        )
        ax_p.plot(
            self.flight_test.df_spiral.index,
            self.flight_test.df_spiral["p"] * 180 / pi,
            label="Measured",
        )
        ax_p.set_ylabel("Roll rate $p$ [°/s]")

        ax_r.plot(
            self.simulated_spiral.index,
            self.simulated_spiral["r"] * 180 / pi,
        )
        ax_r.plot(
            self.flight_test.df_spiral.index,
            self.flight_test.df_spiral["r"] * 180 / pi,
        )
        ax_r.set_xlabel("Time [s]")
        ax_r.set_ylabel("Yaw rate $r$ [°/s]")

        format_plot()
        save_plot("C:\SVV\Results_final", "spiral_full")
        plt.show()

    def plot_phugoid(self):
        fig, (ax_p, ax_r) = plt.subplots(2, 1, figsize=(12, 6))
        ax_p.plot(
            self.simulated_phugoid.index,
            self.simulated_phugoid["u_hat"] * self.flight_test.df_phugoid["tas"].iloc[0]
            + self.flight_test.df_phugoid["tas"].iloc[0],
            label="Simulated",
        )
        ax_p.plot(
            self.flight_test.df_phugoid.index,
            self.flight_test.df_phugoid["tas"],
            label="Measured",
        )
        ax_p.set_ylabel("True airspeed $V_{TAS}$ [m/s]")
        ax_p.legend()

        ax_r.plot(
            self.simulated_phugoid.index,
            self.simulated_phugoid["q"] * 180 / pi,
        )
        ax_r.plot(
            self.flight_test.df_phugoid.index,
            self.flight_test.df_phugoid["q"] * 180 / pi,
        )
        ax_r.set_xlabel("Time [s]")
        ax_r.set_ylabel("Pitch rate $q$ [°/s]")

        format_plot()
        save_plot("C:\SVV\Results_final", "phugoid")
        plt.show()

    def plot_phugoid_full(self):
        fig, (ax_u, ax_alpha, ax_theta, ax_q) = plt.subplots(4, 1, figsize=(12, 12))
        ax_u.plot(
            self.simulated_phugoid.index,
            self.simulated_phugoid["u_hat"] * self.flight_test.df_phugoid["tas"].iloc[0]
            + self.flight_test.df_phugoid["tas"].iloc[0],
            label="Simulated",
        )
        ax_u.plot(
            self.flight_test.df_phugoid.index,
            self.flight_test.df_phugoid["tas"],
            label="Measured",
        )
        ax_u.set_ylabel("True airspeed $V_{TAS}$ [m/s]")
        ax_u.legend()

        ax_alpha.plot(
            self.simulated_phugoid.index,
            self.simulated_phugoid["alpha"] * 180 / pi,
        )
        ax_alpha.plot(
            self.flight_test.df_phugoid.index,
            self.flight_test.df_phugoid["alpha"] * 180 / pi,
        )
        ax_alpha.set_ylabel("Angle of attack $alpha$ [°]")

        ax_theta.plot(
            self.simulated_phugoid.index,
            self.simulated_phugoid["theta"] * 180 / pi,
        )
        ax_theta.plot(
            self.flight_test.df_phugoid.index,
            self.flight_test.df_phugoid["theta"] * 180 / pi,
        )
        ax_theta.set_ylabel("Pitch angle $theta$ [°]")

        ax_q.plot(
            self.simulated_phugoid.index,
            self.simulated_phugoid["q"] * 180 / pi,
        )
        ax_q.plot(
            self.flight_test.df_phugoid.index,
            self.flight_test.df_phugoid["q"] * 180 / pi,
        )
        ax_q.set_xlabel("Time [s]")
        ax_q.set_ylabel("Pitch rate $q$ [°/s]")

        format_plot()
        save_plot("C:\SVV\Results_final", "phugoid_full")
        plt.show()

    def plot_short_period(self):
        fig, (ax_u, ax_q) = plt.subplots(2, 1, figsize=(12, 6))
        ax_u.plot(
            self.simulated_short_period.index,
            self.simulated_short_period["u_hat"] * self.flight_test.df_phugoid["tas"].iloc[0]
            + self.flight_test.df_phugoid["tas"].iloc[0],
            label="Simulated",
        )
        ax_u.plot(
            self.flight_test.df_short_period.index,
            self.flight_test.df_short_period["tas"],
            label="Measured",
        )
        ax_u.set_ylabel("True airspeed $V_{TAS}$ [m/s]")
        ax_u.legend()

        ax_q.plot(
            self.simulated_short_period.index,
            self.simulated_short_period["q"] * 180 / pi,
        )
        ax_q.plot(
            self.flight_test.df_short_period.index,
            self.flight_test.df_short_period["q"] * 180 / pi,
        )
        ax_q.set_xlabel("Time [s]")
        ax_q.set_ylabel("Pitch rate $q$ [°/s]")

        format_plot()
        save_plot("C:\SVV\Results_final", "short_period")
        plt.show()

    def plot_short_period_full(self):
        fig, (ax_u, ax_alpha, ax_theta, ax_q) = plt.subplots(4, 1, figsize=(12, 12))
        ax_u.plot(
            self.simulated_short_period.index,
            self.simulated_short_period["u_hat"] * self.flight_test.df_phugoid["tas"].iloc[0]
            + self.flight_test.df_phugoid["tas"].iloc[0],
            label="Simulated",
        )
        ax_u.plot(
            self.flight_test.df_short_period.index,
            self.flight_test.df_short_period["tas"],
            label="Measured",
        )
        ax_u.set_ylabel("True airspeed $V_{TAS}$ [m/s]")
        ax_u.legend()

        ax_alpha.plot(
            self.simulated_short_period.index,
            self.simulated_short_period["alpha"] * 180 / pi,
        )
        ax_alpha.plot(
            self.flight_test.df_short_period.index,
            self.flight_test.df_short_period["alpha"] * 180 / pi,
        )
        ax_alpha.set_ylabel("Angle of attack $alpha$ [°]")

        ax_theta.plot(
            self.simulated_short_period.index,
            self.simulated_short_period["theta"] * 180 / pi,
        )
        ax_theta.plot(
            self.flight_test.df_short_period.index,
            self.flight_test.df_short_period["theta"] * 180 / pi,
        )
        ax_theta.set_ylabel("Pitch angle $theta$ [°]")

        ax_q.plot(
            self.simulated_short_period.index,
            self.simulated_short_period["q"] * 180 / pi,
        )
        ax_q.plot(
            self.flight_test.df_short_period.index,
            self.flight_test.df_short_period["q"] * 180 / pi,
        )
        ax_q.set_xlabel("Time [s]")
        ax_q.set_ylabel("Pitch rate $q$ [°/s]")

        format_plot()
        save_plot("C:\SVV\Results_final", "short_period_full")
        plt.show()
