from fd.analysis.flight_test import FlightTest
from fd.simulation.aircraft_model import AircraftModel
from fd.simulation.simulation import Simulation
from fd.structs import AerodynamicParameters
from fd.validation.comparison import SimulatedMeasuredComparison
import control.matlab as ml
import numpy as np
import matplotlib.pyplot as plt
from fd.plotting import *

test = "pulse_aileron"

if test == "pulse_elevator":
    flight_test = FlightTest("data/B24")
    aero_params = AerodynamicParameters(
        C_L_alpha=4.758556374647304,
        alpha_0=-0.02312478307006348,
        C_D_0=0.023439123324849084,
        C_m_alpha=-0.5554065208385275,
        C_m_delta=-1.3380975545274032,
        e=1.0713238368125688,
    )

    aircraft_model = AircraftModel(aero_params)
    A, B, C, D = aircraft_model.get_state_space_matrices_symmetric(4500, 150, 0.8, 0)
    # print(np.linalg.eig(A)[0])
    sys = ml.ss(A, B, C, D)
    t = np.linspace(0, 10, 10000)
    x0 = [[0], [0], [0], [0]]
    u = np.zeros([len(t), 1])
    u[0] = 0.1
    yout, t, xout = ml.lsim(sys, u, t, x0)
    plt.figure(figsize=(12, 3))
    plt.plot(t, xout[:, 1])
    plt.ylabel("$alpha$ [rad]")
    plt.xlabel("Time [s]")
    format_plot()
    save_plot("data/", "int_test_pulse_elev")
    plt.show()
    # aircraft_model.get_response_plots_symmetric(sys, x0, t, u, 150)

elif test == "step_elevator":
    flight_test = FlightTest("data/B24")
    aero_params = AerodynamicParameters(
        C_L_alpha=4.758556374647304,
        alpha_0=-0.02312478307006348,
        C_D_0=0.023439123324849084,
        C_m_alpha=-0.5554065208385275,
        C_m_delta=-1.3380975545274032,
        e=1.0713238368125688,
    )

    aircraft_model = AircraftModel(aero_params)
    A, B, C, D = aircraft_model.get_state_space_matrices_symmetric(4500, 150, 0.8, 0)
    sys = ml.ss(A, B, C, D)
    t = np.linspace(0, 400, 10000)
    x0 = [[0], [0], [0], [0]]
    u = np.ones([len(t), 1])
    u = u * 0.01
    yout, t, xout = ml.lsim(sys, u, t, x0)
    plt.figure(figsize=(6, 3))
    plt.plot(t, xout[:, 3])
    plt.ylabel("$q$ [rad]")
    plt.xlabel("Time [s]")
    format_plot()
    save_plot("data/", "int_test_step_elev_q")
    plt.show()
    # aircraft_model.get_response_plots_symmetric(sys, x0, t, u, 150)

elif test == "pulse_rudder":
    flight_test = FlightTest("data/B24")
    aero_params = AerodynamicParameters(
        C_L_alpha=4.758556374647304,
        alpha_0=-0.02312478307006348,
        C_D_0=0.023439123324849084,
        C_m_alpha=-0.5554065208385275,
        C_m_delta=-1.3380975545274032,
        e=1.0713238368125688,
    )

    aircraft_model = AircraftModel(aero_params)
    A, B, C, D = aircraft_model.get_state_space_matrices_asymmetric(4500, 150, 0.8, 0, 0.8)
    # print(np.linalg.eig(A)[0])
    sys = ml.ss(A, B, C, D)
    t = np.linspace(0, 10, 1000)
    x0 = [[0], [0], [0], [0]]
    u = np.zeros([len(t), 2])
    inp = np.ones([10, 1])
    u[0, 1] = 0.01
    # u[1,1] = -0.01
    # u[0:10, 1:] = inp*0.1
    # u[10:20, 1:] = inp * -0.1
    yout, t, xout = ml.lsim(sys, u, t, x0)
    plt.figure(figsize=(12, 3))
    plt.plot(t, xout[:, 3])
    plt.ylabel("$r$ [rad/sec]")
    plt.xlabel("Time [s]")
    format_plot()
    save_plot("data/", "int_test_pulse_rudder")
    plt.show()
    # aircraft_model.get_response_plots_asymmetric(sys, x0, t, u, 150)

elif test == "pulse_aileron":
    flight_test = FlightTest("data/B24")
    aero_params = AerodynamicParameters(
        C_L_alpha=4.758556374647304,
        alpha_0=-0.02312478307006348,
        C_D_0=0.023439123324849084,
        C_m_alpha=-0.5554065208385275,
        C_m_delta=-1.3380975545274032,
        e=1.0713238368125688,
    )

    aircraft_model = AircraftModel(aero_params)
    A, B, C, D = aircraft_model.get_state_space_matrices_asymmetric(4500, 150, 0.8, 0, 0.8)
    # print(np.linalg.eig(A)[0])
    sys = ml.ss(A, B, C, D)
    t = np.linspace(0, 10, 1000)
    x0 = [[0], [0], [0], [0]]
    u = np.zeros([len(t), 2])
    inp = np.ones([10, 1])
    u[0, 0] = 0.01
    # u[1,1] = -0.01
    # u[0:10, 1:] = inp*0.1
    # u[10:20, 1:] = inp * -0.1
    yout, t, xout = ml.lsim(sys, u, t, x0)
    plt.figure(figsize=(12, 3))
    plt.plot(t, xout[:, 2])
    plt.ylabel("$p$ [rad/sec]")
    plt.xlabel("Time [s]")
    format_plot()
    # save_plot("data/", "int_test_pulse_aileron")
    plt.show()
    # aircraft_model.get_response_plots_asymmetric(sys, x0, t, u, 150)

elif test == "step_aileron":
    flight_test = FlightTest("data/B24")
    aero_params = AerodynamicParameters(
        C_L_alpha=4.758556374647304,
        alpha_0=-0.02312478307006348,
        C_D_0=0.023439123324849084,
        C_m_alpha=-0.5554065208385275,
        C_m_delta=-1.3380975545274032,
        e=1.0713238368125688,
    )

    aircraft_model = AircraftModel(aero_params)
    A, B, C, D = aircraft_model.get_state_space_matrices_asymmetric(4500, 150, 0.8, 0, 0.8)
    # print(np.linalg.eig(A)[0])
    sys = ml.ss(A, B, C, D)
    t = np.linspace(0, 100, 1000)
    x0 = [[0], [0], [0], [0]]
    u = np.zeros([len(t), 2])
    inp = np.ones([len(t), 1])
    u[:, :1] = 0.01 * inp
    # print(u)
    # u[1,1] = -0.01
    # u[0:10, 1:] = inp*0.1
    # u[10:20, 1:] = inp * -0.1
    yout, t, xout = ml.lsim(sys, u, t, x0)
    plt.figure(figsize=(12, 3))
    plt.plot(t, xout[:, 2])
    plt.ylabel("$p$ [rad/sec]")
    plt.xlabel("Time [s]")
    format_plot()
    save_plot("data/", "int_test_step_aileron")
    plt.show()
    # aircraft_model.get_response_plots_asymmetric(sys, x0, t, u, 150)
