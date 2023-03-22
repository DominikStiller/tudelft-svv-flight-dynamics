from pathlib import Path

from fd.analysis.aerodynamic_plots import plot_cl_alpha, plot_cl_cd
from fd.analysis.aerodynamics import (
    estimate_CL_alpha,
    estimate_CD0_e,
    estimate_Cmalpha,
    calc_Cmdelta,
)
from fd.analysis.data_sheet import DataSheet, AveragedDataSheet
from fd.analysis.ftis_measurements import FTISMeasurements
from fd.simulation.constants import (
    duration_phugoid,
    duration_dutch_roll,
    duration_short_period,
    duration_dutch_roll_yd,
    duration_aperiodic_roll,
    duration_spiral,
)
from fd.structs import AerodynamicParameters


class FlightTest:
    """Stores raw measurements, data sheet and estimated parameters"""

    ftis_measurements: FTISMeasurements
    data_sheet: AveragedDataSheet
    aerodynamic_parameters: AerodynamicParameters

    def __init__(self, data_path: str):
        self.data_sheet = AveragedDataSheet(
            {p.name: DataSheet(str(p)) for p in Path(data_path).glob("**/*.xlsx")}
        )
        self.ftis_measurements = FTISMeasurements(data_path, self.data_sheet.mass_initial)
        self._estimate_aerodynamic_parameters()
        self.data_sheet.add_reduced_elevator_deflection_timeseries(
            self.aerodynamic_parameters.C_m_delta
        )

    def _estimate_aerodynamic_parameters(self):
        C_L_alpha, _, alpha_0 = estimate_CL_alpha(self.df_clcd["C_L"], self.df_clcd["alpha"])
        C_D_0, e = estimate_CD0_e(self.df_clcd["C_D"], self.df_clcd["C_L"])

        # TODO check which parts should use reduced qties
        cg_aft = self.df_cg_shift.iloc[0]
        cg_front = self.df_cg_shift.iloc[1]
        C_m_delta = calc_Cmdelta(
            cg_aft["x_cg"],
            cg_front["x_cg"],
            cg_aft["delta_e"],
            cg_front["delta_e"],
            cg_aft["W"],
            cg_front["W"],
            self.df_cg_shift["tas"].mean(),
            self.df_cg_shift["rho"].mean(),
        )

        C_m_alpha = estimate_Cmalpha(
            self.df_elevator_trim["alpha"], self.df_elevator_trim["delta_e"], C_m_delta
        )

        self.aerodynamic_parameters = AerodynamicParameters(
            C_L_alpha, alpha_0, C_D_0, C_m_alpha, C_m_delta, e
        )

    def make_aerodynamic_plots(self):
        plot_cl_alpha(
            self.df_clcd["C_L"],
            self.df_clcd["alpha"],
            self.aerodynamic_parameters.C_L_alpha,
            self.aerodynamic_parameters.alpha_0,
        )

        plot_cl_cd(
            self.df_clcd["C_L"],
            self.df_clcd["C_D"],
            self.aerodynamic_parameters.C_D_0,
            self.aerodynamic_parameters.e,
        )

    @property
    def df_ftis(self):
        return self.ftis_measurements.df

    @property
    def df_clcd(self):
        return self.data_sheet.df_clcd

    @property
    def df_elevator_trim(self):
        return self.data_sheet.df_elevator_trim

    @property
    def df_cg_shift(self):
        return self.data_sheet.df_cg_shift

    @property
    def df_phugoid(self):
        return self._get_maneuver_df(self.data_sheet.timestamp_phugoid, duration_phugoid)

    @property
    def df_short_period(self):
        return self._get_maneuver_df(self.data_sheet.timestamp_short_period, duration_short_period)

    @property
    def df_dutch_roll(self):
        return self._get_maneuver_df(self.data_sheet.timestamp_dutch_roll, duration_dutch_roll)

    @property
    def df_dutch_roll_yd(self):
        return self._get_maneuver_df(
            self.data_sheet.timestamp_dutch_roll_yd, duration_dutch_roll_yd
        )

    @property
    def df_aperiodic_roll(self):
        return self._get_maneuver_df(
            self.data_sheet.timestamp_aperiodic_roll, duration_aperiodic_roll
        )

    @property
    def df_spiral(self):
        return self._get_maneuver_df(self.data_sheet.timestamp_spiral, duration_spiral)

    def _get_maneuver_df(self, timestamp: float, duration: float):
        """Get rows from FTIS data corresponding to a certain maneuver, identified by start timestamp and duration"""
        # Give 1 second lead time and select by duration
        df = self.df_ftis.loc[timestamp - 1 : timestamp + duration].copy()
        df.index -= df.index[0]
        df["time_min"] -= df["time_min"][0]
        return df
