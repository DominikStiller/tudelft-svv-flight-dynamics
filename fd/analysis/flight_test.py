from pathlib import Path

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

    def __init__(self, data_path: str):
        self.data_sheet = AveragedDataSheet(
            {p.name: DataSheet(str(p)) for p in Path(data_path).glob("**/*.xlsx")}
        )
        self.ftis_measurements = FTISMeasurements(data_path, self.data_sheet.mass_initial)

    def analyze(self) -> AerodynamicParameters:
        # Estimate aerodynamic and eigenmotion parameters
        # Plot aerodynamic curves
        pass

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
