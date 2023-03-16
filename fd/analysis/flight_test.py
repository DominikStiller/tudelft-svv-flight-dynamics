from pathlib import Path

import pandas as pd

from fd.analysis.data_sheet import DataSheet, AveragedDataSheet
from fd.analysis.ftis_measurements import process_ftis_measurements
from fd.analysis.thrust import calculate_thrust_from_df
from fd.io import load_ftis_measurements
from fd.simulation import constants
from fd.structs import AerodynamicParameters


class FlightTest:
    """Stores raw measurements, data sheet and estimated parameters"""

    df: pd.DataFrame
    data_sheet: AveragedDataSheet

    def __init__(self, data_path: str):
        self._load(data_path)
        self._add_derived_timeseries()
        # TODO compare with FTIS, only visually by plotting both in Jupyter notebook

    def _load(self, data_path: str):
        self.df = process_ftis_measurements(load_ftis_measurements(data_path))
        self.data_sheet = AveragedDataSheet(
            {p.name: DataSheet(str(p)) for p in Path(data_path).glob("**/*.xlsx")}
        )

    def _add_derived_timeseries(self):
        self.df["time_min"] = self.df.index / 60
        self.df["m"] = (
            self.data_sheet.mass_initial - self.df["fuel_used_left"] - self.df["fuel_used_right"]
        )
        self.df["W"] = self.df["m"] * constants.g
        self._add_thrust_timeseries()

    def _add_thrust_timeseries(self):
        def _add_thrusts(df):
            # Calculate thrust for stationary measurements
            df[["T_left", "T_right"]] = calculate_thrust_from_df(df)
            # df[["T_left", "T_right"]] = calculate_thrust_from_df_exe(df)
            df["T"] = df["T_left"] + df["T_right"]

        _add_thrusts(self.data_sheet.df_clcd)
        _add_thrusts(self.data_sheet.df_elevator_trim)
        _add_thrusts(self.data_sheet.df_cg_shift)

    def _add_reduced_velocity_timeseries(self):
        pass

    def _add_reduced_delta_timeseries(self):
        pass

    def analyze(self) -> AerodynamicParameters:
        # Estimate aerodynamic and eigenmotion parameters
        # Plot aerodynamic curves
        pass


if __name__ == "__main__":
    test = FlightTest("data/ref_2023")
