from pathlib import Path

import pandas as pd

from fd.analysis.data_sheet import DataSheet, AveragedDataSheet
from fd.analysis.ftis_measurements import FTISMeasurements
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
        # TODO compare with FTIS, only visually by plotting both in Jupyter notebook

    def _add_reduced_velocity_timeseries(self):
        pass

    def _add_reduced_delta_timeseries(self):
        pass

    def analyze(self) -> AerodynamicParameters:
        # Estimate aerodynamic and eigenmotion parameters
        # Plot aerodynamic curves
        pass

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
        return self.ftis_measurements.df[self.data_sheet.timestamp_phugoid]


if __name__ == "__main__":
    test = FlightTest("data/ref_2023")
    # test = FlightTest("data/B24")
    print(test)
