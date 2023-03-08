import pandas as pd

from fd.analysis.aerodynamic_plots import plot_cl_alpha
from fd.analysis.data_sheet import DataSheet
from fd.analysis.ftis_measurements import FTISMeasurements
from fd.structs import AerodynamicParameters


class FlightTest:
    measurements: FTISMeasurements
    data_sheet: DataSheet

    """Stores raw measurements, data sheet and estimated parameters"""

    def __init__(self, data_path: str):
        self._load(data_path)

    def _load(self, data_path: str):
        # Load measurements file and data sheet
        # Add time-dependent columns
        pass

    def _add_mass_timeseries(self):
        pass

    def _add_thrust_timeseries(self):
        pass

    def _add_reduced_velocity_timeseries(self):
        pass

    def _add_reduced_delta_timeseries(self):
        pass

    def analyze(self) -> AerodynamicParameters:
        # Estimate aerodynamic and eigenmotion parameters
        # Plot aerodynamic curves
        plot_cl_alpha(self)
        pass

    def get_timeseries(self, column_name: str) -> pd.Series:
        pass

    @property
    def get_clcd_measurements(self) -> pd.DataFrame:
        # TODO return data from sheet or FTIS?
        return self.measurements.df.iloc[self.data_sheet.clcd_timestamps]

    @property
    def get_elevator_trim_timestamps(self) -> pd.DataFrame:
        return []

    @property
    def get_cg_shift_timestamps(self) -> pd.DataFrame:
        return []
