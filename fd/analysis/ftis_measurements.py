from fd.analysis.util import add_common_derived_timeseries
from fd.conversion import lbshr_to_kgs, lbs_to_kg, ft_to_m, kts_to_ms, C_to_K
from fd.io import load_ftis_measurements
from fd.simulation import constants

COLUMNS = {
    "vane_AOA": "alpha",
    "elevator_dte": "delta_t_e",
    "lh_engine_FMF": "fuel_flow_left",
    "rh_engine_FMF": "fuel_flow_right",
    "lh_engine_FU": "fuel_used_left",
    "rh_engine_FU": "fuel_used_right",
    "column_Se": "s_e",
    "column_fe": "F_e",
    "delta_a": "delta_a",
    "delta_e": "delta_e",
    "delta_r": "delta_r",
    "Dadc1_bcAlt": "h",
    "Dadc1_mach": "M",
    "Dadc1_cas": "cas",
    "Dadc1_tas": "tas",
    "Dadc1_sat": "T_static",
    "Dadc1_tat": "T_total",
    "Ahrs1_Roll": "phi",
    "Ahrs1_Pitch": "theta",
    "Ahrs1_bRollRate": "p",
    "Ahrs1_bPitchRate": "q",
    "Ahrs1_bYawRate": "r",
    "Fms1_trueHeading": "chi",
}


class FTISMeasurements:
    def __init__(self, data_path: str, mass_initial: float):
        self.df = load_ftis_measurements(data_path)
        self._process_ftis_measurements()
        self._add_derived_timeseries(mass_initial)

    def _process_ftis_measurements(self):
        self.df = self.df[COLUMNS.keys()].rename(columns=COLUMNS)

        # Convert to SI units
        self.df["fuel_flow_left"] = lbshr_to_kgs(self.df["fuel_flow_left"])
        self.df["fuel_flow_right"] = lbshr_to_kgs(self.df["fuel_flow_right"])
        self.df["fuel_used_left"] = lbs_to_kg(self.df["fuel_used_left"])
        self.df["fuel_used_right"] = lbs_to_kg(self.df["fuel_used_right"])
        self.df["h"] = ft_to_m(self.df["h"])
        self.df["cas"] = kts_to_ms(self.df["cas"])
        self.df["tas"] = kts_to_ms(self.df["tas"])
        self.df["T_static"] = C_to_K(self.df["T_static"])
        self.df["T_total"] = C_to_K(self.df["T_total"])

        return self.df

    def _add_derived_timeseries(self, mass_initial: float):
        self.df["time_min"] = self.df.index / 60
        self.df["m"] = mass_initial - self.df["fuel_used_left"] - self.df["fuel_used_right"]
        self.df = add_common_derived_timeseries(self.df)
