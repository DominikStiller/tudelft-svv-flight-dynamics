from typing import Any

import pandas as pd

from fd.conversion import lbs_to_kg, timestamp_to_s, ft_to_m, kts_to_ms, lbshr_to_kgs, C_to_K
from fd.io import load_data_sheet
from fd.simulation.constants import mass_basic_empty


class DataSheet:
    def __init__(self, data_path: str):
        self._extract(load_data_sheet(data_path))

    def _extract(self, ws: list[list[Any]]):
        """
        Extract parameters from PFDS into variables.

        Args:
            ws: worksheet from Excel file as list of lists
        """
        self._extract_mass(ws)

        self.df_clcd = extract_single_timeseries(ws, 27, 33)
        self.df_elevator_trim = extract_single_timeseries(ws, 58, 64)
        self.df_cg_shift = extract_single_timeseries(ws, 74, 75)

        self.timestamp_phugoid = timestamp_to_s(ws[82][3])
        self.timestamp_short_period = timestamp_to_s(ws[83][3])
        self.timestamp_dutch_roll = timestamp_to_s(ws[82][6])
        self.timestamp_dutch_roll_yd = timestamp_to_s(ws[83][6])
        self.timestamp_aperiodic_roll = timestamp_to_s(ws[82][9])
        self.timestamp_spiral = timestamp_to_s(ws[82][9])

    def _extract_mass(self, ws: list[list[Any]]):
        self.mass_pilot_1 = ws[7][7]
        self.mass_pilot_2 = ws[8][7]
        self.mass_coordinator = ws[9][7]
        self.mass_observer_1l = ws[10][7]
        self.mass_observer_1r = ws[11][7]
        self.mass_observer_2l = ws[12][7]
        self.mass_observer_2r = ws[13][7]
        self.mass_observer_3l = ws[14][7]
        self.mass_observer_3r = ws[15][7]
        self.mass_block_fuel = lbs_to_kg(ws[17][3])

        self.mass_initial = (
            mass_basic_empty
            + self.mass_block_fuel
            + self.mass_pilot_1
            + self.mass_pilot_2
            + self.mass_coordinator
            + self.mass_observer_1l
            + self.mass_observer_1r
            + self.mass_observer_2l
            + self.mass_observer_2r
            + self.mass_observer_3l
            + self.mass_observer_3r
        )


def extract_single_timeseries(ws: list[list[Any]], row_start: int, row_end: int) -> pd.DataFrame:
    if ws[row_start - 1][1] is None:
        column_names = ws[row_start - 3]
    else:
        # Empty row between header and data missing for cg shift data
        column_names = ws[row_start - 2]
    rows = ws[row_start : row_end + 1]

    df = pd.DataFrame(rows, columns=column_names).drop(
        columns=["nr.", "ET*", None], errors="ignore"
    )
    df = df.dropna(subset="time").reset_index(drop=True)
    df["time"] = df["time"].apply(timestamp_to_s)
    # Force as float since Excel sheet stores them as string
    df = df.astype("float64")
    df = df.set_index("time", drop=True)

    df = df.rename(
        columns={
            "hp": "h",
            "IAS": "ias",
            "a": "alpha",
            "de": "delta_e",
            "detr": "delta_e_t",
            "Fe": "F_e",
            "FFl": "fuel_flow_left",
            "FFr": "fuel_flow_right",
            "F. used": "fuel_used",
            "TAT": "T_total",
        }
    )

    df["h"] = ft_to_m(df["h"])
    df["ias"] = kts_to_ms(df["ias"])
    df["fuel_flow_left"] = lbshr_to_kgs(df["fuel_flow_left"])
    df["fuel_flow_right"] = lbshr_to_kgs(df["fuel_flow_right"])
    df["fuel_used"] = lbs_to_kg(df["fuel_used"])
    df["T_total"] = C_to_K(df["T_total"])

    return df
