from typing import Any

import numpy as np
import pandas as pd

from fd.analysis.aerodynamic_analysis import calc_mach, calc_static_temp
from fd.conversion import lbs_to_kg, timestamp_to_s, ft_to_m, kts_to_ms, lbshr_to_kgs, C_to_K
from fd.io import load_data_sheet
from fd.simulation.constants import mass_basic_empty
from fd.util import mean_not_none


class DataSheet:
    def __init__(self, data_path: str):
        self._extract(load_data_sheet(data_path))
        self._add_derived_timeseries()

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

    def _add_derived_timeseries(self):
        for df in [self.df_clcd, self.df_elevator_trim, self.df_cg_shift]:
            df["M"] = df.apply(lambda row: calc_mach(row["h"], row["ias"]), axis=1)
            df["T_static"] = df.apply(
                lambda row: calc_static_temp(row["T_total"], row["M"]), axis=1
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
    df["time_min"] = df["time"] / 60
    # Force as float since Excel sheet stores them as string
    df = df.astype("float64")

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


class AveragedDataSheet:
    def __init__(self, data_sheets: dict[str, DataSheet]):
        self.data_sheet_names = list(data_sheets.keys())
        self.data_sheets = list(data_sheets.values())
        self._average()

    def _average(self):
        self.df_clcd = self._average_dataframe_and_check([ds.df_clcd for ds in self.data_sheets])
        self.df_elevator_trim = self._average_dataframe_and_check(
            [ds.df_elevator_trim for ds in self.data_sheets]
        )
        self.df_cg_shift = self._average_dataframe_and_check(
            [ds.df_cg_shift for ds in self.data_sheets]
        )

        self.timestamp_phugoid = mean_not_none([ds.timestamp_phugoid for ds in self.data_sheets])
        self.timestamp_short_period = mean_not_none(
            [ds.timestamp_short_period for ds in self.data_sheets]
        )
        self.timestamp_dutch_roll = mean_not_none(
            [ds.timestamp_dutch_roll for ds in self.data_sheets]
        )
        self.timestamp_dutch_roll_yd = mean_not_none(
            [ds.timestamp_dutch_roll_yd for ds in self.data_sheets]
        )
        self.timestamp_aperiodic_roll = mean_not_none(
            [ds.timestamp_aperiodic_roll for ds in self.data_sheets]
        )
        self.timestamp_spiral = mean_not_none([ds.timestamp_spiral for ds in self.data_sheets])
        self.mass_initial = mean_not_none([ds.mass_initial for ds in self.data_sheets])

    def _average_dataframe_and_check(
        self, dfs: list[pd.DataFrame], threshold_pct=5
    ) -> pd.DataFrame:
        # Calculate mean of non-NA values
        df_mean = dfs[0]
        for ds in dfs[1:]:
            df_mean = df_mean.add(ds, fill_value=0)
        df_mean /= len(dfs)

        for df_idx, df in enumerate(dfs):
            # Calculate mean absolute percentage error
            error: pd.DataFrame = ((df - df_mean) / df_mean).abs() * 100
            exceeds_error_threshold = np.argwhere(error.to_numpy() > threshold_pct)
            if len(exceeds_error_threshold) > 0:
                data_sheet_name = self.data_sheet_names[df_idx]
                print(f"Some values in {data_sheet_name} exceed the {threshold_pct} % threshold")
                for i in range(exceeds_error_threshold.shape[0]):
                    row, col = exceeds_error_threshold[i]
                    column_name = df.columns[col]
                    print(
                        f"Row {row}, column {column_name} ({error.iloc[row, col]:.3} %): "
                        f"{data_sheet_name} = {df.iloc[row, col]:.3}, mean = {df_mean.iloc[row, col]:.3}"
                    )
                print()

        return df_mean
