import pandas as pd

from fd.analysis.thermodynamics import (
    calc_mach,
    calc_static_temperature,
    calc_static_pressure,
    calc_density,
)
from fd.simulation import constants


def add_common_derived_timeseries(df: pd.DataFrame) -> pd.DataFrame:
    df["W"] = df["m"] * constants.g

    df["M"] = df.apply(lambda row: calc_mach(row["h"], row["cas"]), axis=1)
    df["T_static"] = df.apply(lambda row: calc_static_temperature(row["T_total"], row["M"]), axis=1)
    df["p"] = df.apply(lambda row: calc_static_pressure(row["h"]), axis=1)
    df["rho"] = df.apply(lambda row: calc_density(row["p"], row["T_static"]), axis=1)

    return df
