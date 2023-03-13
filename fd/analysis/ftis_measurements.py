import pandas as pd

from fd.conversion import lbshr_to_kgs, lbs_to_kg, ft_to_m, kts_to_ms, C_to_K

COLUMNS = {
    "vane_AOA": "alpha",
    "elevator_dte": "delta_t_e",
    "lh_engine_FMF": "fuel_flow_left",
    "rh_engine_FMF": "fuel_flow_right",
    "lh_engine_FU": "fuel_used_left",
    "rh_engine_FU": "fuel_used_right",
    "column_Se": "se",
    "delta_a": "delta_a",
    "delta_e": "delta_e",
    "delta_r": "delta_r",
    "Dadc1_bcAlt": "h",
    "Dadc1_mach": "M",
    "Dadc1_cas": "cas",
    "Dadc1_tas": "tas",
    "Dadc1_sat": "T_static",
    "Dadc1_tat": "T_total",
}


def process_ftis_measurements(df: pd.DataFrame):
    df = df[COLUMNS.keys()].rename(columns=COLUMNS)

    # Convert to SI units
    df["fuel_flow_left"] = lbshr_to_kgs(df["fuel_flow_left"])
    df["fuel_flow_right"] = lbshr_to_kgs(df["fuel_flow_right"])
    df["fuel_used_left"] = lbs_to_kg(df["fuel_used_left"])
    df["fuel_used_right"] = lbs_to_kg(df["fuel_used_right"])
    df["h"] = ft_to_m(df["h"])
    df["cas"] = kts_to_ms(df["cas"])
    df["tas"] = kts_to_ms(df["tas"])
    df["T_static"] = C_to_K(df["T_static"])
    df["T_total"] = C_to_K(df["T_total"])

    return df
