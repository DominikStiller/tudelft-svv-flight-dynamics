import warnings
from typing import Any

import pandas as pd
import scipy
from openpyxl.reader.excel import load_workbook


def load_measurements(path: str) -> pd.DataFrame:
    raw = scipy.io.loadmat(f"{path}/measurements.mat", simplify_cells=True)["flightdata"]
    data = {}
    for column_name, values in raw.items():
        data[column_name] = values["data"]
    data = pd.DataFrame(data).set_index("time")
    return data


def extract_column_descriptions(path: str):
    raw = scipy.io.loadmat(f"{path}/measurements.mat", simplify_cells=True)["flightdata"]
    metadata = []
    for column_name, values in raw.items():
        metadata.append(
            {"Column": column_name, "Description": values["description"], "Units": values["units"]}
        )
    metadata = pd.DataFrame(metadata)
    metadata.to_excel("data/column_descriptions.xlsx", index=False)


def load_data_sheet(path: str) -> list[list[Any]]:
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UserWarning, module="openpyxl")
        wb = load_workbook(filename=f"{path}/sheet.xlsx")
    ws = wb.worksheets[0]
    return [[cell.value for cell in row] for row in ws.rows]


if __name__ == "__main__":
    import sys

    extract_column_descriptions(sys.argv[1])
