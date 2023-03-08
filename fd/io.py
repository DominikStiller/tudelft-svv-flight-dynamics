import pandas as pd
import scipy


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


if __name__ == "__main__":
    import sys

    extract_column_descriptions(sys.argv[1])
