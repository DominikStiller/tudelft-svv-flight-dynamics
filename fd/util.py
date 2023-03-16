from statistics import mean

import pandas as pd


def get_closest(df: pd.DataFrame, time) -> pd.DataFrame:
    """
    Gets the rows in df that is closest after the time given. If time is past the
    last timestamp, return the last row.

    df's index should be float timestamps.

    Args:
        time: single or multiple timestamps
        df: DataFrame

    Returns:
        Rows corresponding to the closest next time
    """
    closest_idx = df.index.searchsorted(time)
    closest_idx = closest_idx.clip(0, len(df.index) - 1)
    return df.iloc[closest_idx]


def mean_not_none(x: list[float]) -> float:
    """
    Calculate the mean of all non-None values in x.

    Args:
        x: List of elements

    Returns:
        Mean
    """
    return mean(filter(lambda e: e is not None, x))
