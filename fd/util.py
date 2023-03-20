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


def mean_not_none(l: list[float]) -> float:
    """
    Calculate the mean of all non-None values in x.

    Args:
        l: List of elements

    Returns:
        Mean
    """
    return mean(filter(lambda e: e is not None, l))


def mean_not_nan_df(dfs: list[pd.DataFrame]) -> pd.DataFrame:
    """
    Calculate the cell-wise mean of all non-NAN values in dfs.

    Args:
        dfs: List of DataFrames

    Returns:
        Means as DataFrame
    """
    df_mean = pd.concat(dfs).groupby(level=0).mean()
    return df_mean.astype(dfs[0].dtypes)
