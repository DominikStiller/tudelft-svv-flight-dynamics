def get_closest(df, time):
    """
    Gets the row in df that is closest before time given. df's index should be timestamps.

    Args:
        time: time
        df: DataFrame

    Returns:
        Row corresponding to the closest time
    """
    closest_idx = df.index.searchsorted(time)
    return df.iloc[closest_idx]
