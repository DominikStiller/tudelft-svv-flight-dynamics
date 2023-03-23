import datetime
import math
import re
from typing import Union

import numpy as np


def deg_to_rad(deg):
    """Convert value from degree to radians"""
    return math.pi / 180 * deg


def lbshr_to_kgs(lbshr):
    """Convert value from lbs/hr to kg/s"""
    return 0.45359237 / 3600 * lbshr


def psi_to_Pa(psi):
    """Convert value from psi to Pa"""
    return 6894.757 * psi


def ftmin_to_ms(ftmin):
    """Convert value from ft/min to m/s"""
    return 0.3048 / 60 * ftmin


def lbs_to_kg(lbs):
    """Convert mass in pounds to kilograms"""
    return 0.45359237 * lbs


def kts_to_ms(kts):
    """Convert speed in knots to meters per second"""
    return 1852 / 3600 * kts


def ft_to_m(ft):
    """Convert distance in feet to meters"""
    return 0.3048 * ft


def C_to_K(C):
    """Convert temperature in Celcius to Kelvin"""
    return 273.15 + C


def deg_to_rad(deg):
    """Convert angle in degrees to radians"""
    return np.deg2rad(deg)


def degs_to_rads(degs):
    """Convert rotation in degrees/s to radians/s"""
    return np.deg2rad(degs)


def timestamp_to_s(timestamp: Union[str, datetime.time]):
    """
    Convert datetime.time or timestamp string (both from Excel sheet) to seconds.
    There is no consistent format:
     - datetime.time: mm:ss timestamp is parsed but incorrectly as hh:mm
     - h.mm:ss
     - h.mm
     - h:mm:ss
     - mm:ss
     - mm
    """
    if not timestamp:
        return None

    if isinstance(timestamp, str):
        hour = 0
        minute = 0
        second = 0

        timestamp = timestamp.strip()

        if match := re.fullmatch(r"(\d+)\.(\d+):(\d+)", timestamp):
            # h.mm:ss
            hour, minute, second = match.groups()
        elif match := re.fullmatch(r"(\d+)\.(\d+)", timestamp):
            # h.mm
            hour, minute = match.groups()
        elif match := re.fullmatch(r"(\d+):(\d+):(\d+)", timestamp):
            # h:mm:ss
            hour, minute, second = match.groups()
        elif match := re.fullmatch(r"(\d+):(\d+)", timestamp):
            # mm:ss
            minute, second = match.groups()
        elif match := re.fullmatch(r"(\d+)", timestamp):
            # mm
            minute = match.group(0)
        else:
            raise "Invalid time format"
    elif isinstance(timestamp, datetime.time):
        hour = 0
        # This is not a mistake, the datetime is interpreted incorrectly
        minute = timestamp.hour
        second = timestamp.minute
    else:
        raise "Unsupported time type"

    hour = float(hour)
    minute = float(minute)
    second = float(second)

    assert 0 <= hour
    assert 0 <= minute < 60
    assert 0 <= second < 60

    return hour * 3600 + minute * 60 + second


def inchpound_to_kgm(inchpound):
    """

    Args:
        inchpound (float): Massmoment expressed in inchpounds

    Returns (float): Massmoment expressed in kilogram meters

    """
    return inchpound * 0.45359237 * 0.0254
