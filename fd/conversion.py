import datetime
import math
from typing import Union


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


def timestamp_to_s(timestamp: Union[str, datetime.time]):
    """Convert datetime or [h.]mm[:ss] string (both from Excel sheet) to seconds"""
    # TODO check if this format is consistent across data sheets
    if not timestamp:
        return None

    if isinstance(timestamp, str):
        if "." in timestamp:
            # Format is h.mm:ss
            hour, minute_second = timestamp.split(".")
            hour = float(hour)
        else:
            # Format is mm:ss
            minute_second = timestamp
            hour = 0

        minute_second = minute_second.split(":")[:2]
        if len(minute_second) == 2:
            minute, second = minute_second
            second = float(second)
        else:
            # Did not include ss part
            minute = minute_second[0]
            second = 0
        minute = float(minute)

        assert minute < 60
        assert second < 60
    elif isinstance(timestamp, datetime.time):
        hour = 0
        # This is not a mistake, the datetime is interpreted incorrectly
        minute = timestamp.hour
        second = timestamp.minute
    else:
        raise "Unsupported time type"

    return hour * 3600 + minute * 60 + second
