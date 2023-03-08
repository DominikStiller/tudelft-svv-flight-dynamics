class DataSheet:
    def __init__(self, data_path: str):
        # Extract mass and distribution
        # Extract stationary measurements
        # Extract timestamps of maneuvers
        pass

    @property
    def clcd_timestamps(self) -> list[float]:
        return []

    @property
    def elevator_trim_timestamps(self) -> list[float]:
        return []

    @property
    def cg_shift_timestamps(self) -> list[float]:
        return []

    @property
    def phugoid_timestamp(self) -> float:
        return 0

    @property
    def short_period_timestamp(self) -> float:
        return 0

    @property
    def dutch_roll_timestamp(self) -> float:
        return 0

    @property
    def dutch_roll_yd_timestamp(self) -> float:
        return 0

    @property
    def aperiodic_roll_timestamp(self) -> float:
        return 0

    @property
    def spiral_timestamp(self) -> float:
        return 0
