from dataclasses import dataclass


@dataclass()
class SimulationOutput:
    pass


@dataclass
class AerodynamicParameters:
    C_L_alpha: float
    C_D_0: float
    C_m_alpha: float
    C_m_delta: float
    e: float
