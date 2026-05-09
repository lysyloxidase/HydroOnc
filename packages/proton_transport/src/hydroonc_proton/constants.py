"""Reference constants for HydroOnc proton transport."""

from dataclasses import dataclass


# Experimental proton diffusion in water at room temperature.
D_H_PLUS_WATER_CM2_S = 9.31e-5
D_H_PLUS_WATER_A2_PS = D_H_PLUS_WATER_CM2_S * 1.0e4

# Useful context values from the research report.
WATER_SELF_DIFFUSION_CM2_S = 2.3e-5
DEFAULT_HOP_RATE_PER_PS = 1.0
DEFAULT_WATER_WIRE_SPACING_A = 2.65
DEFAULT_STRUCTURAL_DIFFUSION_FRACTION = 0.90
DEFAULT_HOP_CORRELATION_FACTOR = 0.884

TUMOR_EXTRACELLULAR_PH_RANGE = (6.4, 7.0)
TUMOR_INTRACELLULAR_PH_RANGE = (7.0, 7.2)
NORMAL_EXTRACELLULAR_PH_RANGE = (7.2, 7.5)
NORMAL_INTRACELLULAR_PH_RANGE = (6.8, 7.2)


@dataclass(frozen=True)
class ResiduePKA:
    """Intrinsic pKa reference for a titratable residue."""

    residue: str
    pka: float
    acidic: bool


RESIDUE_PKA = {
    "ASP": ResiduePKA("ASP", 3.9, True),
    "GLU": ResiduePKA("GLU", 4.2, True),
    "HIS": ResiduePKA("HIS", 6.45, False),
    "CYS": ResiduePKA("CYS", 8.3, True),
    "TYR": ResiduePKA("TYR", 10.1, True),
    "LYS": ResiduePKA("LYS", 10.5, False),
    "ARG": ResiduePKA("ARG", 12.5, False),
    "NTERM": ResiduePKA("NTERM", 8.0, False),
    "CTERM": ResiduePKA("CTERM", 3.6, True),
}
