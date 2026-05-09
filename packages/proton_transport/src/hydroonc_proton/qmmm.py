"""QM/MM and AIMD configuration helpers for proton transport studies."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Optional

import numpy as np


@dataclass(frozen=True)
class QMMMRegion:
    """Definition of a QM region embedded in a classical environment."""

    qm_atom_indices: tuple[int, ...]
    link_atom_indices: tuple[int, ...] = ()
    charge: int = 1
    spin_multiplicity: int = 1


@dataclass(frozen=True)
class QMMMConfig:
    """Configuration for electrostatically embedded QM/MM proton transport."""

    region: QMMMRegion
    qm_method: str = "BLYP-D3"
    qm_basis: str = "DZVP-MOLOPT-SR-GTH"
    mm_forcefield: str = "AMBER ff19SB + OPC"
    embedding: str = "electrostatic"
    timestep_fs: float = 0.5
    temperature_K: float = 300.0


@dataclass(frozen=True)
class AIMDConfig:
    """CP2K AIMD validation settings."""

    functional: str = "BLYP"
    dispersion: str = "D3"
    timestep_fs: float = 0.5
    ensemble: str = "NVT"
    temperature_K: float = 300.0
    steps: int = 2000


class QMMMProtonTransport:
    """Build QM/MM regions and CP2K inputs for proton-transfer validation."""

    def select_water_wire_region(
        self,
        oxygen_positions_A: np.ndarray,
        proton_position_A: np.ndarray,
        radius_A: float = 4.0,
        max_waters: Optional[int] = None,
    ) -> QMMMRegion:
        """Select waters near an excess proton for a QM region."""

        oxygen_positions = np.asarray(oxygen_positions_A, dtype=float)
        proton_position = np.asarray(proton_position_A, dtype=float)
        if oxygen_positions.ndim != 2 or oxygen_positions.shape[1] != 3:
            raise ValueError("oxygen_positions_A must have shape (n, 3)")
        if proton_position.shape != (3,):
            raise ValueError("proton_position_A must have shape (3,)")
        distances = np.linalg.norm(oxygen_positions - proton_position, axis=1)
        selected = np.where(distances <= radius_A)[0]
        if max_waters is not None and selected.size > max_waters:
            order = np.argsort(distances[selected])
            selected = selected[order[:max_waters]]
        if selected.size == 0:
            selected = np.array([int(np.argmin(distances))])
        return QMMMRegion(qm_atom_indices=tuple(int(idx) for idx in selected), charge=1)

    def build_config(
        self,
        region: QMMMRegion,
        qm_method: str = "BLYP-D3",
        embedding: str = "electrostatic",
    ) -> QMMMConfig:
        """Create a QM/MM configuration with HydroOnc defaults."""

        if embedding not in {"electrostatic", "mechanical"}:
            raise ValueError("embedding must be electrostatic or mechanical")
        return QMMMConfig(region=region, qm_method=qm_method, embedding=embedding)

    def cp2k_aimd_input(
        self,
        project_name: str,
        coordinates_xyz: str,
        config: Optional[AIMDConfig] = None,
    ) -> str:
        """Generate a CP2K input file for AIMD validation."""

        config = config or AIMDConfig()
        lines = [
            "&GLOBAL",
            f"  PROJECT {project_name}",
            "  RUN_TYPE MD",
            "&END GLOBAL",
            "&MOTION",
            "  &MD",
            f"    ENSEMBLE {config.ensemble}",
            f"    TIMESTEP {config.timestep_fs}",
            f"    STEPS {config.steps}",
            f"    TEMPERATURE {config.temperature_K}",
            "  &END MD",
            "&END MOTION",
            "&FORCE_EVAL",
            "  METHOD Quickstep",
            "  &DFT",
            "    BASIS_SET_FILE_NAME BASIS_MOLOPT",
            "    POTENTIAL_FILE_NAME GTH_POTENTIALS",
            "    &XC",
            "      &XC_FUNCTIONAL",
            f"        &{config.functional}",
            f"        &END {config.functional}",
            "      &END XC_FUNCTIONAL",
            "      &VDW_POTENTIAL",
            "        POTENTIAL_TYPE PAIR_POTENTIAL",
            "        &PAIR_POTENTIAL",
            f"          TYPE {config.dispersion}",
            "        &END PAIR_POTENTIAL",
            "      &END VDW_POTENTIAL",
            "    &END XC",
            "  &END DFT",
            "  &SUBSYS",
            "    &COORD",
        ]
        for line in coordinates_xyz.strip().splitlines():
            lines.append(f"      {line}")
        lines.extend(["    &END COORD", "  &END SUBSYS", "&END FORCE_EVAL"])
        return "\n".join(lines) + "\n"

    @staticmethod
    def validate_recommended_aimd(config: AIMDConfig) -> bool:
        """Return whether settings match the Phase 2 AIMD recommendation."""

        return (
            config.functional.upper() in {"BLYP", "RPBE"}
            and config.dispersion.upper() == "D3"
            and abs(config.timestep_fs - 0.5) < 1.0e-12
            and config.ensemble.upper() == "NVT"
            and abs(config.temperature_K - 300.0) < 1.0e-12
        )

    @staticmethod
    def qm_charge_from_residues(residue_names: Iterable[str]) -> int:
        """Estimate QM-region charge from a simple residue-name list."""

        charge = 1
        for name in residue_names:
            upper = name.upper()
            if upper in {"ASP", "GLU"}:
                charge -= 1
            elif upper in {"LYS", "ARG", "HIS", "HIP"}:
                charge += 1
        return charge
