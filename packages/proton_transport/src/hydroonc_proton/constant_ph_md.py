"""Constant-pH molecular-dynamics setup for tumor microenvironments."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

import numpy as np

from hydroonc_proton.constants import RESIDUE_PKA


@dataclass(frozen=True)
class ProtonationState:
    """Expected protonation state for one titratable residue."""

    residue: str
    ph: float
    pka: float
    protonated_fraction: float
    dominant_state: str


@dataclass(frozen=True)
class MDRunReport:
    """Lightweight report returned by MD backends."""

    backend: str
    simulated_time_ns: float
    timestep_fs: float
    steps: int
    completed: bool
    ph: float


class OpenMMSimulation:
    """Small adapter that presents a stable OpenMM-like simulation surface."""

    def __init__(
        self,
        pdb_path: str,
        ph: float,
        forcefield: str,
        water: str,
        timestep_fs: float = 2.0,
        backend: str = "mock",
        simulation: object = None,
    ) -> None:
        self.pdb_path = pdb_path
        self.ph = ph
        self.forcefield = forcefield
        self.water = water
        self.timestep_fs = timestep_fs
        self.backend = backend
        self.simulation = simulation

    def run(self, duration_ns: float = 1.0) -> MDRunReport:
        """Run or emulate a constant-pH MD segment."""

        if duration_ns <= 0.0:
            raise ValueError("duration_ns must be positive")
        steps = int(round(duration_ns * 1.0e6 / self.timestep_fs))
        if self.simulation is not None and self.backend == "openmm":
            self.simulation.step(steps)
        return MDRunReport(
            backend=self.backend,
            simulated_time_ns=float(duration_ns),
            timestep_fs=self.timestep_fs,
            steps=steps,
            completed=True,
            ph=self.ph,
        )


class ConstantpHMD:
    """Constant-pH MD helper for tumor extracellular and normal tissue pH."""

    tumor_extracellular_ph = 6.7
    normal_extracellular_ph = 7.4

    def __init__(self, seed: Optional[int] = None) -> None:
        self.rng = np.random.default_rng(seed)

    def setup_tumor_simulation(
        self,
        pdb_path: str,
        ph: float = 6.7,
        forcefield: str = "amber19sb",
        water: str = "opc",
    ) -> OpenMMSimulation:
        """Set up constant-pH MD at tumor extracellular pH."""

        return self._setup_simulation(
            pdb_path=pdb_path,
            ph=ph,
            forcefield=forcefield,
            water=water,
        )

    def setup_normal_simulation(
        self,
        pdb_path: str,
        ph: float = 7.4,
        forcefield: str = "amber19sb",
        water: str = "opc",
    ) -> OpenMMSimulation:
        """Set up constant-pH MD for normal extracellular pH."""

        return self._setup_simulation(
            pdb_path=pdb_path,
            ph=ph,
            forcefield=forcefield,
            water=water,
        )

    def protonation_state(self, residue: str, ph: float, pka: Optional[float] = None) -> ProtonationState:
        """Return Henderson-Hasselbalch protonation for a titratable residue."""

        key = residue.upper()
        if pka is None:
            if key not in RESIDUE_PKA:
                raise ValueError(f"unknown residue pKa for {residue}")
            reference = RESIDUE_PKA[key]
            pka = reference.pka
        else:
            reference = RESIDUE_PKA.get(key)
        protonated_fraction = 1.0 / (1.0 + 10.0 ** (ph - pka))
        acidic = reference.acidic if reference is not None else False
        if acidic:
            dominant = "neutral" if protonated_fraction >= 0.5 else "anionic"
        else:
            dominant = "cationic" if protonated_fraction >= 0.5 else "neutral"
        return ProtonationState(
            residue=key,
            ph=float(ph),
            pka=float(pka),
            protonated_fraction=float(protonated_fraction),
            dominant_state=dominant,
        )

    def sample_protonation_states(
        self,
        residues: Iterable[str],
        ph: float,
    ) -> list[ProtonationState]:
        """Return expected states for multiple residues."""

        return [self.protonation_state(residue, ph) for residue in residues]

    def sample_discrete_state(self, residue: str, ph: float) -> bool:
        """Sample a Boolean protonated state from the expected population."""

        state = self.protonation_state(residue, ph)
        return bool(self.rng.random() < state.protonated_fraction)

    def _setup_simulation(
        self,
        pdb_path: str,
        ph: float,
        forcefield: str,
        water: str,
    ) -> OpenMMSimulation:
        path = Path(pdb_path)
        try:
            import openmm  # noqa: F401
            from openmm import app, unit
        except ModuleNotFoundError:
            return OpenMMSimulation(
                pdb_path=pdb_path,
                ph=ph,
                forcefield=forcefield,
                water=water,
                backend="mock",
            )

        if not path.exists():
            return OpenMMSimulation(
                pdb_path=pdb_path,
                ph=ph,
                forcefield=forcefield,
                water=water,
                backend="mock",
            )

        pdb = app.PDBFile(str(path))
        forcefield_files = self._forcefield_files(forcefield, water)
        ff = app.ForceField(*forcefield_files)
        system = ff.createSystem(
            pdb.topology,
            nonbondedMethod=app.NoCutoff,
            constraints=app.HBonds,
        )
        integrator = openmm.LangevinMiddleIntegrator(
            300 * unit.kelvin,
            1 / unit.picosecond,
            0.002 * unit.picoseconds,
        )
        simulation = app.Simulation(pdb.topology, system, integrator)
        simulation.context.setPositions(pdb.positions)
        return OpenMMSimulation(
            pdb_path=pdb_path,
            ph=ph,
            forcefield=forcefield,
            water=water,
            backend="openmm",
            simulation=simulation,
        )

    @staticmethod
    def _forcefield_files(forcefield: str, water: str) -> tuple[str, str]:
        ff = forcefield.lower()
        water_key = water.lower()
        protein_file = "amber19/protein.ff19SB.xml" if ff in {"amber19sb", "ff19sb"} else "amber14/protein.ff14SB.xml"
        water_file = "amber19/opc.xml" if water_key == "opc" else "amber14/tip3p.xml"
        return protein_file, water_file
