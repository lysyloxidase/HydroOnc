"""Hydrogen-bond geometry analysis for biomolecular structures."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np


@dataclass(frozen=True)
class AtomSite:
    """Minimal atom representation used by the fallback analyzer."""

    name: str
    element: str
    position_A: np.ndarray
    residue: str = ""


@dataclass(frozen=True)
class HBond:
    """One geometric hydrogen bond."""

    donor: AtomSite
    hydrogen: AtomSite
    acceptor: AtomSite
    distance_DA_A: float
    angle_DHA_deg: float
    strength: str


class HBondAnalyzer:
    """Hydrogen-bond analysis in MD trajectories and static structures."""

    def __init__(
        self,
        distance_range_A: tuple[float, float] = (2.5, 3.5),
        angle_cutoff_deg: float = 120.0,
        strong_angle_cutoff_deg: float = 150.0,
    ) -> None:
        self.distance_range_A = distance_range_A
        self.angle_cutoff_deg = angle_cutoff_deg
        self.strong_angle_cutoff_deg = strong_angle_cutoff_deg

    def identify(
        self,
        donors: Iterable[tuple[AtomSite, AtomSite]],
        acceptors: Iterable[AtomSite],
    ) -> list[HBond]:
        """Identify hydrogen bonds from donor-hydrogen pairs and acceptors."""

        bonds = []
        for donor, hydrogen in donors:
            for acceptor in acceptors:
                if acceptor is donor:
                    continue
                distance = self.distance(donor.position_A, acceptor.position_A)
                angle = self.angle(donor.position_A, hydrogen.position_A, acceptor.position_A)
                if self._passes(distance, angle):
                    strength = "strong" if angle >= self.strong_angle_cutoff_deg else "weak"
                    bonds.append(HBond(donor, hydrogen, acceptor, distance, angle, strength))
        return bonds

    def analyze_base_pair(self, base_pair: str) -> list[HBond]:
        """Return geometric H-bonds for reference DNA base-pair templates."""

        donors, acceptors = reference_base_pair_sites(base_pair)
        return self.identify(donors, acceptors)

    def analyze_mdanalysis_universe(
        self,
        universe: object,
        selection: str = "all",
    ) -> list[HBond]:
        """Analyze an MDAnalysis universe when MDAnalysis is available."""

        try:
            import MDAnalysis as mda  # noqa: F401
        except ModuleNotFoundError as exc:
            raise RuntimeError("Install hydroonc-proton-transport[md] for MDAnalysis support") from exc

        atoms = universe.select_atoms(selection)
        donors = []
        acceptors = []
        for atom in atoms:
            element = (getattr(atom, "element", "") or atom.name[0]).upper()
            if element in {"O", "N"}:
                site = AtomSite(atom.name, element, np.asarray(atom.position, dtype=float), atom.resname)
                acceptors.append(site)
                for bonded in getattr(atom, "bonded_atoms", []):
                    bonded_element = (getattr(bonded, "element", "") or bonded.name[0]).upper()
                    if bonded_element == "H":
                        donors.append(
                            (
                                site,
                                AtomSite(
                                    bonded.name,
                                    "H",
                                    np.asarray(bonded.position, dtype=float),
                                    bonded.resname,
                                ),
                            )
                        )
        return self.identify(donors, acceptors)

    @staticmethod
    def distance(a: np.ndarray, b: np.ndarray) -> float:
        return float(np.linalg.norm(np.asarray(a, dtype=float) - np.asarray(b, dtype=float)))

    @staticmethod
    def angle(donor: np.ndarray, hydrogen: np.ndarray, acceptor: np.ndarray) -> float:
        donor = np.asarray(donor, dtype=float)
        hydrogen = np.asarray(hydrogen, dtype=float)
        acceptor = np.asarray(acceptor, dtype=float)
        vector_hd = donor - hydrogen
        vector_ha = acceptor - hydrogen
        denom = np.linalg.norm(vector_hd) * np.linalg.norm(vector_ha)
        if denom == 0.0:
            return 0.0
        cosine = np.dot(vector_hd, vector_ha) / denom
        return float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))

    def _passes(self, distance_A: float, angle_deg: float) -> bool:
        return (
            self.distance_range_A[0] <= distance_A <= self.distance_range_A[1]
            and angle_deg >= self.angle_cutoff_deg
        )


def reference_base_pair_sites(base_pair: str) -> tuple[list[tuple[AtomSite, AtomSite]], list[AtomSite]]:
    """Return reference donor/hydrogen pairs and acceptors for A-T or G-C."""

    key = base_pair.upper().replace("-", "")
    if key in {"AT", "TA"}:
        donors = [
            _donor("A:N6", "N", (0.0, 0.0, 0.0), "A:H61", (1.0, 0.0, 0.0)),
            _donor("T:N3", "N", (0.0, 3.0, 0.0), "T:H3", (1.0, 3.0, 0.0)),
        ]
        acceptors = [
            AtomSite("T:O4", "O", np.array((2.90, 0.0, 0.0)), "THY"),
            AtomSite("A:N1", "N", np.array((2.84, 3.0, 0.0)), "ADE"),
        ]
        return donors, acceptors
    if key in {"GC", "CG"}:
        donors = [
            _donor("G:N1", "N", (0.0, 0.0, 0.0), "G:H1", (1.0, 0.0, 0.0)),
            _donor("G:N2", "N", (0.0, 3.0, 0.0), "G:H21", (1.0, 3.0, 0.0)),
            _donor("C:N4", "N", (0.0, 6.0, 0.0), "C:H41", (1.0, 6.0, 0.0)),
        ]
        acceptors = [
            AtomSite("C:N3", "N", np.array((2.86, 0.0, 0.0)), "CYT"),
            AtomSite("C:O2", "O", np.array((2.95, 3.0, 0.0)), "CYT"),
            AtomSite("G:O6", "O", np.array((2.91, 6.0, 0.0)), "GUA"),
        ]
        return donors, acceptors
    raise ValueError("base_pair must be A-T or G-C")


def _donor(
    donor_name: str,
    donor_element: str,
    donor_position: tuple[float, float, float],
    hydrogen_name: str,
    hydrogen_position: tuple[float, float, float],
) -> tuple[AtomSite, AtomSite]:
    donor = AtomSite(donor_name, donor_element, np.array(donor_position, dtype=float), donor_name[0])
    hydrogen = AtomSite(hydrogen_name, "H", np.array(hydrogen_position, dtype=float), hydrogen_name[0])
    return donor, hydrogen
