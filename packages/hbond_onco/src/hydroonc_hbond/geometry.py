"""Shared H-bond geometry and graph primitives."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Sequence

import numpy as np


@dataclass(frozen=True)
class AtomFeature:
    """Minimal atom feature for protein-ligand interaction graphs."""

    atom_id: str
    atom_type: str
    charge: float
    is_donor: bool
    is_acceptor: bool
    hybridization: str
    coord_A: np.ndarray
    residue: str = ""

    def feature_vector(self) -> np.ndarray:
        atom_codes = {"H": 1.0, "C": 6.0, "N": 7.0, "O": 8.0, "S": 16.0, "P": 15.0, "ZN": 30.0}
        hybrid_codes = {"sp": 1.0, "sp2": 2.0, "sp3": 3.0, "metal": 4.0, "unknown": 0.0}
        return np.array(
            [
                atom_codes.get(self.atom_type.upper(), 0.0) / 30.0,
                self.charge,
                float(self.is_donor),
                float(self.is_acceptor),
                hybrid_codes.get(self.hybridization.lower(), 0.0) / 4.0,
            ],
            dtype=float,
        )


@dataclass(frozen=True)
class HBondInteraction:
    """Hydrogen-bond or near-contact interaction."""

    identifier: str
    donor: str
    acceptor: str
    distance_A: float
    angle_deg: float
    occupancy: float
    category: str = "sidechain"
    buried: bool = False
    water_mediated: bool = False

    @property
    def is_geometric_hbond(self) -> bool:
        return hbond_geometry_passes(self.distance_A, self.angle_deg)

    @property
    def is_occupied(self) -> bool:
        return self.occupancy >= 0.30 and self.is_geometric_hbond


@dataclass(frozen=True)
class InteractionGraph:
    """Protein-ligand graph with H-bond-aware edge features."""

    nodes: tuple[AtomFeature, ...]
    edges: tuple[tuple[int, int], ...]
    edge_features: np.ndarray
    node_features: np.ndarray

    @property
    def n_hbond_edges(self) -> int:
        if self.edge_features.size == 0:
            return 0
        return int(np.sum(self.edge_features[:, 1] > 0.5))


def hbond_geometry_passes(distance_A: float, angle_deg: float) -> bool:
    """Return whether standard H-bond geometry criteria are met."""

    return 2.5 <= distance_A <= 3.5 and angle_deg >= 120.0


def automated_hbond_count(interactions: Iterable[HBondInteraction], occupancy_cutoff: float = 0.30) -> int:
    """Count occupied H-bonds from an interaction list."""

    count = 0
    for interaction in interactions:
        if interaction.occupancy >= occupancy_cutoff and interaction.is_geometric_hbond:
            count += 1
    return count


def hbond_occupancy_map(interactions: Iterable[HBondInteraction]) -> dict[str, float]:
    """Return an identifier to occupancy map."""

    return {interaction.identifier: interaction.occupancy for interaction in interactions}


def build_interaction_graph(
    protein_atoms: Sequence[AtomFeature],
    ligand_atoms: Sequence[AtomFeature],
    distance_cutoff_A: float = 4.5,
) -> InteractionGraph:
    """Build a H-bond-aware protein-ligand interaction graph."""

    nodes = tuple(protein_atoms) + tuple(ligand_atoms)
    node_features = np.vstack([node.feature_vector() for node in nodes]) if nodes else np.empty((0, 5))
    edges = []
    edge_features = []
    for i, protein_atom in enumerate(protein_atoms):
        for j, ligand_atom in enumerate(ligand_atoms, start=len(protein_atoms)):
            distance = float(np.linalg.norm(protein_atom.coord_A - ligand_atom.coord_A))
            if distance > distance_cutoff_A:
                continue
            donor_acceptor = (
                protein_atom.is_donor
                and ligand_atom.is_acceptor
                or ligand_atom.is_donor
                and protein_atom.is_acceptor
            )
            angle_proxy = 155.0 if donor_acceptor and 2.5 <= distance <= 3.5 else 90.0
            is_hbond = donor_acceptor and hbond_geometry_passes(distance, angle_proxy)
            hydrophobic = float(protein_atom.atom_type.upper() == "C" and ligand_atom.atom_type.upper() == "C")
            covalent = float(distance < 1.9)
            edges.append((i, j))
            edge_features.append([distance, float(is_hbond), 1.0 if covalent else 0.0, covalent, hydrophobic])
    return InteractionGraph(
        nodes=nodes,
        edges=tuple(edges),
        edge_features=np.asarray(edge_features, dtype=float).reshape((-1, 5)) if edge_features else np.empty((0, 5)),
        node_features=node_features,
    )
