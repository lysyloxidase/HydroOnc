"""p53 hotspot mutation H-bond disruption analysis."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Union

import numpy as np

from hydroonc_hbond.fixtures import P53_MUTANT_HBONDS, P53_WT_2XWR_HBONDS
from hydroonc_hbond.geometry import HBondInteraction, automated_hbond_count, hbond_occupancy_map


@dataclass(frozen=True)
class MutationHBondComparison:
    """Wild-type versus mutant H-bond comparison."""

    mutation: str
    mutation_type: str
    effect: str
    wild_type_count: int
    mutant_count: int
    lost_hbonds: tuple[str, ...]
    gained_hbonds: tuple[str, ...]
    dna_contact_occupancy: float
    tumor_ph: float
    drug_target: str

    @property
    def hbond_loss_count(self) -> int:
        return len(self.lost_hbonds)


class P53MutantAnalyzer:
    """Analyze hydrogen bond disruption in p53 cancer mutations."""

    MUTATIONS = {
        "R175H": {
            "type": "structural",
            "effect": "zinc_loss",
            "pdb_wt": "2XWR",
            "drug": "ZMC1 (zinc metallochaperone)",
        },
        "R273H": {
            "type": "contact",
            "effect": "dna_contact_loss",
            "drug": "PC14586 (Rezatapopt, Phase 2)",
        },
        "R248W": {"type": "contact+structural", "effect": "minor_groove_loss", "drug": "DNA-contact rescue"},
        "G245S": {"type": "structural", "effect": "loop_displacement", "drug": "L3-loop stabilizer"},
        "R249S": {"type": "structural", "effect": "L3_destabilization", "drug": "core-domain stabilizer"},
        "R282W": {"type": "structural", "effect": "H2_helix_disruption", "drug": "helix rescue binder"},
    }

    DRUG_HISTORY = {
        "APR-246 (eprenetapopt)": {
            "mechanism": "Michael acceptor MQ covalently modifies Cys124/Cys277",
            "phase2": "71% ORR, 44% CR with azacitidine in TP53-mutant MDS",
            "phase3": "FAILED: 33.3% vs 22.4% CR, p=0.13, Aprea 2021",
            "lesson": "Mechanism-driven optimism does not guarantee Phase 3 success",
        },
        "Sotorasib (KRAS G12C)": {
            "mechanism": "Covalent binding to Cys12 switch-II pocket in GDP-bound KRAS",
            "trial": "CodeBreaK 200: mPFS 5.6 vs 4.5 months, HR 0.66",
            "caveat": "PFS benefit without overall-survival benefit",
        },
    }

    def __init__(self, occupancy_cutoff: float = 0.30) -> None:
        self.occupancy_cutoff = occupancy_cutoff

    def md_occupancy_map(self, mutation: str = "WT", ph: float = 7.4) -> dict[str, float]:
        """Return a deterministic MD-style H-bond occupancy map."""

        interactions = self._interactions_for(mutation)
        ph_penalty = 0.04 if ph <= 6.8 else 0.0
        occupancies = {}
        for interaction in interactions:
            penalty = ph_penalty if interaction.category in {"sidechain", "dna_contact"} else ph_penalty / 2.0
            occupancies[interaction.identifier] = float(np.clip(interaction.occupancy - penalty, 0.0, 1.0))
        return occupancies

    def analyze_mutation(self, mutation: str, ph: float = 7.4) -> MutationHBondComparison:
        """Compare wild-type 2XWR H-bond occupancy against a p53 mutant."""

        key = mutation.upper()
        if key not in self.MUTATIONS:
            raise ValueError(f"unknown p53 hotspot mutation: {mutation}")
        wild_type = tuple(P53_WT_2XWR_HBONDS)
        mutant = self._interactions_for(key)
        wt_present = self._present_identifiers(wild_type)
        mutant_present = self._present_identifiers(mutant)
        lost = tuple(identifier for identifier in wt_present if identifier not in mutant_present)
        gained = tuple(identifier for identifier in mutant_present if identifier not in wt_present)
        meta = self.MUTATIONS[key]
        return MutationHBondComparison(
            mutation=key,
            mutation_type=meta["type"],
            effect=meta["effect"],
            wild_type_count=automated_hbond_count(wild_type, self.occupancy_cutoff),
            mutant_count=automated_hbond_count(mutant, self.occupancy_cutoff),
            lost_hbonds=lost,
            gained_hbonds=gained,
            dna_contact_occupancy=self.dna_contact_occupancy(key),
            tumor_ph=ph,
            drug_target=meta.get("drug", "not established"),
        )

    def dna_contact_occupancy(self, mutation: str, contact_id: str = "R273_DNA_phosphate") -> float:
        """Return occupancy for a named p53-DNA contact."""

        occupancy = hbond_occupancy_map(self._interactions_for(mutation))
        return float(occupancy.get(contact_id, 0.0))

    def identify_druggable_cavities(self, mutation: str) -> list[dict[str, Union[float, str]]]:
        """Identify curated cavities near disrupted H-bond networks."""

        comparison = self.analyze_mutation(mutation)
        if comparison.mutation == "R175H":
            return [
                {
                    "site": "L2/L3 zinc pocket",
                    "priority": 0.92,
                    "strategy": "zinc metallochaperone rescue",
                    "example": "ZMC1",
                }
            ]
        if comparison.mutation in {"R273H", "R248W"}:
            return [
                {
                    "site": "DNA-contact surface",
                    "priority": 0.78,
                    "strategy": "restore phosphate/minor-groove H-bond network",
                    "example": comparison.drug_target,
                }
            ]
        return [
            {
                "site": "core-domain stabilizing pocket",
                "priority": 0.66,
                "strategy": "stabilize disrupted loop or helix H-bond network",
                "example": comparison.drug_target,
            }
        ]

    def hbond_count(self, mutation: str = "WT") -> int:
        """Count occupied H-bonds for wild type or mutant p53."""

        return automated_hbond_count(self._interactions_for(mutation), self.occupancy_cutoff)

    def _present_identifiers(self, interactions: tuple[HBondInteraction, ...]) -> tuple[str, ...]:
        return tuple(
            interaction.identifier
            for interaction in interactions
            if interaction.occupancy >= self.occupancy_cutoff and interaction.is_geometric_hbond
        )

    @staticmethod
    def _interactions_for(mutation: str) -> tuple[HBondInteraction, ...]:
        key = mutation.upper()
        if key in {"WT", "WILD_TYPE", "2XWR"}:
            return tuple(P53_WT_2XWR_HBONDS)
        if key not in P53_MUTANT_HBONDS:
            raise ValueError(f"unknown p53 mutation fixture: {mutation}")
        return tuple(P53_MUTANT_HBONDS[key])
