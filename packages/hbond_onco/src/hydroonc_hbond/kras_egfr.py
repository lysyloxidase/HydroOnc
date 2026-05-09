"""H-bond analysis for KRAS and EGFR resistance mutations."""

from __future__ import annotations

from dataclasses import dataclass

from hydroonc_hbond.geometry import HBondInteraction, hbond_geometry_passes


@dataclass(frozen=True)
class KRASHBondResult:
    """KRAS inhibitor H-bond geometry result."""

    mutation: str
    inhibitor: str
    residue: str
    distance_A: float
    angle_deg: float
    occupancy: float
    resistance_flag: bool

    @property
    def matches_cocrystal(self) -> bool:
        return hbond_geometry_passes(self.distance_A, self.angle_deg) and self.occupancy >= 0.50


@dataclass(frozen=True)
class EGFRGatekeeperResult:
    """EGFR gatekeeper H-bond result."""

    mutation: str
    inhibitor: str
    gatekeeper: str
    distance_A: float
    angle_deg: float
    occupancy: float
    abolished: bool
    bypass: str


class KRASMutantAnalyzer:
    """H-bond analysis for KRAS G12C/D/V and EGFR resistance mutations."""

    KRAS_GEOMETRY = {
        ("G12C", "sotorasib"): KRASHBondResult("G12C", "sotorasib", "His95", 2.96, 142.0, 0.62, False),
        ("G12C", "adagrasib"): KRASHBondResult("G12C", "adagrasib", "His95", 2.82, 158.0, 0.78, False),
        ("G12C+Y96D", "sotorasib"): KRASHBondResult("G12C+Y96D", "sotorasib", "His95", 3.72, 103.0, 0.08, True),
        ("G12C+H95Q", "sotorasib"): KRASHBondResult("G12C+H95Q", "sotorasib", "His95", 3.64, 110.0, 0.05, True),
    }

    EGFR_GEOMETRY = {
        ("WT", "gefitinib"): EGFRGatekeeperResult("WT", "gefitinib", "T790", 2.88, 156.0, 0.74, False, "none"),
        ("T790M", "gefitinib"): EGFRGatekeeperResult("T790M", "gefitinib", "M790", 4.18, 82.0, 0.02, True, "none"),
        ("T790M", "osimertinib"): EGFRGatekeeperResult("T790M", "osimertinib", "M790", 3.92, 95.0, 0.03, True, "C797 covalent bond"),
        ("C797S", "osimertinib"): EGFRGatekeeperResult("C797S", "osimertinib", "M790/C797S", 4.05, 91.0, 0.02, True, "lost covalent partner"),
    }

    def kras_hbond_geometry(self, mutation: str = "G12C", inhibitor: str = "sotorasib") -> KRASHBondResult:
        """Return KRAS switch-II H-bond geometry for an inhibitor."""

        key = (mutation.upper(), inhibitor.lower())
        if key not in self.KRAS_GEOMETRY:
            raise ValueError(f"unknown KRAS/inhibitor pair: {mutation}/{inhibitor}")
        return self.KRAS_GEOMETRY[key]

    def egfr_gatekeeper_hbond(self, mutation: str = "T790M", inhibitor: str = "gefitinib") -> EGFRGatekeeperResult:
        """Return EGFR gatekeeper H-bond geometry."""

        key = (mutation.upper(), inhibitor.lower())
        if key not in self.EGFR_GEOMETRY:
            raise ValueError(f"unknown EGFR/inhibitor pair: {mutation}/{inhibitor}")
        return self.EGFR_GEOMETRY[key]

    def resistance_network(self, mutation: str) -> list[str]:
        """Return resistance mutations disrupting KRAS H-bond networks."""

        key = mutation.upper()
        if key == "G12C":
            return ["Y96D/S", "R68S", "H95Q/R"]
        if key in {"G12D", "G12V"}:
            return ["loss of glycine-loop geometry", "altered switch-II water network"]
        if key == "T790M":
            return ["gatekeeper H-bond loss", "steric block", "C797S covalent escape"]
        raise ValueError(f"unknown resistance network: {mutation}")

    def interaction_fixture(self, target: str) -> tuple[HBondInteraction, ...]:
        """Return H-bond interactions for automated count checks."""

        key = target.lower()
        if key == "kras_g12c_sotorasib":
            result = self.kras_hbond_geometry("G12C", "sotorasib")
            return (
                HBondInteraction("KRAS_H95_sotorasib", "H95:NE2", "Sotorasib:O", result.distance_A, result.angle_deg, result.occupancy),
            )
        if key == "egfr_wt_gefitinib":
            result = self.egfr_gatekeeper_hbond("WT", "gefitinib")
            return (
                HBondInteraction("EGFR_T790_gefitinib", "T790:OG1", "Gefitinib:N", result.distance_A, result.angle_deg, result.occupancy),
            )
        if key == "egfr_t790m_gefitinib":
            result = self.egfr_gatekeeper_hbond("T790M", "gefitinib")
            return (
                HBondInteraction("EGFR_M790_gefitinib", "M790:SD", "Gefitinib:N", result.distance_A, result.angle_deg, result.occupancy),
            )
        raise ValueError(f"unknown interaction fixture: {target}")
