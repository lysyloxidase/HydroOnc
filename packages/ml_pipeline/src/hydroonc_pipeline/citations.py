"""Literature and reference registry for biological and physical parameters."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Citation:
    """One traceable parameter citation."""

    key: str
    parameter: str
    value: str
    source: str
    module: str


CITATIONS = [
    Citation("codata_2018", "CODATA constants", "2018 SI constants", "CODATA 2018 recommended values", "quantum"),
    Citation("nist_pstar", "proton water stopping power", "70-250 MeV PSTAR-style references", "NIST PSTAR tables", "therapy"),
    Citation("icru_78", "clinical proton RBE", "1.1", "ICRU Report 78", "therapy"),
    Citation("tumor_ph", "tumor extracellular pH", "6.4-7.0", "CEST MRI tumor acidosis literature", "tumor_ph"),
    Citation("h_proton_diffusion", "proton diffusion in water", "9.31e-5 cm2/s", "standard electrochemical transport data", "proton_transport"),
    Citation("h2_oh_rate", "H2 + hydroxyl radical", "4.2e7 M^-1 s^-1", "Ohsawa molecular hydrogen kinetics literature", "h2"),
    Citation("h2_onoo_rate", "H2 + peroxynitrite", "3.6e4 M^-1 s^-1", "molecular hydrogen ROS kinetics literature", "h2"),
    Citation("gsh_oh_rate", "GSH + hydroxyl radical", "1e10 M^-1 s^-1", "radiation chemistry scavenger references", "h2"),
    Citation("akagi_h2", "H2 inhalation oncology", "3 h/day", "Akagi 2019/2020 small clinical series", "h2"),
    Citation("chen_h2", "H2 exhaustion reversal", "2 week inhalation", "Chen 2020 NSCLC immune marker study", "h2"),
    Citation("tp53_hotspots", "p53 hotspot residues", "R175/G245/R248/R249/R273/R282", "TP53 structural oncology literature", "hbond"),
    Citation("kras_g12c", "KRAS G12C His95 H-bond", "sotorasib/adagrasib switch-II pocket", "KRAS G12C co-crystal literature", "hbond"),
    Citation("egfr_t790m", "EGFR T790M gatekeeper", "gefitinib H-bond loss", "EGFR resistance structural literature", "hbond"),
    Citation("fast01", "FLASH FAST-01", "n=10, 8 Gy single fraction", "JAMA Oncology 2023 FAST-01", "therapy"),
]


def citations_by_module(module: str) -> list[Citation]:
    """Return citations for one HydroOnc module."""

    return [citation for citation in CITATIONS if citation.module == module]
