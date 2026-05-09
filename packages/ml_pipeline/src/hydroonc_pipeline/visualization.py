"""Visualization page manifest and validation helpers."""

from __future__ import annotations

from dataclasses import dataclass

from hydroonc_pipeline.caveats import CAVEATS


@dataclass(frozen=True)
class WebPageSpec:
    """One HydroOnc visualization page."""

    slug: str
    title: str
    modules: tuple[str, ...]
    checks: tuple[str, ...]
    caveat_indices: tuple[int, ...]

    def as_dict(self) -> dict[str, object]:
        return {
            "slug": self.slug,
            "title": self.title,
            "modules": list(self.modules),
            "checks": list(self.checks),
            "caveats": [CAVEATS[index] for index in self.caveat_indices],
        }


WEB_PAGES = [
    WebPageSpec(
        "quantum",
        "Quantum Explorer",
        ("3D orbital viewer", "emission spectrum", "H2 PEC", "TDSE wavepacket"),
        ("1s spherical", "2p dumbbell", "3d multi-lobed"),
        (),
    ),
    WebPageSpec(
        "proton-transport",
        "Proton Transport",
        ("Grotthuss animation", "H-bond network", "pH gradient volume"),
        ("1 ps^-1 hop rate", "water-wire proton transfer"),
        (),
    ),
    WebPageSpec(
        "tumor-ph",
        "Tumor pH Dashboard",
        ("Warburg sliders", "reaction-diffusion animation", "CEST overlay", "inhibitor playground"),
        ("pH_e 6.4-7.0", "CAIX inhibition raises pH"),
        (),
    ),
    WebPageSpec(
        "proton-therapy",
        "Proton Therapy",
        ("Bragg peak", "SOBP builder", "3D dose", "RBE heatmap", "FLASH comparison"),
        ("150 MeV peak near 15.6 cm", "PSTAR within 2%"),
        (2, 3),
    ),
    WebPageSpec(
        "hbond",
        "Oncoprotein H-Bonds",
        ("Mol* protein viewer", "H-bond difference map", "KRAS drug geometry", "GNN affinity"),
        ("p53 R175H lost H-bonds", "KRAS His95 geometry"),
        (1,),
    ),
    WebPageSpec(
        "h2-therapy",
        "H2 Therapy",
        ("ROS dynamics", "immune ODE", "clinical dashboard"),
        ("H2 + anti-PD-1 synergy", "GSH dominates OH quenching"),
        (0, 4, 5),
    ),
]


def validate_visualization_manifest() -> dict[str, object]:
    """Validate the six-page visualization manifest."""

    slugs = {page.slug for page in WEB_PAGES}
    expected = {"quantum", "proton-transport", "tumor-ph", "proton-therapy", "hbond", "h2-therapy"}
    pages = [page.as_dict() for page in WEB_PAGES]
    return {
        "valid": slugs == expected and all(page["caveats"] is not None for page in pages),
        "page_count": len(WEB_PAGES),
        "slugs": sorted(slugs),
        "pages": pages,
    }
