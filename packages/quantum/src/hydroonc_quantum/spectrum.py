"""Hydrogen emission spectrum utilities."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Optional

from hydroonc_quantum.hydrogen_atom import HydrogenAtom


SPECTRAL_SERIES = {
    "Lyman": {
        "n_lower": 1,
        "range": "UV",
        "lines_nm": [121.567, 102.572, 97.254, 94.974, 93.780],
    },
    "Balmer": {
        "n_lower": 2,
        "range": "visible/UV",
        "lines_nm": [656.281, 486.135, 434.047, 410.174, 397.007],
    },
    "Paschen": {
        "n_lower": 3,
        "range": "infrared",
        "lines_nm": [1875.1, 1281.8, 1093.8, 1005.0],
    },
    "Brackett": {"n_lower": 4, "range": "infrared"},
    "Pfund": {"n_lower": 5, "range": "far infrared"},
}


@dataclass(frozen=True)
class EmissionLine:
    """One hydrogen emission line."""

    series: str
    n_upper: int
    n_lower: int
    wavelength_nm: float
    energy_eV: float
    range: str


class HydrogenSpectrum:
    """Generate hydrogen spectral series and validate dipole selection rules."""

    def __init__(self, atom: Optional[HydrogenAtom] = None) -> None:
        self.atom = atom or HydrogenAtom()

    def line(self, n_upper: int, n_lower: int, *, series: Optional[str] = None) -> EmissionLine:
        """Return one emission line."""

        name = series or self.series_name(n_lower)
        meta = SPECTRAL_SERIES.get(name, {"range": "unknown"})
        return EmissionLine(
            series=name,
            n_upper=n_upper,
            n_lower=n_lower,
            wavelength_nm=self.atom.transition_wavelength(n_upper, n_lower),
            energy_eV=self.atom.transition_energy(n_upper, n_lower),
            range=meta.get("range", "unknown"),
        )

    def series(self, name: str, n_max: Optional[int] = None) -> list[EmissionLine]:
        """Return emission lines for a named spectral series."""

        if name not in SPECTRAL_SERIES:
            raise ValueError(f"unknown spectral series: {name}")
        n_lower = int(SPECTRAL_SERIES[name]["n_lower"])
        if n_max is None:
            known = SPECTRAL_SERIES[name].get("lines_nm", [])
            n_max = n_lower + max(len(known), 5)
        return [self.line(n_upper, n_lower, series=name) for n_upper in range(n_lower + 1, n_max + 1)]

    @staticmethod
    def series_name(n_lower: int) -> str:
        """Return the conventional series name for a lower level."""

        for name, meta in SPECTRAL_SERIES.items():
            if meta["n_lower"] == n_lower:
                return name
        return f"n={n_lower}"

    @staticmethod
    def allowed_dipole_transition(
        l_initial: int,
        m_initial: int,
        l_final: int,
        m_final: int,
        spin_initial: float = 0.5,
        spin_final: float = 0.5,
    ) -> bool:
        """Check electric-dipole selection rules: dl = +/-1, dm = 0,+/-1, ds = 0."""

        return (
            abs(l_final - l_initial) == 1
            and abs(m_final - m_initial) <= 1
            and spin_initial == spin_final
        )

    def wavelengths(self, name: str, n_upper_values: Iterable[int]) -> list[float]:
        """Return wavelengths in nm for selected upper levels of one series."""

        if name not in SPECTRAL_SERIES:
            raise ValueError(f"unknown spectral series: {name}")
        n_lower = int(SPECTRAL_SERIES[name]["n_lower"])
        return [self.atom.transition_wavelength(n_upper, n_lower) for n_upper in n_upper_values]
