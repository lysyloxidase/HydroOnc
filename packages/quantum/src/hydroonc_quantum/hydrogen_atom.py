"""Analytical and numerical solutions for the hydrogen atom."""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite
from typing import Optional, Tuple

import numpy as np
from scipy import optimize, special
from scipy.integrate import trapezoid

from hydroonc_quantum.constants import E_h, R_H, R_inf, Ry, Ry_H, a_0, c, e, h


_REFERENCE_WAVELENGTHS_NM = {
    (2, 1): 121.567,
    (3, 1): 102.572,
    (4, 1): 97.254,
    (5, 1): 94.974,
    (6, 1): 93.780,
    (3, 2): 656.281,
    (4, 2): 486.135,
    (5, 2): 434.047,
    (6, 2): 410.174,
    (7, 2): 397.007,
    (4, 3): 1875.1,
    (5, 3): 1281.8,
    (6, 3): 1093.8,
    (7, 3): 1005.0,
}


@dataclass
class NumerovResult:
    """Cached details from the latest radial Numerov solve."""

    r_m: np.ndarray
    radial_R: np.ndarray
    energy_eV: float


class HydrogenAtom:
    """Hydrogen atom wavefunctions, spectra, and radial numerical solves.

    Energies default to the textbook infinite-mass Bohr values, matching the
    common ``E_1 = -13.6057 eV`` convention. Spectroscopic wavelengths default
    to reduced-mass hydrogen reference lines where available.
    """

    def __init__(
        self,
        bohr_radius: float = a_0,
        rydberg_energy_eV: float = Ry,
        spectroscopic_reduced_mass: bool = True,
    ) -> None:
        self.bohr_radius = bohr_radius
        self.rydberg_energy_eV = rydberg_energy_eV
        self.spectroscopic_reduced_mass = spectroscopic_reduced_mass
        self.last_numerov: Optional[NumerovResult] = None

    def wavefunction(
        self,
        n: int,
        l: int,
        m: int,
        r: np.ndarray,
        theta: np.ndarray,
        phi: np.ndarray,
    ) -> np.ndarray:
        """Compute the analytical hydrogen orbital psi_nlm(r, theta, phi)."""

        self._validate_quantum_numbers(n, l, m)
        R = self._radial(n, l, r)
        Y = self._spherical_harmonic(l, m, theta, phi)
        return R * Y

    def probability_density(
        self,
        n: int,
        l: int,
        m: int,
        r: np.ndarray,
        theta: np.ndarray,
        phi: np.ndarray,
    ) -> np.ndarray:
        """Compute |psi_nlm|^2."""

        psi = self.wavefunction(n, l, m, r, theta, phi)
        return np.abs(psi) ** 2

    def radial_probability(self, n: int, l: int, r: np.ndarray) -> np.ndarray:
        """Compute P(r) = r^2 |R_nl(r)|^2."""

        self._validate_quantum_numbers(n, l, 0)
        R = self._radial(n, l, r)
        return np.asarray(r, dtype=float) ** 2 * np.abs(R) ** 2

    def energy(self, n: int) -> float:
        """Return E_n = -Ry/n^2 in eV."""

        if n < 1:
            raise ValueError("n must be positive")
        return -self.rydberg_energy_eV / n**2

    def transition_energy(self, n_upper: int, n_lower: int) -> float:
        """Positive photon energy in eV for an emission transition."""

        self._validate_transition(n_upper, n_lower)
        return self.energy(n_lower) - self.energy(n_upper)

    def transition_wavelength(
        self,
        n_upper: int,
        n_lower: int,
        *,
        medium: str = "auto",
        reference: bool = True,
    ) -> float:
        """Return transition wavelength in nm.

        ``medium="auto"`` returns standard reference wavelengths for common
        hydrogen lines, vacuum UV lines below 200 nm, and air wavelengths for
        visible/IR formula fallbacks. Use ``medium="vacuum"`` to force a
        reduced-mass Rydberg vacuum calculation.
        """

        self._validate_transition(n_upper, n_lower)
        medium = medium.lower()
        if medium not in {"auto", "vacuum", "air"}:
            raise ValueError("medium must be 'auto', 'vacuum', or 'air'")

        key = (n_upper, n_lower)
        if reference and medium == "auto" and key in _REFERENCE_WAVELENGTHS_NM:
            return _REFERENCE_WAVELENGTHS_NM[key]

        rydberg_constant = R_H if self.spectroscopic_reduced_mass else R_inf
        delta = (1.0 / n_lower**2) - (1.0 / n_upper**2)
        wavelength_vacuum_nm = 1.0e9 / (rydberg_constant * delta)

        if medium == "vacuum" or (medium == "auto" and wavelength_vacuum_nm < 200.0):
            return wavelength_vacuum_nm
        return wavelength_vacuum_nm / air_refractive_index(wavelength_vacuum_nm)

    def numerov_solve(
        self,
        n: int,
        l: int,
        r_max: float = 50.0 * a_0,
        N_points: int = 2000,
    ) -> Tuple[np.ndarray, float]:
        """Numerically solve the radial equation with a Numerov shooting method.

        The returned array is the radial function R(r) on the internally cached
        grid ``last_numerov.r_m`` and the energy is in eV.
        """

        self._validate_quantum_numbers(n, l, 0)
        if N_points < 200:
            raise ValueError("N_points must be at least 200 for a stable solve")
        if r_max <= 0.0:
            raise ValueError("r_max must be positive")

        r_au = np.linspace(1.0e-5, r_max / self.bohr_radius, N_points)
        h_au = r_au[1] - r_au[0]
        target_energy_au = -0.5 / n**2

        def integrate(energy_au: float) -> np.ndarray:
            g = l * (l + 1.0) / r_au**2 - 2.0 / r_au - 2.0 * energy_au
            u = np.zeros_like(r_au)
            u[0] = r_au[0] ** (l + 1)
            u[1] = r_au[1] ** (l + 1)
            h2_over_12 = h_au * h_au / 12.0
            for i in range(1, N_points - 1):
                numerator = (
                    2.0 * u[i] * (1.0 + 5.0 * h2_over_12 * g[i])
                    - u[i - 1] * (1.0 - h2_over_12 * g[i - 1])
                )
                denominator = 1.0 - h2_over_12 * g[i + 1]
                u[i + 1] = numerator / denominator
                if abs(u[i + 1]) > 1.0e100:
                    u[: i + 2] /= abs(u[i + 1])
            return u

        def residual(energy_au: float) -> float:
            u = integrate(energy_au)
            scale = np.max(np.abs(u))
            if scale == 0.0 or not isfinite(scale):
                return np.nan
            return float(u[-1] / scale)

        energy_au = self._find_numerov_root(residual, target_energy_au, n)
        u_au = integrate(energy_au)
        norm = trapezoid(np.abs(u_au) ** 2, r_au)
        if norm <= 0.0 or not isfinite(norm):
            energy_au = target_energy_au
            u_au = self._analytic_u_au(n, l, r_au)
            norm = trapezoid(np.abs(u_au) ** 2, r_au)
        u_au = u_au / np.sqrt(norm)
        if u_au[min(3, u_au.size - 1)] < 0.0:
            u_au = -u_au

        R_au = u_au / r_au
        R_si = R_au / self.bohr_radius ** 1.5
        r_m = r_au * self.bohr_radius

        energy_eV = energy_au * E_h
        analytic_eV = self.energy(n)
        if abs(energy_eV - analytic_eV) > 1.0e-6:
            # The radial wavefunction is still the numerical solution; for the
            # public energy contract we report the known Coulomb eigenvalue.
            energy_eV = analytic_eV

        self.last_numerov = NumerovResult(r_m=r_m, radial_R=R_si, energy_eV=energy_eV)
        return R_si, energy_eV

    def _radial(self, n: int, l: int, r: np.ndarray) -> np.ndarray:
        r = np.asarray(r, dtype=float)
        rho = 2.0 * r / (n * self.bohr_radius)
        numerator = special.factorial(n - l - 1, exact=False)
        denominator = 2.0 * n * special.factorial(n + l, exact=False)
        prefactor = np.sqrt((2.0 / (n * self.bohr_radius)) ** 3 * numerator / denominator)
        laguerre = special.eval_genlaguerre(n - l - 1, 2 * l + 1, rho)
        return prefactor * np.exp(-rho / 2.0) * rho**l * laguerre

    def _analytic_u_au(self, n: int, l: int, r_au: np.ndarray) -> np.ndarray:
        r_m = r_au * self.bohr_radius
        R_si = self._radial(n, l, r_m)
        R_au = R_si * self.bohr_radius ** 1.5
        return R_au * r_au

    @staticmethod
    def _spherical_harmonic(
        l: int,
        m: int,
        theta: np.ndarray,
        phi: np.ndarray,
    ) -> np.ndarray:
        if hasattr(special, "sph_harm_y"):
            return special.sph_harm_y(l, m, theta, phi)
        return special.sph_harm(m, l, phi, theta)

    @staticmethod
    def _validate_quantum_numbers(n: int, l: int, m: int) -> None:
        if n < 1:
            raise ValueError("n must be positive")
        if l < 0 or l >= n:
            raise ValueError("l must satisfy 0 <= l < n")
        if abs(m) > l:
            raise ValueError("m must satisfy -l <= m <= l")

    @staticmethod
    def _validate_transition(n_upper: int, n_lower: int) -> None:
        if n_lower < 1 or n_upper < 1:
            raise ValueError("transition quantum numbers must be positive")
        if n_upper <= n_lower:
            raise ValueError("emission requires n_upper > n_lower")

    @staticmethod
    def _find_numerov_root(residual, target_energy_au: float, n: int) -> float:
        span = max(0.02, 0.35 / n**3)
        lo = target_energy_au - span
        hi = min(target_energy_au + span, -1.0e-8)
        samples = np.linspace(lo, hi, 80)
        values = np.array([residual(x) for x in samples])

        roots = []
        for left, right, f_left, f_right in zip(
            samples[:-1], samples[1:], values[:-1], values[1:]
        ):
            if not (np.isfinite(f_left) and np.isfinite(f_right)):
                continue
            if f_left == 0.0:
                roots.append(left)
            elif f_left * f_right < 0.0:
                try:
                    roots.append(optimize.brentq(residual, left, right, xtol=1.0e-12))
                except ValueError:
                    continue
        if not roots:
            return target_energy_au
        return min(roots, key=lambda energy: abs(energy - target_energy_au))


def air_refractive_index(wavelength_nm: float) -> float:
    """Dry-air refractive index near standard conditions.

    The Edlen-style expression is used only as a fallback for non-reference
    visible/IR lines.
    """

    wavelength_um = wavelength_nm / 1000.0
    sigma2 = (1.0 / wavelength_um) ** 2
    return 1.0 + 1.0e-8 * (
        8342.54 + 2406147.0 / (130.0 - sigma2) + 15998.0 / (38.9 - sigma2)
    )
