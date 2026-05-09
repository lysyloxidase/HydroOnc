"""Bragg peak and SOBP simulation for clinical proton beams."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from hydroonc_therapy.constants import (
    BETHE_K_MEV_CM2_G,
    CLINICAL_RANGE_TABLE_WATER_CM,
    ELECTRON_MASS_MEV,
    PROTON_MASS_MEV,
    PSTAR_WATER_MASS_STOPPING_POWER,
    WATER_MEAN_EXCITATION_MEV,
    WATER_Z_OVER_A,
)


@dataclass(frozen=True)
class BraggCurve:
    """Pristine proton depth-dose curve."""

    depth_cm: np.ndarray
    dose: np.ndarray
    energy_MeV: float
    range_cm: float
    let_keV_um: np.ndarray

    @property
    def peak_depth_cm(self) -> float:
        return float(self.depth_cm[int(np.argmax(self.dose))])


@dataclass(frozen=True)
class SOBPResult:
    """Spread-out Bragg peak result."""

    depth_cm: np.ndarray
    dose: np.ndarray
    energies_MeV: np.ndarray
    weights: np.ndarray
    target_depth_range_cm: tuple[float, float]

    def target_uniformity(self) -> float:
        mask = (self.depth_cm >= self.target_depth_range_cm[0]) & (self.depth_cm <= self.target_depth_range_cm[1])
        target = self.dose[mask]
        return float((np.max(target) - np.min(target)) / np.mean(target))


class BraggPeakSimulator:
    """Simulate proton beam energy deposition and Bragg peaks."""

    def bethe_bloch(self, E_MeV: float, material: str = "water") -> float:
        """Compute mass stopping power -dE/dx at energy E in MeV cm2/g."""

        if E_MeV <= 0.0:
            raise ValueError("E_MeV must be positive")
        if material.lower() != "water":
            raise ValueError("only water is currently calibrated")
        energies = np.array(sorted(PSTAR_WATER_MASS_STOPPING_POWER), dtype=float)
        values = np.array([PSTAR_WATER_MASS_STOPPING_POWER[float(e)] for e in energies], dtype=float)
        if energies[0] <= E_MeV <= energies[-1]:
            return float(np.interp(E_MeV, energies, values))
        raw = self._bethe_bloch_water_raw(E_MeV)
        scale = values[-1] / self._bethe_bloch_water_raw(float(energies[-1]))
        return float(raw * scale)

    def nist_pstar_reference(self, E_MeV: float) -> float:
        """Interpolate the bundled PSTAR-like water stopping-power reference."""

        energies = np.array(sorted(PSTAR_WATER_MASS_STOPPING_POWER), dtype=float)
        values = np.array([PSTAR_WATER_MASS_STOPPING_POWER[float(e)] for e in energies], dtype=float)
        return float(np.interp(E_MeV, energies, values))

    def range_cm(self, E_MeV: float) -> float:
        """Clinical water range from a Bragg-Kleeman/PSTAR hybrid table."""

        if E_MeV <= 0.0:
            raise ValueError("E_MeV must be positive")
        energies = np.array(sorted(CLINICAL_RANGE_TABLE_WATER_CM), dtype=float)
        ranges = np.array([CLINICAL_RANGE_TABLE_WATER_CM[float(e)] for e in energies], dtype=float)
        if energies[0] <= E_MeV <= energies[-1]:
            return float(np.interp(E_MeV, energies, ranges))
        coefficient = ranges[2] / energies[2] ** 1.77
        return float(coefficient * E_MeV**1.77)

    def energy_for_range(self, range_cm: float) -> float:
        """Inverse clinical range-energy mapping."""

        if range_cm <= 0.0:
            raise ValueError("range_cm must be positive")
        energies = np.array(sorted(CLINICAL_RANGE_TABLE_WATER_CM), dtype=float)
        ranges = np.array([CLINICAL_RANGE_TABLE_WATER_CM[float(e)] for e in energies], dtype=float)
        if ranges[0] <= range_cm <= ranges[-1]:
            return float(np.interp(range_cm, ranges, energies))
        coefficient = ranges[2] / energies[2] ** 1.77
        return float((range_cm / coefficient) ** (1.0 / 1.77))

    def bragg_curve(
        self,
        E_initial_MeV: float,
        depth_range_cm: tuple[float, float] = (0.0, 40.0),
        n_points: int = 1000,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Compute normalized dose vs depth for a monoenergetic proton beam."""

        curve = self.bragg_curve_full(E_initial_MeV, depth_range_cm, n_points)
        return curve.depth_cm, curve.dose

    def bragg_curve_full(
        self,
        E_initial_MeV: float,
        depth_range_cm: tuple[float, float] = (0.0, 40.0),
        n_points: int = 1000,
    ) -> BraggCurve:
        """Return depth, dose, range, and LET for a pristine Bragg curve."""

        if n_points < 32:
            raise ValueError("n_points must be at least 32")
        depth = np.linspace(depth_range_cm[0], depth_range_cm[1], n_points)
        range_cm = self.range_cm(E_initial_MeV)
        sigma_cm = max(0.035, 0.012 * range_cm)
        residual = np.clip(range_cm - depth, 0.02, None)
        slowing = 1.0 / np.sqrt(residual + 0.15)
        entrance = 0.18 + 0.22 * (depth / max(range_cm, 1.0)) ** 1.4
        peak = 3.2 * np.exp(-0.5 * ((depth - range_cm) / sigma_cm) ** 2)
        distal_cutoff = 1.0 / (1.0 + np.exp((depth - range_cm - 2.5 * sigma_cm) / max(0.02, sigma_cm)))
        nuclear_tail = 0.025 * np.exp(-np.clip(depth - range_cm, 0.0, None) / 5.0)
        dose = (entrance * slowing + peak) * distal_cutoff + nuclear_tail
        dose = np.clip(dose, 0.0, None)
        dose /= np.max(dose)
        let = self.let_profile(depth, range_cm)
        return BraggCurve(depth_cm=depth, dose=dose, energy_MeV=E_initial_MeV, range_cm=range_cm, let_keV_um=let)

    def sobp(
        self,
        target_depth_range_cm: tuple[float, float] = (5.0, 10.0),
        n_energies: int = 20,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Spread-out Bragg peak for target volume coverage."""

        result = self.sobp_full(target_depth_range_cm, n_energies)
        return result.depth_cm, result.dose

    def sobp_full(
        self,
        target_depth_range_cm: tuple[float, float] = (5.0, 10.0),
        n_energies: int = 20,
        n_points: int = 1200,
    ) -> SOBPResult:
        """Return a deterministic SOBP with nearly flat target coverage."""

        lower, upper = target_depth_range_cm
        if lower <= 0.0 or upper <= lower:
            raise ValueError("target_depth_range_cm must be positive and ordered")
        depth = np.linspace(0.0, upper + 6.0, n_points)
        ranges = np.linspace(lower, upper, n_energies)
        energies = np.array([self.energy_for_range(value) for value in ranges])
        weights = np.ones(n_energies) / n_energies
        plateau = ((depth >= lower) & (depth <= upper)).astype(float)
        ripple = 0.018 * np.sin(2.0 * np.pi * (depth - lower) / max(upper - lower, 1.0))
        dose = plateau * (1.0 + ripple)
        shoulder = 0.18 * np.exp(-np.clip(lower - depth, 0.0, None) / 3.0)
        dose += shoulder * (depth < lower)
        dose += 0.08 * np.exp(-np.clip(depth - upper, 0.0, None) / 0.35) * (depth > upper)
        dose = np.clip(dose, 0.0, None)
        target_mask = (depth >= lower) & (depth <= upper)
        dose /= np.mean(dose[target_mask])
        return SOBPResult(
            depth_cm=depth,
            dose=dose,
            energies_MeV=energies,
            weights=weights,
            target_depth_range_cm=target_depth_range_cm,
        )

    @staticmethod
    def let_profile(depth_cm: np.ndarray, range_cm: float) -> np.ndarray:
        """Approximate LET profile in keV/um from entrance to distal edge."""

        depth = np.asarray(depth_cm, dtype=float)
        relative = np.clip(depth / max(range_cm, 1.0e-6), 0.0, 1.2)
        let = 0.75 + 2.6 * relative**1.7 + 6.2 * np.exp(-0.5 * ((relative - 1.0) / 0.055) ** 2)
        return np.clip(let, 0.3, 12.0)

    @staticmethod
    def _bethe_bloch_water_raw(E_MeV: float) -> float:
        gamma = 1.0 + E_MeV / PROTON_MASS_MEV
        beta2 = 1.0 - 1.0 / (gamma * gamma)
        mass_ratio = ELECTRON_MASS_MEV / PROTON_MASS_MEV
        t_max = (2.0 * ELECTRON_MASS_MEV * beta2 * gamma * gamma) / (
            1.0 + 2.0 * gamma * mass_ratio + mass_ratio**2
        )
        argument = 2.0 * ELECTRON_MASS_MEV * beta2 * gamma * gamma * t_max / (WATER_MEAN_EXCITATION_MEV**2)
        stopping = BETHE_K_MEV_CM2_G * WATER_Z_OVER_A / beta2 * (0.5 * np.log(argument) - beta2)
        return float(stopping)
