"""Dose-volume histogram utilities."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class DVHCurve:
    """Cumulative dose-volume histogram."""

    dose_bins_Gy: np.ndarray
    volume_fraction: np.ndarray
    structure_name: str

    def d_at_volume(self, volume_fraction: float) -> float:
        """Return dose received by at least the requested volume fraction."""

        if not 0.0 <= volume_fraction <= 1.0:
            raise ValueError("volume_fraction must be in [0, 1]")
        order = np.argsort(self.volume_fraction)
        return float(np.interp(volume_fraction, self.volume_fraction[order], self.dose_bins_Gy[order]))

    def v_at_dose(self, dose_Gy: float) -> float:
        """Return volume fraction receiving at least dose_Gy."""

        return float(np.interp(dose_Gy, self.dose_bins_Gy, self.volume_fraction))


def compute_dvh(
    dose_grid_Gy: np.ndarray,
    structure_mask: np.ndarray,
    n_bins: int = 100,
    structure_name: str = "structure",
) -> DVHCurve:
    """Compute a cumulative DVH from a 3D dose grid and Boolean mask."""

    dose = np.asarray(dose_grid_Gy, dtype=float)
    mask = np.asarray(structure_mask, dtype=bool)
    if dose.shape != mask.shape:
        raise ValueError("dose_grid_Gy and structure_mask must have the same shape")
    values = dose[mask]
    if values.size == 0:
        raise ValueError("structure_mask selects no voxels")
    max_dose = float(np.max(values))
    bins = np.linspace(0.0, max_dose, n_bins)
    volume_fraction = np.array([np.mean(values >= value) for value in bins], dtype=float)
    return DVHCurve(dose_bins_Gy=bins, volume_fraction=volume_fraction, structure_name=structure_name)
