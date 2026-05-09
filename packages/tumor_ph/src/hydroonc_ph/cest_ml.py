"""CEST-MRI pH map reconstruction with lightweight NumPy models."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np


@dataclass(frozen=True)
class CESTArchitecture:
    """Compiled CEST pH mapper architecture."""

    input_features: int
    hidden_units: int
    residual_blocks: int
    compiled: bool
    weights: tuple[np.ndarray, ...]
    biases: tuple[np.ndarray, ...]


class CESTpHMapper:
    """Deep-learning style CEST-MRI pH map reconstruction.

    This class avoids a hard dependency on PyTorch/JAX for CI. It provides a
    deterministic residual MLP surface over sparse Z-spectra and can be swapped
    for a heavier backend later without changing the public API.
    """

    def __init__(
        self,
        ph_range: tuple[float, float] = (6.0, 7.5),
        physics_penalty_weight: float = 0.05,
        seed: int = 13,
    ) -> None:
        self.ph_range = ph_range
        self.physics_penalty_weight = physics_penalty_weight
        self.rng = np.random.default_rng(seed)
        self.architecture: Optional[CESTArchitecture] = None

    def compile(
        self,
        input_features: int = 64,
        hidden_units: int = 32,
        residual_blocks: int = 2,
    ) -> CESTArchitecture:
        """Compile a residual MLP architecture for Z-spectrum inputs."""

        if input_features < 4:
            raise ValueError("input_features must be at least 4")
        if hidden_units < 4:
            raise ValueError("hidden_units must be at least 4")
        weights = [
            self.rng.normal(0.0, 0.08, size=(input_features, hidden_units)),
            self.rng.normal(0.0, 0.05, size=(hidden_units, hidden_units)),
            self.rng.normal(0.0, 0.05, size=(hidden_units, hidden_units)),
            self.rng.normal(0.0, 0.04, size=(hidden_units, 1)),
        ]
        biases = [
            np.zeros(hidden_units),
            np.zeros(hidden_units),
            np.zeros(hidden_units),
            np.zeros(1),
        ]
        self.architecture = CESTArchitecture(
            input_features=input_features,
            hidden_units=hidden_units,
            residual_blocks=residual_blocks,
            compiled=True,
            weights=tuple(weights),
            biases=tuple(biases),
        )
        return self.architecture

    def predict(self, z_spectrum: np.ndarray) -> np.ndarray:
        """Predict extracellular pH from one or more CEST Z-spectra."""

        z = np.asarray(z_spectrum, dtype=float)
        if z.ndim == 1:
            z = z[None, :]
        if z.ndim != 2:
            raise ValueError("z_spectrum must be 1D or 2D")
        if self.architecture is None or self.architecture.input_features != z.shape[1]:
            self.compile(input_features=z.shape[1])

        arch = self.architecture
        assert arch is not None
        z_norm = np.clip(z, 0.0, 1.2)
        hidden = np.tanh(_dense(z_norm, arch.weights[0]) + arch.biases[0])
        for block in range(arch.residual_blocks):
            w = arch.weights[1 + (block % 2)]
            b = arch.biases[1 + (block % 2)]
            hidden = np.tanh(_dense(hidden, w) + b) + 0.35 * hidden
        learned_residual = (_dense(hidden, arch.weights[-1]) + arch.biases[-1]).reshape(-1)

        cest_contrast = np.mean(1.0 - z_norm, axis=1)
        asymmetry = np.mean(z_norm[:, : z.shape[1] // 2], axis=1) - np.mean(z_norm[:, z.shape[1] // 2 :], axis=1)
        ph = 7.45 - 2.15 * cest_contrast + 0.25 * asymmetry + 0.04 * learned_residual
        return np.clip(ph, self.ph_range[0], self.ph_range[1])

    def physics_penalty(self, predicted_ph: np.ndarray) -> float:
        """Penalty for pH predictions outside the acidoCEST measurement range."""

        predicted = np.asarray(predicted_ph, dtype=float)
        low, high = self.ph_range
        below = np.clip(low - predicted, 0.0, None)
        above = np.clip(predicted - high, 0.0, None)
        return float(self.physics_penalty_weight * np.mean(below**2 + above**2))

    def fit_linear_calibration(self, z_spectra: np.ndarray, reference_ph: np.ndarray) -> dict[str, float]:
        """Fit a simple CEST-contrast-to-pH calibration for reporting."""

        z = np.asarray(z_spectra, dtype=float)
        y = np.asarray(reference_ph, dtype=float)
        if z.ndim != 2:
            raise ValueError("z_spectra must be a 2D array")
        if y.shape[0] != z.shape[0]:
            raise ValueError("reference_ph length must match number of spectra")
        contrast = np.mean(1.0 - np.clip(z, 0.0, 1.2), axis=1)
        slope, intercept = np.polyfit(contrast, y, deg=1)
        predictions = slope * contrast + intercept
        rmse = float(np.sqrt(np.mean((predictions - y) ** 2)))
        return {"slope": float(slope), "intercept": float(intercept), "rmse": rmse}

    @staticmethod
    def overlay_tumor_mask(ph_map: np.ndarray, tumor_mask: np.ndarray) -> np.ndarray:
        """Return pH values inside a tumor segmentation mask, NaN outside."""

        ph = np.asarray(ph_map, dtype=float)
        mask = np.asarray(tumor_mask, dtype=bool)
        if ph.shape != mask.shape:
            raise ValueError("ph_map and tumor_mask must have the same shape")
        return np.where(mask, ph, np.nan)


def _dense(inputs: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Small deterministic dense layer avoiding platform BLAS warning noise."""

    return np.sum(inputs[:, :, None] * weights[None, :, :], axis=1)
