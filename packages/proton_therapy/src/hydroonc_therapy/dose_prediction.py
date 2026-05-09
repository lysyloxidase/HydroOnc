"""Lightweight 3D U-Net-compatible proton dose prediction."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np


@dataclass(frozen=True)
class UNetArchitecture:
    """Compiled 3D U-Net architecture descriptor."""

    input_shape: tuple[int, int, int, int]
    output_shape: tuple[int, int, int]
    levels: int
    base_filters: int
    compiled: bool
    parameter_count: int


class DosePredictionUNet:
    """3D U-Net surface for proton dose prediction.

    The class provides the shape and inference contract of a clinical dose
    predictor without requiring a deep-learning backend in CI.
    """

    def __init__(self, levels: int = 5, base_filters: int = 32, seed: int = 17) -> None:
        self.levels = levels
        self.base_filters = base_filters
        self.rng = np.random.default_rng(seed)
        self.architecture: Optional[UNetArchitecture] = None

    def compile(
        self,
        input_shape: tuple[int, int, int, int] = (64, 64, 64, 3),
    ) -> UNetArchitecture:
        """Compile a 3D U-Net architecture descriptor."""

        if len(input_shape) != 4:
            raise ValueError("input_shape must be (x, y, z, channels)")
        if min(input_shape[:3]) < 16:
            raise ValueError("spatial dimensions must be at least 16")
        channels = input_shape[-1]
        params = 0
        in_channels = channels
        for level in range(self.levels):
            filters = self.base_filters * (2**level)
            params += (3 * 3 * 3 * in_channels * filters + filters)
            params += (3 * 3 * 3 * filters * filters + filters)
            in_channels = filters
        for level in reversed(range(self.levels - 1)):
            filters = self.base_filters * (2**level)
            params += (2 * 2 * 2 * in_channels * filters + filters)
            params += (3 * 3 * 3 * (filters * 2) * filters + filters)
            in_channels = filters
        params += in_channels + 1
        self.architecture = UNetArchitecture(
            input_shape=input_shape,
            output_shape=input_shape[:3],
            levels=self.levels,
            base_filters=self.base_filters,
            compiled=True,
            parameter_count=int(params),
        )
        return self.architecture

    def predict(self, inputs: np.ndarray) -> np.ndarray:
        """Predict a 3D dose grid from CT, beam mask, and ROI channels."""

        array = np.asarray(inputs, dtype=float)
        if array.ndim != 4:
            raise ValueError("inputs must have shape (x, y, z, channels)")
        if self.architecture is None or self.architecture.input_shape != array.shape:
            self.compile(input_shape=array.shape)
        ct = array[..., 0]
        beam = array[..., 1] if array.shape[-1] > 1 else np.ones_like(ct)
        roi = array[..., 2] if array.shape[-1] > 2 else np.ones_like(ct)
        z_axis = np.linspace(0.0, 1.0, array.shape[2])
        depth_kernel = np.exp(-0.5 * ((z_axis - 0.62) / 0.18) ** 2)
        density_modifier = np.clip(1.0 + 0.15 * (ct - np.mean(ct)), 0.65, 1.35)
        dose = 70.0 * beam * density_modifier * depth_kernel[None, None, :]
        dose *= 0.85 + 0.15 * np.clip(roi, 0.0, 1.0)
        return np.clip(dose, 0.0, None)

    def loss(
        self,
        predicted: np.ndarray,
        reference: np.ndarray,
        roi_mask: Optional[np.ndarray] = None,
        dvh_weight: float = 0.25,
    ) -> float:
        """MAE plus a DVH-aware weighted MSE term."""

        pred = np.asarray(predicted, dtype=float)
        ref = np.asarray(reference, dtype=float)
        if pred.shape != ref.shape:
            raise ValueError("predicted and reference must have the same shape")
        mae = np.mean(np.abs(pred - ref))
        if roi_mask is None:
            weights = np.ones_like(pred)
        else:
            weights = np.where(np.asarray(roi_mask, dtype=bool), 3.0, 1.0)
        mse = np.mean(weights * (pred - ref) ** 2)
        return float(mae + dvh_weight * mse)

    def mc_dropout_uncertainty(self, inputs: np.ndarray, n_samples: int = 8) -> tuple[np.ndarray, np.ndarray]:
        """Simple Monte-Carlo dropout-style uncertainty estimate."""

        base = self.predict(inputs)
        samples = []
        for _ in range(n_samples):
            noise = self.rng.normal(1.0, 0.025, size=base.shape)
            samples.append(np.clip(base * noise, 0.0, None))
        stack = np.stack(samples, axis=0)
        return np.mean(stack, axis=0), np.std(stack, axis=0)
