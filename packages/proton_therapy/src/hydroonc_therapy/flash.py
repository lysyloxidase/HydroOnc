"""FLASH proton therapy dose-rate and oxygen-response models."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from hydroonc_therapy.constants import FLASH_DOSE_RATE_THRESHOLD_GY_S


@dataclass(frozen=True)
class FLASHResponse:
    """Predicted FLASH biological response summary."""

    dose_Gy: float
    delivery_time_s: float
    dose_rate_Gy_s: float
    is_flash: bool
    oxygen_before_mmHg: float
    oxygen_after_mmHg: float
    normal_tissue_sparing_factor: float


class FLASHProtonTherapy:
    """FLASH proton therapy simulation above 40 Gy/s."""

    def __init__(self, threshold_Gy_s: float = FLASH_DOSE_RATE_THRESHOLD_GY_S) -> None:
        self.threshold_Gy_s = threshold_Gy_s
        self.fast01 = {
            "trial": "FAST-01",
            "patients": 10,
            "lesions": 12,
            "single_fraction_Gy": 8.0,
            "pain_response_fraction": 0.67,
        }

    def dose_rate(self, dose_Gy: float, delivery_time_s: float) -> float:
        """Return dose rate in Gy/s."""

        if dose_Gy < 0.0:
            raise ValueError("dose_Gy must be non-negative")
        if delivery_time_s <= 0.0:
            raise ValueError("delivery_time_s must be positive")
        return float(dose_Gy / delivery_time_s)

    def is_flash(self, dose_Gy: float, delivery_time_s: float) -> bool:
        """Return whether delivery meets the FLASH threshold."""

        return self.dose_rate(dose_Gy, delivery_time_s) > self.threshold_Gy_s

    def oxygen_depletion(
        self,
        dose_Gy: float,
        delivery_time_s: float,
        oxygen_initial_mmHg: float = 30.0,
    ) -> float:
        """Estimate post-irradiation oxygen tension under rapid delivery."""

        rate = self.dose_rate(dose_Gy, delivery_time_s)
        depletion_per_Gy = 0.42 * (rate / (rate + self.threshold_Gy_s))
        return float(max(0.0, oxygen_initial_mmHg - depletion_per_Gy * dose_Gy))

    def response(
        self,
        dose_Gy: float = 8.0,
        delivery_time_s: float = 0.1,
        oxygen_initial_mmHg: float = 30.0,
    ) -> FLASHResponse:
        """Return a compact FLASH response summary."""

        rate = self.dose_rate(dose_Gy, delivery_time_s)
        oxygen_after = self.oxygen_depletion(dose_Gy, delivery_time_s, oxygen_initial_mmHg)
        flash = rate > self.threshold_Gy_s
        oxygen_ratio = oxygen_after / max(oxygen_initial_mmHg, 1.0e-9)
        sparing = 1.0 - (0.18 * (1.0 - oxygen_ratio) if flash else 0.0)
        return FLASHResponse(
            dose_Gy=float(dose_Gy),
            delivery_time_s=float(delivery_time_s),
            dose_rate_Gy_s=rate,
            is_flash=flash,
            oxygen_before_mmHg=float(oxygen_initial_mmHg),
            oxygen_after_mmHg=oxygen_after,
            normal_tissue_sparing_factor=float(np.clip(sparing, 0.70, 1.0)),
        )

    def ros_recombination_factor(self, dose_rate_Gy_s: float) -> float:
        """Phenomenological radical-radical recombination factor."""

        if dose_rate_Gy_s < 0.0:
            raise ValueError("dose_rate_Gy_s must be non-negative")
        return float(1.0 + 0.35 * dose_rate_Gy_s / (dose_rate_Gy_s + self.threshold_Gy_s))
