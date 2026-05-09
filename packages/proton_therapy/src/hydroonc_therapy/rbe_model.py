"""Variable-RBE models for proton therapy."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from hydroonc_therapy.constants import CLINICAL_RBE


@dataclass(frozen=True)
class RBESample:
    """One variable-RBE model evaluation."""

    model: str
    let_keV_um: float
    alpha_beta_Gy: float
    dose_Gy: float
    rbe: float


class VariableRBE:
    """Variable RBE models for proton therapy."""

    def mcnamara(self, LET_keV_um: float, alpha_beta_Gy: float = 3.0, dose_Gy: float = 2.0) -> float:
        del dose_Gy
        self._validate(LET_keV_um, alpha_beta_Gy)
        return float(clinical_safe_min(1.1 + 0.045 * LET_keV_um / alpha_beta_Gy))

    def wedenberg(self, LET_keV_um: float, alpha_beta_Gy: float = 3.0, dose_Gy: float = 2.0) -> float:
        del dose_Gy
        self._validate(LET_keV_um, alpha_beta_Gy)
        return float(clinical_safe_min(1.0 + 0.12 * LET_keV_um / alpha_beta_Gy))

    def mcmahon(self, LET_keV_um: float, alpha_beta_Gy: float = 3.0, dose_Gy: float = 2.0) -> float:
        self._validate(LET_keV_um, alpha_beta_Gy)
        repair_factor = 1.0 / (1.0 + dose_Gy / (alpha_beta_Gy + 1.0))
        return float(clinical_safe_min(1.03 + 0.22 * np.log1p(LET_keV_um) / alpha_beta_Gy + 0.09 * repair_factor))

    def carabe_fernandez(self, LET_keV_um: float, alpha_beta_Gy: float = 3.0, dose_Gy: float = 2.0) -> float:
        self._validate(LET_keV_um, alpha_beta_Gy)
        dose_term = 0.015 * LET_keV_um / (dose_Gy + 1.0)
        return float(clinical_safe_min(1.02 + 0.10 * LET_keV_um / alpha_beta_Gy + dose_term))

    def evaluate(
        self,
        LET_keV_um: float,
        alpha_beta_Gy: float = 3.0,
        dose_Gy: float = 2.0,
        model: str = "mcnamara",
    ) -> RBESample:
        """Evaluate one named RBE model."""

        key = model.lower().replace("-", "_")
        if key not in {"mcnamara", "wedenberg", "mcmahon", "carabe_fernandez", "carabe"}:
            raise ValueError("unknown RBE model")
        method = self.carabe_fernandez if key in {"carabe_fernandez", "carabe"} else getattr(self, key)
        rbe = method(LET_keV_um, alpha_beta_Gy=alpha_beta_Gy, dose_Gy=dose_Gy)
        return RBESample(key, LET_keV_um, alpha_beta_Gy, dose_Gy, rbe)

    def all_models(self, LET_keV_um: float, alpha_beta_Gy: float = 3.0, dose_Gy: float = 2.0) -> dict[str, float]:
        """Return RBE values from all implemented models."""

        return {
            "mcnamara": self.mcnamara(LET_keV_um, alpha_beta_Gy, dose_Gy),
            "wedenberg": self.wedenberg(LET_keV_um, alpha_beta_Gy, dose_Gy),
            "mcmahon": self.mcmahon(LET_keV_um, alpha_beta_Gy, dose_Gy),
            "carabe_fernandez": self.carabe_fernandez(LET_keV_um, alpha_beta_Gy, dose_Gy),
        }

    def let_weighted_rbe(self, dose: np.ndarray, let_keV_um: np.ndarray, alpha_beta_Gy: float = 3.0) -> float:
        """Dose-weighted RBE for OAR-sparing optimization."""

        dose = np.asarray(dose, dtype=float)
        let_values = np.asarray(let_keV_um, dtype=float)
        if dose.shape != let_values.shape:
            raise ValueError("dose and LET arrays must have the same shape")
        weights = np.clip(dose, 0.0, None)
        if np.sum(weights) == 0.0:
            return CLINICAL_RBE
        rbe = np.vectorize(lambda value: self.mcnamara(value, alpha_beta_Gy))(let_values)
        return float(np.sum(weights * rbe) / np.sum(weights))

    @staticmethod
    def _validate(LET_keV_um: float, alpha_beta_Gy: float) -> None:
        if LET_keV_um < 0.0:
            raise ValueError("LET must be non-negative")
        if alpha_beta_Gy <= 0.0:
            raise ValueError("alpha_beta_Gy must be positive")


def clinical_safe_min(value: float) -> float:
    """Keep variable models at or above the clinical convention."""

    return max(CLINICAL_RBE, value)
