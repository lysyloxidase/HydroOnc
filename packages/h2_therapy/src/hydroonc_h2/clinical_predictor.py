"""Research-grade clinical response predictor for H2 therapy."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Union

import numpy as np

from hydroonc_h2.constants import RESEARCH_GRADE_SAMPLE_SIZE_UPPER_BOUND


@dataclass(frozen=True)
class H2PatientFeatures:
    """Multi-omic features for H2 response estimation."""

    pd1_cd8_percent: float
    tumor_ph_e: float
    nrf2_keap1_activity: float
    coq10_index: float
    complex_i_activity: float
    tmb_mut_per_mb: float
    tumor_type: str
    stage: str


@dataclass(frozen=True)
class H2ResponsePrediction:
    """Research-grade response prediction."""

    probability_pfs_benefit: float
    risk_group: str
    caveat: str
    feature_score: float


class H2ResponsePredictor:
    """ML-style predictor for H2 therapy response.

    The model is intentionally transparent because published H2 oncology data
    remain too small for clinical-grade prediction.
    """

    TUMOR_TYPE_PRIORS = {
        "nsclc": 0.08,
        "lung": 0.08,
        "crc": 0.04,
        "colorectal": 0.04,
        "glioma": 0.02,
        "brain_mets": 0.03,
    }

    def predict(self, features: H2PatientFeatures) -> H2ResponsePrediction:
        """Estimate probability of >3 month PFS benefit."""

        score = self._linear_score(features)
        probability = 1.0 / (1.0 + np.exp(-score))
        if probability >= 0.66:
            group = "high"
        elif probability >= 0.40:
            group = "intermediate"
        else:
            group = "low"
        return H2ResponsePrediction(
            probability_pfs_benefit=float(probability),
            risk_group=group,
            caveat=f"RESEARCH-GRADE ONLY: published H2 oncology n < {RESEARCH_GRADE_SAMPLE_SIZE_UPPER_BOUND}.",
            feature_score=float(score),
        )

    def predict_from_dict(self, features: dict[str, Union[float, str]]) -> H2ResponsePrediction:
        """Predict from a dictionary of feature values."""

        return self.predict(
            H2PatientFeatures(
                pd1_cd8_percent=float(features.get("pd1_cd8_percent", 35.0)),
                tumor_ph_e=float(features.get("tumor_ph_e", 6.7)),
                nrf2_keap1_activity=float(features.get("nrf2_keap1_activity", 0.5)),
                coq10_index=float(features.get("coq10_index", 0.5)),
                complex_i_activity=float(features.get("complex_i_activity", 0.5)),
                tmb_mut_per_mb=float(features.get("tmb_mut_per_mb", 8.0)),
                tumor_type=str(features.get("tumor_type", "nsclc")),
                stage=str(features.get("stage", "IV")),
            )
        )

    def _linear_score(self, features: H2PatientFeatures) -> float:
        pd1 = np.clip(features.pd1_cd8_percent / 100.0, 0.0, 1.0)
        acid = np.clip((7.2 - features.tumor_ph_e) / 0.8, 0.0, 1.5)
        nrf2 = np.clip(features.nrf2_keap1_activity, 0.0, 1.5)
        mito = 0.5 * np.clip(features.coq10_index, 0.0, 1.5) + 0.5 * np.clip(features.complex_i_activity, 0.0, 1.5)
        tmb = np.clip(features.tmb_mut_per_mb / 20.0, 0.0, 1.5)
        tumor_prior = self.TUMOR_TYPE_PRIORS.get(features.tumor_type.lower(), 0.0)
        stage_penalty = -0.08 if features.stage.upper() in {"IV", "4"} else 0.02
        return float(-0.65 + 1.10 * pd1 + 0.42 * acid + 0.34 * nrf2 + 0.38 * mito + 0.22 * tmb + tumor_prior + stage_penalty)
