"""Molecular hydrogen ROS chemistry and signaling models."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from hydroonc_h2.constants import (
    GSH_PHYSIOLOGICAL_MM,
    H2_SATURATED_MM,
    K_GSH_OH_M_INV_S,
    K_H2_H2O2_M_INV_S,
    K_H2_OH_M_INV_S,
    K_H2_ONOO_M_INV_S,
    K_H2_SUPEROXIDE_M_INV_S,
)


@dataclass(frozen=True)
class ROSResult:
    """Selective ROS scavenging result for molecular hydrogen."""

    h2_concentration_mM: float
    oh_initial_uM: float
    oh_remaining_uM: float
    superoxide_initial_uM: float
    superoxide_remaining_uM: float
    h2o2_initial_uM: float
    h2o2_remaining_uM: float
    peroxynitrite_initial_uM: float
    peroxynitrite_remaining_uM: float
    k_h2_oh_M_inv_s: float
    k_h2_onoo_M_inv_s: float
    h2_oh_pseudo_first_order_s: float
    gsh_oh_pseudo_first_order_s: float
    h2_fraction_of_oh_scavenging: float
    caveat: str

    @property
    def oh_reduction_fraction(self) -> float:
        return 1.0 - self.oh_remaining_uM / self.oh_initial_uM

    @property
    def superoxide_reduction_fraction(self) -> float:
        return 1.0 - self.superoxide_remaining_uM / self.superoxide_initial_uM

    @property
    def h2o2_reduction_fraction(self) -> float:
        return 1.0 - self.h2o2_remaining_uM / self.h2o2_initial_uM


@dataclass(frozen=True)
class SignalingEffect:
    """Phenomenological H2 signaling modulation summary."""

    h2_concentration_mM: float
    nrf2_activation: float
    nfkb_suppression: float
    mitochondrial_support: float
    p53_context_score: float


class MolecularH2Model:
    """Model molecular hydrogen effects in cancer."""

    def simulate_ros_scavenging(
        self,
        h2_concentration_mM: float = H2_SATURATED_MM,
        oh_radical_concentration_uM: float = 1.0,
        superoxide_concentration_uM: float = 10.0,
        h2o2_concentration_uM: float = 50.0,
        peroxynitrite_concentration_uM: float = 1.0,
        gsh_concentration_mM: float = GSH_PHYSIOLOGICAL_MM,
        exposure_time_s: float = 1.0e-6,
    ) -> ROSResult:
        """Simulate selective ROS scavenging by H2 with competing GSH."""

        if min(
            h2_concentration_mM,
            oh_radical_concentration_uM,
            superoxide_concentration_uM,
            h2o2_concentration_uM,
            peroxynitrite_concentration_uM,
            gsh_concentration_mM,
            exposure_time_s,
        ) < 0.0:
            raise ValueError("concentrations and exposure_time_s must be non-negative")

        h2_M = h2_concentration_mM * 1.0e-3
        gsh_M = gsh_concentration_mM * 1.0e-3
        h2_oh_pseudo = K_H2_OH_M_INV_S * h2_M
        gsh_oh_pseudo = K_GSH_OH_M_INV_S * gsh_M
        total_oh_pseudo = h2_oh_pseudo + gsh_oh_pseudo
        h2_fraction = h2_oh_pseudo / total_oh_pseudo if total_oh_pseudo > 0.0 else 0.0

        total_oh_decay = 1.0 - np.exp(-total_oh_pseudo * exposure_time_s)
        oh_removed_by_h2 = oh_radical_concentration_uM * total_oh_decay * h2_fraction
        oh_remaining = max(0.0, oh_radical_concentration_uM - oh_removed_by_h2)

        onoo_decay = 1.0 - np.exp(-K_H2_ONOO_M_INV_S * h2_M * exposure_time_s)
        onoo_remaining = peroxynitrite_concentration_uM * (1.0 - onoo_decay)

        return ROSResult(
            h2_concentration_mM=float(h2_concentration_mM),
            oh_initial_uM=float(oh_radical_concentration_uM),
            oh_remaining_uM=float(oh_remaining),
            superoxide_initial_uM=float(superoxide_concentration_uM),
            superoxide_remaining_uM=float(superoxide_concentration_uM),
            h2o2_initial_uM=float(h2o2_concentration_uM),
            h2o2_remaining_uM=float(h2o2_concentration_uM),
            peroxynitrite_initial_uM=float(peroxynitrite_concentration_uM),
            peroxynitrite_remaining_uM=float(onoo_remaining),
            k_h2_oh_M_inv_s=K_H2_OH_M_INV_S,
            k_h2_onoo_M_inv_s=K_H2_ONOO_M_INV_S,
            h2_oh_pseudo_first_order_s=float(h2_oh_pseudo),
            gsh_oh_pseudo_first_order_s=float(gsh_oh_pseudo),
            h2_fraction_of_oh_scavenging=float(h2_fraction),
            caveat="At physiological H2 and GSH, stoichiometric OH scavenging by H2 is <1%; signaling is modeled as dominant.",
        )

    def signaling_modulation(self, h2_concentration_mM: float = H2_SATURATED_MM) -> SignalingEffect:
        """Return phenomenological signaling effects likely dominating biology."""

        if h2_concentration_mM < 0.0:
            raise ValueError("h2_concentration_mM must be non-negative")
        saturation = h2_concentration_mM / (0.25 + h2_concentration_mM) if h2_concentration_mM > 0.0 else 0.0
        return SignalingEffect(
            h2_concentration_mM=float(h2_concentration_mM),
            nrf2_activation=float(0.55 * saturation),
            nfkb_suppression=float(0.42 * saturation),
            mitochondrial_support=float(0.50 * saturation),
            p53_context_score=float(0.10 * saturation),
        )

    @staticmethod
    def selective_rate_constants() -> dict[str, float]:
        """Return core H2 ROS rate constants in M^-1 s^-1."""

        return {
            "H2+OH": K_H2_OH_M_INV_S,
            "H2+ONOO": K_H2_ONOO_M_INV_S,
            "H2+O2_superoxide": K_H2_SUPEROXIDE_M_INV_S,
            "H2+H2O2": K_H2_H2O2_M_INV_S,
            "GSH+OH": K_GSH_OH_M_INV_S,
        }
