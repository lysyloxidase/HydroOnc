"""Immune-restoration ODE model for molecular hydrogen therapy."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from hydroonc_h2.constants import H2_SATURATED_MM


@dataclass(frozen=True)
class ImmuneTrajectory:
    """Time course for effector/exhausted T cells and tumor burden."""

    time_days: np.ndarray
    t_eff: np.ndarray
    t_exh: np.ndarray
    tumor: np.ndarray
    ifn_gamma: np.ndarray
    pd1_exhaustion_index: np.ndarray
    h2_concentration_mM: float
    anti_pd1: bool

    @property
    def final_t_eff(self) -> float:
        return float(self.t_eff[-1])

    @property
    def final_t_exh(self) -> float:
        return float(self.t_exh[-1])

    @property
    def final_tumor(self) -> float:
        return float(self.tumor[-1])

    @property
    def tumor_shrinkage_fraction(self) -> float:
        return float((self.tumor[0] - self.tumor[-1]) / self.tumor[0])


@dataclass(frozen=True)
class TherapyEndpoint:
    """Endpoint comparison for one H2/checkpoint intervention."""

    h2_concentration_mM: float
    anti_pd1: bool
    final_t_eff: float
    final_t_exh: float
    final_tumor: float
    benefit_score: float


class H2ImmuneModel:
    """Model H2 effects on anti-tumor immunity in acidic TME."""

    def __init__(
        self,
        k_prolif: float = 0.095,
        k_exhaust: float = 0.12,
        k_death_eff: float = 0.035,
        k_death_exh: float = 0.028,
        k_h2_rescue: float = 0.34,
        k_grow: float = 0.045,
        k_kill: float = 0.060,
        carrying_capacity: float = 4.0,
    ) -> None:
        self.k_prolif = k_prolif
        self.k_exhaust = k_exhaust
        self.k_death_eff = k_death_eff
        self.k_death_exh = k_death_exh
        self.k_h2_rescue = k_h2_rescue
        self.k_grow = k_grow
        self.k_kill = k_kill
        self.carrying_capacity = carrying_capacity

    def simulate(
        self,
        h2_concentration_mM: float = H2_SATURATED_MM,
        anti_pd1: bool = False,
        ph_e: float = 6.7,
        antigen: float = 1.0,
        t_max_days: float = 30.0,
        dt_days: float = 0.05,
        initial_t_eff: float = 1.0,
        initial_t_exh: float = 0.8,
        initial_tumor: float = 1.0,
    ) -> ImmuneTrajectory:
        """Integrate the T-cell exhaustion/rescue/tumor ODE system."""

        if min(h2_concentration_mM, antigen, t_max_days, dt_days, initial_t_eff, initial_t_exh, initial_tumor) < 0.0:
            raise ValueError("model inputs must be non-negative")
        if dt_days == 0.0 or t_max_days == 0.0:
            raise ValueError("dt_days and t_max_days must be positive")

        steps = int(np.floor(t_max_days / dt_days)) + 1
        time = np.linspace(0.0, dt_days * (steps - 1), steps)
        t_eff = np.empty(steps, dtype=float)
        t_exh = np.empty(steps, dtype=float)
        tumor = np.empty(steps, dtype=float)
        t_eff[0] = initial_t_eff
        t_exh[0] = initial_t_exh
        tumor[0] = initial_tumor

        for idx in range(1, steps):
            eff, exh, burden = t_eff[idx - 1], t_exh[idx - 1], tumor[idx - 1]
            d_eff, d_exh, d_tumor = self._derivatives(
                eff,
                exh,
                burden,
                h2_concentration_mM,
                anti_pd1,
                ph_e,
                antigen,
            )
            t_eff[idx] = max(0.0, eff + dt_days * d_eff)
            t_exh[idx] = max(0.0, exh + dt_days * d_exh)
            tumor[idx] = max(0.0, burden + dt_days * d_tumor)

        h2_signal = self._h2_saturation(h2_concentration_mM)
        ifn_gamma = 0.55 * t_eff * (1.0 + 0.25 * h2_signal + (0.12 if anti_pd1 else 0.0))
        pd1_index = t_exh / np.maximum(t_eff + t_exh, 1.0e-12)
        return ImmuneTrajectory(
            time_days=time,
            t_eff=t_eff,
            t_exh=t_exh,
            tumor=tumor,
            ifn_gamma=ifn_gamma,
            pd1_exhaustion_index=pd1_index,
            h2_concentration_mM=float(h2_concentration_mM),
            anti_pd1=bool(anti_pd1),
        )

    def endpoint(
        self,
        h2_concentration_mM: float = H2_SATURATED_MM,
        anti_pd1: bool = False,
        ph_e: float = 6.7,
    ) -> TherapyEndpoint:
        """Return a compact endpoint summary."""

        trajectory = self.simulate(h2_concentration_mM=h2_concentration_mM, anti_pd1=anti_pd1, ph_e=ph_e)
        benefit = trajectory.final_t_eff / trajectory.t_eff[0] + trajectory.tumor_shrinkage_fraction - trajectory.final_t_exh / 10.0
        return TherapyEndpoint(
            h2_concentration_mM=float(h2_concentration_mM),
            anti_pd1=bool(anti_pd1),
            final_t_eff=trajectory.final_t_eff,
            final_t_exh=trajectory.final_t_exh,
            final_tumor=trajectory.final_tumor,
            benefit_score=float(benefit),
        )

    def dose_response(self, h2_concentrations_mM: np.ndarray, anti_pd1: bool = False) -> list[TherapyEndpoint]:
        """Evaluate monotonic H2 dose response across concentrations."""

        values = np.asarray(h2_concentrations_mM, dtype=float)
        return [self.endpoint(float(value), anti_pd1=anti_pd1) for value in values]

    def synergy_score(self, h2_concentration_mM: float = H2_SATURATED_MM) -> dict[str, float]:
        """Compare H2, anti-PD-1, and combination benefit."""

        control = self.endpoint(0.0, anti_pd1=False)
        h2 = self.endpoint(h2_concentration_mM, anti_pd1=False)
        pd1 = self.endpoint(0.0, anti_pd1=True)
        combo = self.endpoint(h2_concentration_mM, anti_pd1=True)
        return {
            "control_final_tumor": control.final_tumor,
            "h2_final_tumor": h2.final_tumor,
            "anti_pd1_final_tumor": pd1.final_tumor,
            "combo_final_tumor": combo.final_tumor,
            "h2_benefit": control.final_tumor - h2.final_tumor,
            "anti_pd1_benefit": control.final_tumor - pd1.final_tumor,
            "combo_benefit": control.final_tumor - combo.final_tumor,
            "synergy_margin": min(h2.final_tumor, pd1.final_tumor) - combo.final_tumor,
        }

    def _derivatives(
        self,
        t_eff: float,
        t_exh: float,
        tumor: float,
        h2_concentration_mM: float,
        anti_pd1: bool,
        ph_e: float,
        antigen: float,
    ) -> tuple[float, float, float]:
        acid_stress = np.clip((6.9 - ph_e) / 0.45, 0.0, 1.5)
        h2_signal = self._h2_saturation(h2_concentration_mM)
        checkpoint_relief = 0.36 if anti_pd1 else 0.0
        rescue = self.k_h2_rescue * h2_signal + (0.085 if anti_pd1 else 0.0) + (0.075 * h2_signal if anti_pd1 else 0.0)
        exhaust = self.k_exhaust * acid_stress * (1.0 - checkpoint_relief)
        proliferation = self.k_prolif * antigen * (1.0 + 0.22 * h2_signal + (0.13 if anti_pd1 else 0.0))
        d_t_eff = proliferation * t_eff * (1.0 - t_eff / 9.0) - exhaust * t_eff - self.k_death_eff * t_eff + rescue * t_exh
        d_t_exh = exhaust * t_eff - rescue * t_exh - self.k_death_exh * t_exh
        immune_kill = self.k_kill * t_eff * tumor / (0.7 + tumor)
        growth = self.k_grow * tumor * (1.0 - tumor / self.carrying_capacity)
        d_tumor = growth - immune_kill
        return float(d_t_eff), float(d_t_exh), float(d_tumor)

    @staticmethod
    def _h2_saturation(h2_concentration_mM: float) -> float:
        if h2_concentration_mM <= 0.0:
            return 0.0
        return float(h2_concentration_mM / (0.20 + h2_concentration_mM))
