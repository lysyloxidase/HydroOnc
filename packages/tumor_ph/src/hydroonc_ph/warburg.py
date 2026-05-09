"""Warburg-effect tumor acidosis and pH-regulator dynamics."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Optional, Union

import numpy as np

from hydroonc_ph.constants import (
    DEFAULT_REGULATORS,
    INHIBITOR_ALIASES,
    INHIBITOR_IC50_UM,
    NORMAL_EXTRACELLULAR_PH,
)


@dataclass(frozen=True)
class RegulatorState:
    """Normalized pH-regulator activities."""

    v_atpase: float = 1.0
    nhe1: float = 1.0
    mct1: float = 1.0
    mct4: float = 1.0
    caix: float = 1.0
    nbce1: float = 1.0
    nbcn1: float = 1.0
    ae2: float = 1.0

    def as_dict(self) -> dict[str, float]:
        return {
            "v_atpase": self.v_atpase,
            "nhe1": self.nhe1,
            "mct1": self.mct1,
            "mct4": self.mct4,
            "caix": self.caix,
            "nbce1": self.nbce1,
            "nbcn1": self.nbcn1,
            "ae2": self.ae2,
        }


@dataclass(frozen=True)
class pHTimeSeries:
    """Tumor extracellular and intracellular pH trajectory."""

    time_s: np.ndarray
    ph_e: np.ndarray
    ph_i: np.ndarray
    lactate_mM: np.ndarray
    acid_source_rate: np.ndarray
    regulator_fluxes: dict[str, np.ndarray]

    @property
    def final_ph_e(self) -> float:
        return float(self.ph_e[-1])

    @property
    def final_ph_i(self) -> float:
        return float(self.ph_i[-1])

    @property
    def extracellular_drop(self) -> float:
        return float(self.ph_e[0] - self.ph_e[-1])


@dataclass(frozen=True)
class InhibitorResponse:
    """Predicted pH response to regulator inhibition."""

    target: str
    concentration_uM: float
    inhibition_fraction: float
    baseline_final_ph_e: float
    inhibited_final_ph_e: float
    baseline_final_ph_i: float
    inhibited_final_ph_i: float

    @property
    def delta_ph_e(self) -> float:
        return self.inhibited_final_ph_e - self.baseline_final_ph_e

    @property
    def delta_ph_i(self) -> float:
        return self.inhibited_final_ph_i - self.baseline_final_ph_i

    def as_dict(self) -> dict[str, Union[float, str]]:
        return {
            "target": self.target,
            "concentration_uM": self.concentration_uM,
            "inhibition_fraction": self.inhibition_fraction,
            "baseline_final_ph_e": self.baseline_final_ph_e,
            "inhibited_final_ph_e": self.inhibited_final_ph_e,
            "delta_ph_e": self.delta_ph_e,
            "baseline_final_ph_i": self.baseline_final_ph_i,
            "inhibited_final_ph_i": self.inhibited_final_ph_i,
            "delta_ph_i": self.delta_ph_i,
        }


class WarburgModel:
    """Model the Warburg effect and tumor acidosis.

    The model is intentionally low-dimensional but calibrated to the Phase 3
    biology gates: aerobic glycolysis drives extracellular pH from 7.4 toward
    about 6.7 over one hour, while transporters maintain a reversed pH gradient.
    """

    def __init__(
        self,
        initial_ph_e: float = NORMAL_EXTRACELLULAR_PH,
        initial_ph_i: float = 7.15,
        extracellular_tau_s: float = 950.0,
        intracellular_tau_s: float = 700.0,
    ) -> None:
        self.initial_ph_e = initial_ph_e
        self.initial_ph_i = initial_ph_i
        self.extracellular_tau_s = extracellular_tau_s
        self.intracellular_tau_s = intracellular_tau_s

    def simulate_ph_dynamics(
        self,
        cell_count: int = 10000,
        glycolysis_rate: float = 1.0,
        regulators: Optional[Mapping[str, float]] = None,
        dt_s: float = 1.0,
        t_max_s: float = 3600.0,
    ) -> pHTimeSeries:
        """Simulate extracellular and intracellular pH evolution."""

        if cell_count <= 0:
            raise ValueError("cell_count must be positive")
        if glycolysis_rate < 0.0:
            raise ValueError("glycolysis_rate must be non-negative")
        if dt_s <= 0.0 or t_max_s <= 0.0:
            raise ValueError("dt_s and t_max_s must be positive")

        state = self._regulators(regulators)
        steps = int(np.floor(t_max_s / dt_s)) + 1
        time_s = np.linspace(0.0, dt_s * (steps - 1), steps)
        ph_e = np.empty(steps, dtype=float)
        ph_i = np.empty(steps, dtype=float)
        lactate = np.empty(steps, dtype=float)
        acid = np.empty(steps, dtype=float)
        ph_e[0] = self.initial_ph_e
        ph_i[0] = self.initial_ph_i
        lactate[0] = 1.2

        target_e = self._extracellular_target(cell_count, glycolysis_rate, state)
        target_i = self._intracellular_target(glycolysis_rate, state)
        lactate_target = 1.2 + 16.0 * glycolysis_rate * (0.55 + 0.45 * state["mct4"])
        acid_rate = glycolysis_rate * (cell_count / 10000.0)
        e_factor = 1.0 - np.exp(-dt_s / self.extracellular_tau_s)
        i_factor = 1.0 - np.exp(-dt_s / self.intracellular_tau_s)
        lactate_factor = 1.0 - np.exp(-dt_s / 1200.0)

        for idx in range(1, steps):
            ph_e[idx] = ph_e[idx - 1] + (target_e - ph_e[idx - 1]) * e_factor
            ph_i[idx] = ph_i[idx - 1] + (target_i - ph_i[idx - 1]) * i_factor
            lactate[idx] = lactate[idx - 1] + (lactate_target - lactate[idx - 1]) * lactate_factor
        acid.fill(acid_rate)

        regulator_fluxes = self._regulator_fluxes(time_s, glycolysis_rate, state)
        return pHTimeSeries(
            time_s=time_s,
            ph_e=ph_e,
            ph_i=ph_i,
            lactate_mM=lactate,
            acid_source_rate=acid,
            regulator_fluxes=regulator_fluxes,
        )

    def inhibitor_effect(self, target: str, concentration_uM: float) -> dict[str, Union[float, str]]:
        """Predict pH change from inhibiting a specific pH regulator."""

        if concentration_uM < 0.0:
            raise ValueError("concentration_uM must be non-negative")
        canonical = self._canonical_target(target)
        ic50 = INHIBITOR_IC50_UM[canonical]
        inhibition = 0.0 if concentration_uM == 0.0 else concentration_uM / (concentration_uM + ic50)
        regulators = DEFAULT_REGULATORS.copy()
        regulators[canonical] = 1.0 - inhibition
        baseline = self.simulate_ph_dynamics()
        inhibited = self.simulate_ph_dynamics(regulators=regulators)
        response = InhibitorResponse(
            target=canonical,
            concentration_uM=float(concentration_uM),
            inhibition_fraction=float(inhibition),
            baseline_final_ph_e=baseline.final_ph_e,
            inhibited_final_ph_e=inhibited.final_ph_e,
            baseline_final_ph_i=baseline.final_ph_i,
            inhibited_final_ph_i=inhibited.final_ph_i,
        )
        return response.as_dict()

    @staticmethod
    def _regulators(regulators: Optional[Mapping[str, float]]) -> dict[str, float]:
        state = DEFAULT_REGULATORS.copy()
        if regulators:
            for key, value in regulators.items():
                canonical = WarburgModel._canonical_target(key, allow_buffer=True)
                state[canonical] = float(np.clip(value, 0.0, 2.0))
        return state

    @staticmethod
    def _canonical_target(target: str, allow_buffer: bool = False) -> str:
        key = target.lower().replace("_", "-").strip()
        if key in INHIBITOR_ALIASES:
            return INHIBITOR_ALIASES[key]
        normalized = key.replace("-", "_")
        if allow_buffer and normalized in DEFAULT_REGULATORS:
            return normalized
        if normalized in INHIBITOR_IC50_UM:
            return normalized
        raise ValueError(f"unknown pH regulator target: {target}")

    @staticmethod
    def _extracellular_target(
        cell_count: int,
        glycolysis_rate: float,
        state: Mapping[str, float],
    ) -> float:
        cell_scale = np.clip(np.sqrt(cell_count / 10000.0), 0.4, 2.2)
        export_activity = (
            0.30 * state["v_atpase"]
            + 0.25 * state["nhe1"]
            + 0.10 * state["mct1"]
            + 0.20 * state["mct4"]
            + 0.15 * state["caix"]
        )
        caix_acidification = 0.18 * state["caix"]
        buffer_relief = 0.05 * (state["nbce1"] + state["nbcn1"] - 2.0)
        acidification = glycolysis_rate * cell_scale * (0.50 * export_activity + caix_acidification)
        return float(np.clip(NORMAL_EXTRACELLULAR_PH - acidification + buffer_relief, 5.8, 7.6))

    @staticmethod
    def _intracellular_target(
        glycolysis_rate: float,
        state: Mapping[str, float],
    ) -> float:
        alkaline_defense = (
            0.27 * (state["nhe1"] - 1.0)
            + 0.08 * (state["nbce1"] - 1.0)
            + 0.05 * (state["nbcn1"] - 1.0)
            + 0.04 * (state["v_atpase"] - 1.0)
            + 0.04 * (state["mct4"] - 1.0)
            - 0.04 * (state["ae2"] - 1.0)
        )
        metabolic_acid = 0.08 * (glycolysis_rate - 1.0)
        return float(np.clip(7.12 + alkaline_defense - metabolic_acid, 6.5, 7.5))

    @staticmethod
    def _regulator_fluxes(
        time_s: np.ndarray,
        glycolysis_rate: float,
        state: Mapping[str, float],
    ) -> dict[str, np.ndarray]:
        envelope = 1.0 - np.exp(-time_s / 600.0)
        return {
            "v_atpase": 0.30 * glycolysis_rate * state["v_atpase"] * envelope,
            "nhe1": 0.24 * glycolysis_rate * state["nhe1"] * envelope,
            "mct_lactate_h": 0.32 * glycolysis_rate * (state["mct1"] + state["mct4"]) / 2.0 * envelope,
            "caix": 0.18 * glycolysis_rate * state["caix"] * envelope,
            "bicarbonate_buffer": 0.12 * (state["nbce1"] + state["nbcn1"]) / 2.0 * envelope,
        }
