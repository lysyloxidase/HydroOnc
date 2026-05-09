"""PPO-style reinforcement-learning helper for proton PBS planning."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class PPOPlanResult:
    """One RL planning policy output."""

    objective_weights: dict[str, float]
    reward: float
    target_coverage: float
    oar_penalty: float


class PPOProtonPlanner:
    """PPO-based automated proton PBS treatment planning surface."""

    def __init__(self, action_dim: int = 4, seed: int = 23) -> None:
        self.action_dim = action_dim
        self.rng = np.random.default_rng(seed)

    def action_to_weights(self, action: np.ndarray) -> dict[str, float]:
        """Map continuous PPO action to normalized objective weights."""

        values = np.asarray(action, dtype=float)
        if values.shape != (self.action_dim,):
            raise ValueError("action has wrong dimension")
        positive = np.exp(values - np.max(values))
        positive /= np.sum(positive)
        labels = ["target_coverage", "distal_oar", "proximal_oar", "robustness"]
        return {labels[idx]: float(positive[idx]) for idx in range(self.action_dim)}

    def reward(self, target_d95: float, oar_mean_dose: float, robustness_penalty: float = 0.0) -> float:
        """150-point plan-quality score."""

        coverage_score = 100.0 * np.clip(target_d95 / 0.95, 0.0, 1.2)
        oar_score = 50.0 * np.clip(1.0 - oar_mean_dose / 45.0, 0.0, 1.0)
        return float(np.clip(coverage_score + oar_score - robustness_penalty, 0.0, 150.0))

    def propose(self, state: np.ndarray) -> PPOPlanResult:
        """Produce deterministic objective weights from a planning state."""

        state = np.asarray(state, dtype=float)
        summary = np.array(
            [
                np.mean(state),
                np.std(state),
                np.percentile(state, 95),
                np.percentile(state, 5),
            ]
        )
        action = (summary - np.mean(summary)) / (np.std(summary) + 1.0e-9)
        weights = self.action_to_weights(action)
        target_coverage = float(np.clip(0.90 + 0.08 * weights["target_coverage"], 0.0, 1.0))
        oar_penalty = float(35.0 * weights["distal_oar"] + 20.0 * weights["proximal_oar"])
        reward = self.reward(target_coverage, oar_penalty, 5.0 * weights["robustness"])
        return PPOPlanResult(weights, reward, target_coverage, oar_penalty)
