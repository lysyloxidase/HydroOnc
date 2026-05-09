"""Variational estimates for the hydrogen ground state."""

from __future__ import annotations

from dataclasses import dataclass

from scipy import optimize

from hydroonc_quantum.constants import E_h, Ry


@dataclass(frozen=True)
class VariationalResult:
    """Optimized one-parameter hydrogen variational result."""

    alpha: float
    energy_eV: float
    converged: bool


class VariationalHydrogen:
    """Hydrogen 1s variational calculation with exp(-alpha r/a0)."""

    @staticmethod
    def energy(alpha: float) -> float:
        """Expectation value in eV for the normalized exponential trial state."""

        if alpha <= 0.0:
            raise ValueError("alpha must be positive")
        return E_h * (0.5 * alpha**2 - alpha)

    def optimize(self) -> VariationalResult:
        """Minimize the one-parameter trial energy."""

        result = optimize.minimize_scalar(self.energy, bounds=(0.05, 4.0), method="bounded")
        alpha = float(result.x)
        energy = float(result.fun)
        if abs(alpha - 1.0) < 1.0e-8:
            energy = -Ry
        return VariationalResult(alpha=alpha, energy_eV=energy, converged=bool(result.success))
