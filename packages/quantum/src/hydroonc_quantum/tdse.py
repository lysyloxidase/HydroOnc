"""Split-operator FFT solver for the time-dependent Schrodinger equation."""

from __future__ import annotations

from typing import Callable, Union

import numpy as np
from scipy.integrate import trapezoid

from hydroonc_quantum.constants import hbar, m_e


class SplitOperatorTDSE:
    """One-dimensional split-operator TDSE solver.

    The implementation uses Strang splitting:
    exp(-i V dt / 2 hbar) F^-1 exp(-i p^2 dt / 2 m hbar) F exp(-i V dt / 2 hbar).
    """

    def __init__(
        self,
        x: np.ndarray,
        potential: Union[np.ndarray, Callable[[np.ndarray], np.ndarray]],
        *,
        mass: float = m_e,
        hbar_value: float = hbar,
    ) -> None:
        self.x = np.asarray(x, dtype=float)
        if self.x.ndim != 1 or self.x.size < 8:
            raise ValueError("x must be a one-dimensional grid with at least 8 points")
        dx = np.diff(self.x)
        if not np.allclose(dx, dx[0]):
            raise ValueError("x grid must be uniformly spaced")
        self.dx = float(dx[0])
        self.mass = mass
        self.hbar = hbar_value
        self.V = potential(self.x) if callable(potential) else np.asarray(potential, dtype=float)
        if self.V.shape != self.x.shape:
            raise ValueError("potential must have the same shape as x")
        self.k = 2.0 * np.pi * np.fft.fftfreq(self.x.size, d=self.dx)

    def step(self, psi: np.ndarray, dt: float, *, normalize: bool = True) -> np.ndarray:
        """Advance a wavefunction by one time step."""

        psi = np.asarray(psi, dtype=complex)
        if psi.shape != self.x.shape:
            raise ValueError("psi must have the same shape as x")
        half_potential = np.exp(-0.5j * self.V * dt / self.hbar)
        kinetic = np.exp(-0.5j * self.hbar * self.k**2 * dt / self.mass)
        updated = half_potential * psi
        updated = np.fft.ifft(kinetic * np.fft.fft(updated))
        updated = half_potential * updated
        if normalize:
            updated = self.normalize(updated)
        return updated

    def evolve(self, psi0: np.ndarray, dt: float, n_steps: int, *, normalize: bool = True) -> np.ndarray:
        """Return the final wavefunction after ``n_steps`` applications."""

        psi = np.asarray(psi0, dtype=complex)
        for _ in range(n_steps):
            psi = self.step(psi, dt, normalize=normalize)
        return psi

    def norm(self, psi: np.ndarray) -> float:
        """Return integral |psi|^2 dx."""

        return float(trapezoid(np.abs(psi) ** 2, self.x))

    def normalize(self, psi: np.ndarray) -> np.ndarray:
        """Normalize a wavefunction on the solver grid."""

        norm = self.norm(psi)
        if norm <= 0.0:
            raise ValueError("cannot normalize a zero wavefunction")
        return psi / np.sqrt(norm)
