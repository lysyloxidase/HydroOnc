"""Grotthuss proton hopping along water wires."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from hydroonc_proton.constants import (
    D_H_PLUS_WATER_CM2_S,
    DEFAULT_HOP_CORRELATION_FACTOR,
    DEFAULT_HOP_RATE_PER_PS,
    DEFAULT_STRUCTURAL_DIFFUSION_FRACTION,
    DEFAULT_WATER_WIRE_SPACING_A,
)


@dataclass(frozen=True)
class HopEvent:
    """One proton-transfer event on a water wire."""

    time_ps: float
    from_site: int
    to_site: int
    carrier_state: str
    barrier_kj_mol: float


@dataclass(frozen=True)
class WaterWireTrajectory:
    """Trajectory summary for proton hopping on a one-dimensional water wire."""

    n_sites: int
    duration_ps: float
    hop_events: tuple[HopEvent, ...]
    site_history: np.ndarray
    time_ps: np.ndarray

    @property
    def n_hops(self) -> int:
        return len(self.hop_events)

    @property
    def hopping_rate_per_ps(self) -> float:
        if self.duration_ps <= 0.0:
            return 0.0
        return self.n_hops / self.duration_ps

    @property
    def net_displacement_sites(self) -> int:
        if self.site_history.size == 0:
            return 0
        return int(self.site_history[-1] - self.site_history[0])


class GrotthussSimulator:
    """Simulate proton hopping along water wires.

    The model is a simplified MS-EVB-like representation: each water oxygen is a
    possible excess-proton valence state, nearest-neighbor states are coupled,
    and hopping events interconvert Eigen-like and Zundel-like configurations.
    """

    def __init__(
        self,
        n_waters: int = 32,
        hop_rate_per_ps: float = DEFAULT_HOP_RATE_PER_PS,
        oxygen_spacing_A: float = DEFAULT_WATER_WIRE_SPACING_A,
        structural_fraction: float = DEFAULT_STRUCTURAL_DIFFUSION_FRACTION,
        correlation_factor: float = DEFAULT_HOP_CORRELATION_FACTOR,
        barrier_kj_mol: float = 10.0,
        seed: Optional[int] = None,
    ) -> None:
        if n_waters < 2:
            raise ValueError("n_waters must be at least 2")
        if hop_rate_per_ps <= 0.0:
            raise ValueError("hop_rate_per_ps must be positive")
        if oxygen_spacing_A <= 0.0:
            raise ValueError("oxygen_spacing_A must be positive")
        if not 0.0 < structural_fraction <= 1.0:
            raise ValueError("structural_fraction must be in (0, 1]")
        if not 0.0 < correlation_factor <= 1.0:
            raise ValueError("correlation_factor must be in (0, 1]")

        self.n_waters = n_waters
        self.hop_rate_per_ps = hop_rate_per_ps
        self.oxygen_spacing_A = oxygen_spacing_A
        self.structural_fraction = structural_fraction
        self.correlation_factor = correlation_factor
        self.barrier_kj_mol = barrier_kj_mol
        self.rng = np.random.default_rng(seed)

    def ms_evb_hamiltonian(
        self,
        proton_site: int,
        electric_field_strength: float = 0.0,
    ) -> np.ndarray:
        """Return a nearest-neighbor empirical valence-bond Hamiltonian.

        Energies are in kJ/mol. The optional field gives a linear bias along the
        wire and is useful for membrane-potential sensitivity studies.
        """

        self._validate_site(proton_site)
        indices = np.arange(self.n_waters)
        diagonal = 8.0 * (indices - proton_site) ** 2
        diagonal += electric_field_strength * (indices - proton_site) * self.oxygen_spacing_A
        hamiltonian = np.diag(diagonal.astype(float))
        coupling = -18.0 * np.exp(-abs(self.oxygen_spacing_A - 2.5) / 0.35)
        for idx in range(self.n_waters - 1):
            hamiltonian[idx, idx + 1] = coupling
            hamiltonian[idx + 1, idx] = coupling
        return hamiltonian

    def valence_state_probabilities(self, proton_site: int) -> np.ndarray:
        """Ground-state excess-proton probabilities from the EVB Hamiltonian."""

        eigenvalues, eigenvectors = np.linalg.eigh(self.ms_evb_hamiltonian(proton_site))
        del eigenvalues
        ground = eigenvectors[:, 0]
        probabilities = np.abs(ground) ** 2
        return probabilities / probabilities.sum()

    def simulate(
        self,
        duration_ps: float = 100.0,
        start_site: Optional[int] = None,
        stochastic: bool = False,
    ) -> WaterWireTrajectory:
        """Simulate hopping events for ``duration_ps``.

        The deterministic default schedules one hop every ``1/rate`` ps, which
        makes regression tests stable. Set ``stochastic=True`` for a Poisson
        continuous-time random walk.
        """

        if duration_ps <= 0.0:
            raise ValueError("duration_ps must be positive")
        site = self.n_waters // 2 if start_site is None else start_site
        self._validate_site(site)
        event_times = self._event_times(duration_ps, stochastic)
        site_history = [site]
        time_history = [0.0]
        events = []
        for index, time_ps in enumerate(event_times):
            direction = int(self.rng.choice([-1, 1])) if stochastic else (1 if index % 2 == 0 else -1)
            to_site = self._reflect_site(site + direction)
            if to_site == site:
                direction *= -1
                to_site = self._reflect_site(site + direction)
            state = "Zundel" if index % 2 == 0 else "Eigen"
            events.append(
                HopEvent(
                    time_ps=float(time_ps),
                    from_site=int(site),
                    to_site=int(to_site),
                    carrier_state=state,
                    barrier_kj_mol=self.barrier_kj_mol,
                )
            )
            site = to_site
            time_history.append(float(time_ps))
            site_history.append(site)

        if time_history[-1] < duration_ps:
            time_history.append(float(duration_ps))
            site_history.append(site)
        return WaterWireTrajectory(
            n_sites=self.n_waters,
            duration_ps=float(duration_ps),
            hop_events=tuple(events),
            site_history=np.asarray(site_history, dtype=int),
            time_ps=np.asarray(time_history, dtype=float),
        )

    def diffusion_coefficient_A2_ps(self) -> float:
        """Return calibrated bulk proton diffusion in Angstrom^2/ps."""

        return (
            self.structural_fraction
            * self.correlation_factor
            * self.hop_rate_per_ps
            * self.oxygen_spacing_A**2
            / 6.0
        )

    def diffusion_coefficient_cm2_s(self) -> float:
        """Return calibrated bulk proton diffusion in cm^2/s."""

        return self.diffusion_coefficient_A2_ps() * 1.0e-4

    def diffusion_error_fraction(self, reference: float = D_H_PLUS_WATER_CM2_S) -> float:
        """Relative error against an experimental diffusion reference."""

        return abs(self.diffusion_coefficient_cm2_s() - reference) / reference

    def structural_diffusion_share(self) -> float:
        """Fraction of net proton transport attributed to hopping."""

        return self.structural_fraction

    def _event_times(self, duration_ps: float, stochastic: bool) -> np.ndarray:
        if not stochastic:
            interval = 1.0 / self.hop_rate_per_ps
            count = int(np.floor(duration_ps / interval))
            if count == 0:
                return np.empty(0, dtype=float)
            return interval * np.arange(1, count + 1, dtype=float)

        times = []
        time_ps = 0.0
        while True:
            time_ps += float(self.rng.exponential(1.0 / self.hop_rate_per_ps))
            if time_ps > duration_ps:
                break
            times.append(time_ps)
        return np.asarray(times, dtype=float)

    def _reflect_site(self, site: int) -> int:
        if site < 0:
            return 1
        if site >= self.n_waters:
            return self.n_waters - 2
        return site

    def _validate_site(self, site: int) -> None:
        if site < 0 or site >= self.n_waters:
            raise ValueError("proton site is outside the water wire")
