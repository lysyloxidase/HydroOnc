"""Perturbative hydrogen energy corrections."""

from __future__ import annotations

from hydroonc_quantum.constants import Ry, a_0, alpha, e, epsilon_0, mu_B, pi


def fine_structure_correction(n: int, j: float) -> float:
    """Return the leading fine-structure correction in eV.

    This uses the standard order-alpha^2 correction to the Bohr level.
    """

    if n < 1:
        raise ValueError("n must be positive")
    if j <= 0.0:
        raise ValueError("j must be positive")
    return -Ry / n**2 * (alpha**2 / n**2) * (n / (j + 0.5) - 0.75)


def zeeman_shift(B_tesla: float, m_l: int, m_s: float = 0.5, g_s: float = 2.00231930436256) -> float:
    """Linear Zeeman shift in eV for weak magnetic fields."""

    return (mu_B * B_tesla * (m_l + g_s * m_s)) / e


def stark_shift_ground_state(field_V_per_m: float) -> float:
    """Quadratic Stark shift for hydrogen 1s in eV."""

    polarizability_SI = 4.0 * pi * epsilon_0 * 4.5 * a_0**3
    return -0.5 * polarizability_SI * field_V_per_m**2 / e
