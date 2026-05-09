"""Hydrogen orbital helpers for grids and real-valued combinations."""

from __future__ import annotations

from typing import Optional, Tuple

import numpy as np

from hydroonc_quantum.constants import a_0
from hydroonc_quantum.hydrogen_atom import HydrogenAtom


def cartesian_to_spherical(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Convert Cartesian coordinates to r, theta, phi arrays."""

    r = np.sqrt(x * x + y * y + z * z)
    theta = np.zeros_like(r)
    np.divide(z, r, out=theta, where=r > 0.0)
    theta = np.arccos(np.clip(theta, -1.0, 1.0))
    phi = np.mod(np.arctan2(y, x), 2.0 * np.pi)
    return r, theta, phi


def real_orbital(
    n: int,
    l: int,
    m: int,
    r: np.ndarray,
    theta: np.ndarray,
    phi: np.ndarray,
    atom: Optional[HydrogenAtom] = None,
) -> np.ndarray:
    """Return a real-valued hydrogen orbital.

    For ``m = 0`` this is the analytical orbital. For nonzero ``m`` the usual
    real tesseral combinations are used.
    """

    atom = atom or HydrogenAtom()
    if m == 0:
        return np.real(atom.wavefunction(n, l, 0, r, theta, phi))
    abs_m = abs(m)
    psi_pos = atom.wavefunction(n, l, abs_m, r, theta, phi)
    psi_neg = atom.wavefunction(n, l, -abs_m, r, theta, phi)
    if m > 0:
        return np.real((psi_neg + (-1) ** abs_m * psi_pos) / np.sqrt(2.0))
    return np.real((psi_neg - (-1) ** abs_m * psi_pos) / (1j * np.sqrt(2.0)))


def orbital_grid(
    n: int,
    l: int,
    m: int,
    *,
    grid_points: int = 64,
    extent_bohr: float = 10.0,
    real: bool = True,
    atom: Optional[HydrogenAtom] = None,
) -> dict:
    """Sample an orbital on a cubic Cartesian grid.

    Coordinates are returned in meters; ``extent_bohr`` is measured in Bohr
    radii on each side of the origin.
    """

    if grid_points < 8:
        raise ValueError("grid_points must be at least 8")
    atom = atom or HydrogenAtom()
    axis = np.linspace(-extent_bohr * a_0, extent_bohr * a_0, grid_points)
    x, y, z = np.meshgrid(axis, axis, axis, indexing="ij")
    r, theta, phi = cartesian_to_spherical(x, y, z)
    if real:
        psi = real_orbital(n, l, m, r, theta, phi, atom=atom)
    else:
        psi = atom.wavefunction(n, l, m, r, theta, phi)
    return {
        "axis_m": axis,
        "x_m": x,
        "y_m": y,
        "z_m": z,
        "r_m": r,
        "theta": theta,
        "phi": phi,
        "psi": psi,
    }


def probability_density_grid(
    n: int,
    l: int,
    m: int,
    *,
    grid_points: int = 64,
    extent_bohr: float = 10.0,
    real: bool = True,
    atom: Optional[HydrogenAtom] = None,
) -> dict:
    """Sample |psi|^2 on a cubic grid."""

    grid = orbital_grid(
        n,
        l,
        m,
        grid_points=grid_points,
        extent_bohr=extent_bohr,
        real=real,
        atom=atom,
    )
    grid["density"] = np.abs(grid["psi"]) ** 2
    return grid
