"""3D orbital visualization primitives."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from hydroonc_quantum.constants import a_0
from hydroonc_quantum.orbitals import probability_density_grid


@dataclass(frozen=True)
class OrbitalIsosurface:
    """Mesh or point-cloud representation of an orbital isosurface."""

    vertices_m: np.ndarray
    faces: np.ndarray
    values: np.ndarray
    isovalue: float
    n: int
    l: int
    m: int
    shape: str


def orbital_isosurface(
    n: int,
    l: int,
    m: int,
    *,
    grid_points: int = 48,
    extent_bohr: float = 10.0,
    isovalue_fraction: float = 0.1,
    real: bool = True,
) -> OrbitalIsosurface:
    """Build an orbital density isosurface.

    If scikit-image is installed, marching cubes is used. Otherwise the return
    value is a thresholded point cloud with an empty face array.
    """

    if not 0.0 < isovalue_fraction < 1.0:
        raise ValueError("isovalue_fraction must be between 0 and 1")
    grid = probability_density_grid(
        n,
        l,
        m,
        grid_points=grid_points,
        extent_bohr=extent_bohr,
        real=real,
    )
    density = grid["density"]
    isovalue = float(np.max(density) * isovalue_fraction)
    axis = grid["axis_m"]
    spacing = (axis[1] - axis[0], axis[1] - axis[0], axis[1] - axis[0])
    origin = np.array([axis[0], axis[0], axis[0]])

    try:
        from skimage import measure

        vertices, faces, normals, values = measure.marching_cubes(density, level=isovalue, spacing=spacing)
        del normals
        vertices = vertices + origin
    except Exception:
        mask = density >= isovalue
        vertices = np.column_stack(
            [grid["x_m"][mask], grid["y_m"][mask], grid["z_m"][mask]]
        )
        faces = np.empty((0, 3), dtype=int)
        values = density[mask]

    return OrbitalIsosurface(
        vertices_m=np.asarray(vertices, dtype=float),
        faces=np.asarray(faces, dtype=int),
        values=np.asarray(values, dtype=float),
        isovalue=isovalue,
        n=n,
        l=l,
        m=m,
        shape=classify_orbital_shape(n, l, m),
    )


def classify_orbital_shape(n: int, l: int, m: int) -> str:
    """Return the expected qualitative orbital shape."""

    del n, m
    if l == 0:
        return "spherical"
    if l == 1:
        return "dumbbell"
    if l == 2:
        return "clover"
    return "multi-lobed"


def radial_extent_bohr(surface: OrbitalIsosurface) -> tuple[float, float, float]:
    """Return x/y/z half-extents of an isosurface in Bohr radii."""

    if surface.vertices_m.size == 0:
        return (0.0, 0.0, 0.0)
    extents_m = np.max(np.abs(surface.vertices_m), axis=0)
    extents_bohr = extents_m / a_0
    return tuple(float(value) for value in extents_bohr)
