"""Quantum-mechanical foundation for HydroOnc Phase 1."""

from hydroonc_quantum.constants import (
    E_h,
    E_h_J,
    R_H,
    R_inf,
    Ry,
    Ry_H,
    a_0,
    a_H,
    alpha,
    c,
    e,
    epsilon_0,
    h,
    hbar,
    k_B,
    m_e,
    m_p,
    mu_B,
    mu_H,
)
from hydroonc_quantum.h2_molecule import H2Molecule, PotentialEnergyCurve
from hydroonc_quantum.hydrogen_atom import HydrogenAtom
from hydroonc_quantum.orbitals import (
    cartesian_to_spherical,
    orbital_grid,
    probability_density_grid,
    real_orbital,
)
from hydroonc_quantum.perturbation import (
    fine_structure_correction,
    stark_shift_ground_state,
    zeeman_shift,
)
from hydroonc_quantum.spectrum import HydrogenSpectrum, SPECTRAL_SERIES
from hydroonc_quantum.tdse import SplitOperatorTDSE
from hydroonc_quantum.variational import VariationalHydrogen
from hydroonc_quantum.visualization import (
    OrbitalIsosurface,
    classify_orbital_shape,
    orbital_isosurface,
)

__all__ = [
    "E_h",
    "E_h_J",
    "R_H",
    "R_inf",
    "Ry",
    "Ry_H",
    "a_0",
    "a_H",
    "alpha",
    "c",
    "e",
    "epsilon_0",
    "h",
    "hbar",
    "k_B",
    "m_e",
    "m_p",
    "mu_B",
    "mu_H",
    "H2Molecule",
    "HydrogenAtom",
    "HydrogenSpectrum",
    "OrbitalIsosurface",
    "PotentialEnergyCurve",
    "SPECTRAL_SERIES",
    "SplitOperatorTDSE",
    "VariationalHydrogen",
    "cartesian_to_spherical",
    "classify_orbital_shape",
    "fine_structure_correction",
    "orbital_grid",
    "orbital_isosurface",
    "probability_density_grid",
    "real_orbital",
    "stark_shift_ground_state",
    "zeeman_shift",
]
