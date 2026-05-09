import math

import numpy as np
import pytest

from hydroonc_quantum.constants import Ry, a_0
from hydroonc_quantum.h2_molecule import H2Molecule
from hydroonc_quantum.hydrogen_atom import HydrogenAtom
from hydroonc_quantum.spectrum import HydrogenSpectrum
from hydroonc_quantum.variational import VariationalHydrogen
from hydroonc_quantum.visualization import classify_orbital_shape, orbital_isosurface, radial_extent_bohr


def test_hydrogen_energy_levels_match_bohr_values():
    atom = HydrogenAtom()
    assert atom.energy(1) == pytest.approx(-13.6057, abs=5e-5)
    assert atom.energy(2) == pytest.approx(-3.4014, abs=5e-5)
    assert Ry == pytest.approx(13.6057, abs=5e-5)


def test_named_hydrogen_spectral_lines():
    atom = HydrogenAtom()
    spectrum = HydrogenSpectrum(atom)
    assert atom.transition_wavelength(3, 2) == pytest.approx(656.28, abs=0.01)
    assert atom.transition_wavelength(2, 1) == pytest.approx(121.57, abs=0.01)
    balmer = spectrum.series("Balmer", n_max=4)
    assert [round(line.wavelength_nm, 3) for line in balmer] == [656.281, 486.135]


def test_radial_probability_1s_peaks_at_bohr_radius():
    atom = HydrogenAtom()
    r = np.linspace(0.05 * a_0, 5.0 * a_0, 20_000)
    probability = atom.radial_probability(1, 0, r)
    r_peak = r[int(np.argmax(probability))]
    assert r_peak == pytest.approx(a_0, rel=5e-4)


def test_numerov_ground_state_energy_gate():
    atom = HydrogenAtom()
    radial, energy = atom.numerov_solve(1, 0, N_points=2000)
    assert radial.shape == (2000,)
    assert energy == pytest.approx(atom.energy(1), abs=1e-6)
    assert atom.last_numerov is not None


def test_variational_ground_state_recovers_exact_minimum():
    variational = VariationalHydrogen()
    result = variational.optimize()
    assert result.converged
    assert result.alpha == pytest.approx(1.0, abs=1e-4)
    assert result.energy_eV == pytest.approx(-13.6057, abs=5e-5)
    assert variational.energy(1.0) == pytest.approx(-13.6057, abs=5e-5)


def test_h2_ccsd_reference_or_backend_energy_and_pec_minimum():
    h2 = H2Molecule(reference_fallback=True)
    result = h2.compute_energy(method="ccsd", basis="cc-pVTZ")
    assert result["energy_Eh"] == pytest.approx(-1.1723, abs=0.01)

    pec = h2.potential_energy_curve(
        R_values=[0.50, 0.62, 0.7414, 0.90, 1.20],
        method="ccsd",
        basis="cc-pVTZ",
    )
    assert pec.equilibrium_bond_length == pytest.approx(0.7414, abs=0.02)


def test_orbital_visualization_shapes_are_distinguishable():
    assert classify_orbital_shape(1, 0, 0) == "spherical"
    assert classify_orbital_shape(2, 1, 0) == "dumbbell"

    one_s = orbital_isosurface(1, 0, 0, grid_points=24, extent_bohr=5.0, isovalue_fraction=0.15)
    two_p = orbital_isosurface(2, 1, 0, grid_points=24, extent_bohr=8.0, isovalue_fraction=0.15)
    assert one_s.vertices_m.shape[0] > 0
    assert two_p.vertices_m.shape[0] > 0

    x_extent, y_extent, z_extent = radial_extent_bohr(one_s)
    assert math.isclose(x_extent, y_extent, rel_tol=0.25)
    assert math.isclose(y_extent, z_extent, rel_tol=0.25)
    _, _, pz_extent = radial_extent_bohr(two_p)
    assert pz_extent > 2.0
