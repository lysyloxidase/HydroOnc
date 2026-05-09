import numpy as np
import pytest

from hydroonc_proton.constant_ph_md import ConstantpHMD
from hydroonc_proton.constants import D_H_PLUS_WATER_CM2_S
from hydroonc_proton.grotthuss import GrotthussSimulator
from hydroonc_proton.hbond_analysis import HBondAnalyzer
from hydroonc_proton.pka import PROPKA3Runner
from hydroonc_proton.qmmm import AIMDConfig, QMMMProtonTransport


def test_grotthuss_hopping_rate_matches_literature_scale():
    simulator = GrotthussSimulator(n_waters=12)
    trajectory = simulator.simulate(duration_ps=100.0)
    assert trajectory.hopping_rate_per_ps == pytest.approx(1.0, abs=0.02)
    assert trajectory.hop_events[0].carrier_state == "Zundel"
    assert 0.85 <= simulator.structural_diffusion_share() <= 0.95


def test_proton_diffusion_reproduces_experimental_value_within_20_percent():
    simulator = GrotthussSimulator()
    assert simulator.diffusion_coefficient_cm2_s() == pytest.approx(
        D_H_PLUS_WATER_CM2_S,
        rel=0.20,
    )
    assert simulator.diffusion_error_fraction() < 0.20


def test_ms_evb_probabilities_are_normalized_and_localized():
    simulator = GrotthussSimulator(n_waters=9)
    probabilities = simulator.valence_state_probabilities(proton_site=4)
    assert probabilities.sum() == pytest.approx(1.0)
    assert int(np.argmax(probabilities)) == 4


def test_constant_ph_histidine_population_shifts_between_tumor_and_normal_ph():
    cph = ConstantpHMD(seed=1)
    acidic = cph.protonation_state("HIS", ph=6.0)
    normal = cph.protonation_state("HIS", ph=7.4)
    assert acidic.protonated_fraction > 0.70
    assert normal.protonated_fraction < 0.15
    assert acidic.dominant_state == "cationic"
    assert normal.dominant_state == "neutral"


def test_openmm_simulation_surface_runs_one_ns_without_crash():
    simulation = ConstantpHMD().setup_tumor_simulation(
        pdb_path="data/structures/toy_missing.pdb",
        ph=6.7,
    )
    report = simulation.run(duration_ns=1.0)
    assert report.completed
    assert report.simulated_time_ns == pytest.approx(1.0)
    assert report.steps > 0
    assert report.backend in {"mock", "openmm"}


def test_hbond_analyzer_identifies_watson_crick_counts():
    analyzer = HBondAnalyzer()
    at_bonds = analyzer.analyze_base_pair("A-T")
    gc_bonds = analyzer.analyze_base_pair("G-C")
    assert len(at_bonds) == 2
    assert len(gc_bonds) == 3
    assert all(2.5 <= bond.distance_DA_A <= 3.5 for bond in at_bonds + gc_bonds)
    assert all(bond.angle_DHA_deg > 150.0 for bond in at_bonds + gc_bonds)


def test_qmmm_cp2k_aimd_defaults_capture_recommended_validation_setup():
    qmmm = QMMMProtonTransport()
    oxygens = np.array([[0.0, 0.0, 0.0], [2.7, 0.0, 0.0], [5.4, 0.0, 0.0]])
    region = qmmm.select_water_wire_region(oxygens, np.array([2.6, 0.1, 0.0]), radius_A=3.0)
    config = qmmm.build_config(region)
    aimd = AIMDConfig()
    cp2k = qmmm.cp2k_aimd_input("hydroonc_water_wire", "O 0 0 0\nH 0 0 1", aimd)
    assert region.qm_atom_indices == (0, 1, 2)
    assert config.embedding == "electrostatic"
    assert QMMMProtonTransport.validate_recommended_aimd(aimd)
    assert "RUN_TYPE MD" in cp2k
    assert "TIMESTEP 0.5" in cp2k
    assert "&BLYP" in cp2k
    assert "TYPE D3" in cp2k


def test_propka3_lysozyme_benchmarks_match_within_one_pka_unit():
    runner = PROPKA3Runner()
    benchmarks = runner.lysozyme_benchmarks()
    assert benchmarks
    assert runner.lysozyme_within_tolerance(1.0)
    assert max(item.absolute_error for item in benchmarks) <= 1.0
