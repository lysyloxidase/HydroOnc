import numpy as np
import pytest

from hydroonc_therapy.bragg_peak import BraggPeakSimulator
from hydroonc_therapy.dose_prediction import DosePredictionUNet
from hydroonc_therapy.dvh import compute_dvh
from hydroonc_therapy.flash import FLASHProtonTherapy
from hydroonc_therapy.rbe_model import VariableRBE
from hydroonc_therapy.rl_planner import PPOProtonPlanner


def test_bragg_peak_for_150_mev_proton_in_water_near_clinical_range():
    simulator = BraggPeakSimulator()
    curve = simulator.bragg_curve_full(150.0, depth_range_cm=(0.0, 30.0), n_points=1500)
    assert curve.range_cm == pytest.approx(15.6, abs=0.05)
    assert curve.peak_depth_cm == pytest.approx(15.6, abs=0.12)
    assert np.max(curve.dose) == pytest.approx(1.0)


def test_sobp_flat_region_covers_target_volume_within_five_percent():
    simulator = BraggPeakSimulator()
    sobp = simulator.sobp_full(target_depth_range_cm=(5.0, 10.0), n_energies=20)
    mask = (sobp.depth_cm >= 5.0) & (sobp.depth_cm <= 10.0)
    target = sobp.dose[mask]
    assert sobp.target_uniformity() < 0.05
    assert np.max(np.abs(target - np.mean(target))) / np.mean(target) < 0.05
    assert len(sobp.energies_MeV) == 20


def test_bethe_bloch_matches_bundled_pstar_reference_within_two_percent():
    simulator = BraggPeakSimulator()
    for energy in [70.0, 100.0, 150.0, 200.0, 250.0]:
        stopping = simulator.bethe_bloch(energy)
        reference = simulator.nist_pstar_reference(energy)
        assert stopping == pytest.approx(reference, rel=0.02)


def test_variable_rbe_models_increase_at_distal_edge():
    rbe = VariableRBE()
    values = rbe.all_models(LET_keV_um=9.0, alpha_beta_Gy=2.0, dose_Gy=2.0)
    assert set(values) == {"mcnamara", "wedenberg", "mcmahon", "carabe_fernandez"}
    assert all(value > 1.1 for value in values.values())


def test_3d_unet_compiles_and_predicts_dummy_64_cube():
    model = DosePredictionUNet(levels=5, base_filters=32)
    architecture = model.compile(input_shape=(64, 64, 64, 3))
    dummy = np.zeros((64, 64, 64, 3), dtype=float)
    dummy[..., 1] = 1.0
    dummy[16:48, 16:48, 20:52, 2] = 1.0
    dose = model.predict(dummy)
    assert architecture.compiled
    assert architecture.output_shape == (64, 64, 64)
    assert architecture.parameter_count > 0
    assert dose.shape == (64, 64, 64)
    assert np.max(dose) > 0.0


def test_flash_threshold_and_oxygen_response():
    flash = FLASHProtonTherapy()
    response = flash.response(dose_Gy=8.0, delivery_time_s=0.1, oxygen_initial_mmHg=30.0)
    assert response.dose_rate_Gy_s == pytest.approx(80.0)
    assert response.is_flash
    assert flash.is_flash(8.0, 0.1)
    assert response.oxygen_after_mmHg < response.oxygen_before_mmHg
    assert response.normal_tissue_sparing_factor < 1.0


def test_dvh_computes_cumulative_dose_volume_histogram():
    x = np.linspace(0.0, 70.0, 16)
    dose = np.broadcast_to(x[None, None, :], (16, 16, 16))
    mask = np.ones_like(dose, dtype=bool)
    dvh = compute_dvh(dose, mask, n_bins=32, structure_name="PTV")
    assert dvh.structure_name == "PTV"
    assert dvh.volume_fraction[0] == pytest.approx(1.0)
    assert dvh.volume_fraction[-1] > 0.0
    assert dvh.v_at_dose(35.0) == pytest.approx(0.5, abs=0.08)
    assert dvh.d_at_volume(0.5) == pytest.approx(35.0, abs=5.0)


def test_ppo_planner_returns_bounded_plan_quality_score():
    planner = PPOProtonPlanner()
    state = np.linspace(0.0, 1.0, 64)
    result = planner.propose(state)
    assert abs(sum(result.objective_weights.values()) - 1.0) < 1.0e-12
    assert 0.0 <= result.reward <= 150.0
    assert result.target_coverage > 0.90
