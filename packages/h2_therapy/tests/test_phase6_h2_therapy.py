import numpy as np
import pytest

from hydroonc_h2.clinical_predictor import H2PatientFeatures, H2ResponsePredictor
from hydroonc_h2.constants import K_H2_OH_M_INV_S
from hydroonc_h2.h2_model import MolecularH2Model
from hydroonc_h2.immune_model import H2ImmuneModel


def test_h2_hydroxyl_rate_matches_literature_constant():
    model = MolecularH2Model()
    constants = model.selective_rate_constants()
    result = model.simulate_ros_scavenging()
    assert constants["H2+OH"] == pytest.approx(4.2e7)
    assert result.k_h2_oh_M_inv_s == pytest.approx(K_H2_OH_M_INV_S)
    assert result.k_h2_onoo_M_inv_s == pytest.approx(3.6e4)


def test_selective_scavenging_reduces_oh_not_superoxide_or_h2o2():
    result = MolecularH2Model().simulate_ros_scavenging(
        h2_concentration_mM=0.78,
        oh_radical_concentration_uM=1.0,
        superoxide_concentration_uM=10.0,
        h2o2_concentration_uM=50.0,
    )
    assert result.oh_remaining_uM < result.oh_initial_uM
    assert result.oh_reduction_fraction > 0.0
    assert result.superoxide_remaining_uM == pytest.approx(result.superoxide_initial_uM)
    assert result.h2o2_remaining_uM == pytest.approx(result.h2o2_initial_uM)
    assert result.superoxide_reduction_fraction == pytest.approx(0.0)
    assert result.h2o2_reduction_fraction == pytest.approx(0.0)


def test_competing_gsh_dominates_hydroxyl_quenching_over_h2():
    result = MolecularH2Model().simulate_ros_scavenging()
    assert result.gsh_oh_pseudo_first_order_s > 1000.0 * result.h2_oh_pseudo_first_order_s
    assert result.h2_fraction_of_oh_scavenging < 0.01
    assert "signaling" in result.caveat.lower()


def test_immune_ode_h2_rescues_exhaustion_and_shrinks_tumor():
    model = H2ImmuneModel()
    control = model.simulate(h2_concentration_mM=0.0)
    h2 = model.simulate(h2_concentration_mM=0.78)
    assert h2.final_t_eff > control.final_t_eff
    assert h2.final_t_exh < control.final_t_exh
    assert h2.final_tumor < control.final_tumor
    assert h2.final_tumor < h2.tumor[0]


def test_h2_dose_response_monotonically_improves_effector_t_cells():
    model = H2ImmuneModel()
    endpoints = model.dose_response(np.array([0.0, 0.2, 0.5, 0.78]))
    effector_values = [endpoint.final_t_eff for endpoint in endpoints]
    tumor_values = [endpoint.final_tumor for endpoint in endpoints]
    assert effector_values == sorted(effector_values)
    assert tumor_values == sorted(tumor_values, reverse=True)


def test_anti_pd1_h2_combination_beats_either_monotherapy():
    synergy = H2ImmuneModel().synergy_score(h2_concentration_mM=0.78)
    assert synergy["combo_final_tumor"] < synergy["h2_final_tumor"]
    assert synergy["combo_final_tumor"] < synergy["anti_pd1_final_tumor"]
    assert synergy["combo_benefit"] > synergy["h2_benefit"]
    assert synergy["combo_benefit"] > synergy["anti_pd1_benefit"]
    assert synergy["synergy_margin"] > 0.0


def test_research_grade_response_predictor_outputs_probability_and_caveat():
    features = H2PatientFeatures(
        pd1_cd8_percent=42.0,
        tumor_ph_e=6.65,
        nrf2_keap1_activity=0.8,
        coq10_index=0.7,
        complex_i_activity=0.65,
        tmb_mut_per_mb=12.0,
        tumor_type="NSCLC",
        stage="IV",
    )
    prediction = H2ResponsePredictor().predict(features)
    assert 0.0 <= prediction.probability_pfs_benefit <= 1.0
    assert prediction.risk_group in {"low", "intermediate", "high"}
    assert "RESEARCH-GRADE" in prediction.caveat
