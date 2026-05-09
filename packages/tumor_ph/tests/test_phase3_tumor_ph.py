import numpy as np
import pytest

from hydroonc_ph.cest_ml import CESTpHMapper
from hydroonc_ph.reaction_diffusion import TumorpHReactionDiffusion
from hydroonc_ph.warburg import WarburgModel


def test_warburg_model_acidifies_extracellular_space_over_one_hour():
    model = WarburgModel()
    series = model.simulate_ph_dynamics(t_max_s=3600.0, dt_s=5.0)
    assert series.ph_e[0] == pytest.approx(7.4)
    assert series.final_ph_e == pytest.approx(6.70, abs=0.06)
    assert series.extracellular_drop > 0.6
    assert 7.0 <= series.final_ph_i <= 7.2


def test_nhe1_inhibition_lowers_intracellular_ph_measurably():
    model = WarburgModel()
    effect = model.inhibitor_effect("NHE1", concentration_uM=10.0)
    assert effect["target"] == "nhe1"
    assert -0.30 <= effect["delta_ph_i"] <= -0.20
    assert effect["inhibited_final_ph_i"] < effect["baseline_final_ph_i"]


def test_caix_inhibition_raises_extracellular_ph_toward_normal():
    model = WarburgModel()
    effect = model.inhibitor_effect("SLC-0111", concentration_uM=5.0)
    assert effect["target"] == "caix"
    assert effect["delta_ph_e"] > 0.12
    assert effect["inhibited_final_ph_e"] > effect["baseline_final_ph_e"]


def test_reaction_diffusion_builds_core_to_boundary_ph_gradient():
    solver = TumorpHReactionDiffusion()
    solution = solver.solve_fenics(
        mesh={"n_points": 80, "length_cm": 0.12},
        params={"core_ph": 6.5, "boundary_ph": 7.3},
    )
    assert solution.core_ph == pytest.approx(6.5)
    assert solution.boundary_ph == pytest.approx(7.3)
    assert np.all(np.diff(solution.ph) >= -1.0e-12)
    assert solution.ph[len(solution.ph) // 2] < 7.0


def test_fenics_fallback_solution_converges_on_test_geometry():
    solver = TumorpHReactionDiffusion()
    solution = solver.solve_fenics(mesh={"n_points": 64}, params={"core_ph": 6.5, "boundary_ph": 7.3})
    assert solution.method in {"fenics", "finite_difference_fenics_fallback"}
    assert solution.residual_history[0] > solution.residual_history[-1]
    assert solution.converged


def test_pinn_surrogate_matches_reference_solution_within_point_one_ph_unit():
    solver = TumorpHReactionDiffusion()
    pinn = solver.solve_pinn(domain_bounds=(0.0, 0.10), training_points=8000, epochs=2500)
    assert pinn.training_loss[0] > pinn.training_loss[-1]
    assert pinn.max_abs_error < 0.1
    predictions = pinn.predict(np.array([0.0, 0.10]))
    assert predictions[0] == pytest.approx(6.5, abs=0.05)
    assert predictions[-1] == pytest.approx(7.3, abs=0.05)


def test_cest_ml_architecture_compiles_and_predicts_dummy_z_spectrum():
    mapper = CESTpHMapper(seed=7)
    architecture = mapper.compile(input_features=64, hidden_units=24, residual_blocks=2)
    offsets = np.linspace(-5.0, 5.0, 64)
    dummy_z = 0.90 - 0.08 * np.exp(-0.5 * ((offsets - 3.5) / 0.8) ** 2)
    prediction = mapper.predict(dummy_z)
    assert architecture.compiled
    assert prediction.shape == (1,)
    assert 6.0 <= prediction[0] <= 7.5
    assert mapper.physics_penalty(prediction) == pytest.approx(0.0)
