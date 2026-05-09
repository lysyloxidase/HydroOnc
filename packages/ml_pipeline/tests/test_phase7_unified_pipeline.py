import json
from pathlib import Path

import pytest
from hydroonc_h2.immune_model import H2ImmuneModel
from hydroonc_hbond.p53_mutations import P53MutantAnalyzer
from hydroonc_ph.reaction_diffusion import TumorpHReactionDiffusion
from hydroonc_quantum.constants import a_0, alpha, c, e, hbar, m_e
from hydroonc_quantum.visualization import classify_orbital_shape
from hydroonc_therapy.bragg_peak import BraggPeakSimulator

from hydroonc_pipeline.citations import CITATIONS
from hydroonc_pipeline.cli import main
from hydroonc_pipeline.unified import HydroOncPipeline
from hydroonc_pipeline.visualization import validate_visualization_manifest


def test_unified_pipeline_use_case_1_runs_quantum_to_clinical():
    report = HydroOncPipeline().run_use_case(
        "use-case-1",
        {
            "tumor_type": "NSCLC",
            "stage": "IV",
            "tumor_ph_e": 6.7,
            "mutation": "R175H",
            "pd1_cd8_percent": 42.0,
        },
    )
    modules = {stage.module for stage in report.stages}
    assert report.completed
    assert "quantum" in modules
    assert "qmmm" in modules
    assert "constant_ph_md" in modules
    assert "continuum_pH_PINN" in modules
    assert "clinical_response" in modules
    assert report.as_dict()["completed"]


def test_visualization_manifest_has_all_six_pages_and_caveats():
    manifest = validate_visualization_manifest()
    assert manifest["valid"]
    assert manifest["page_count"] == 6
    assert set(manifest["slugs"]) == {"quantum", "proton-transport", "tumor-ph", "proton-therapy", "hbond", "h2-therapy"}

    static_html = Path("apps/web/static/index.html").read_text(encoding="utf-8")
    for phrase in ["Quantum Explorer", "Proton Therapy", "Oncoprotein H-Bonds", "H₂ Therapy"]:
        assert phrase in static_html
    for caveat in ["Variable RBE is not ICRU standard", "Research platform, not clinical decision support"]:
        assert caveat in static_html


def test_quantum_explorer_orbital_shapes_render_correctly():
    assert classify_orbital_shape(1, 0, 0) == "spherical"
    assert classify_orbital_shape(2, 1, 0) == "dumbbell"
    assert classify_orbital_shape(3, 2, 0) == "clover"


def test_bragg_interactive_physics_gate_matches_range_reference():
    simulator = BraggPeakSimulator()
    curve = simulator.bragg_curve_full(150.0, depth_range_cm=(0.0, 30.0), n_points=1200)
    assert curve.peak_depth_cm == pytest.approx(15.6, rel=0.02)
    assert simulator.bethe_bloch(150.0) == pytest.approx(simulator.nist_pstar_reference(150.0), rel=0.02)


def test_tumor_ph_animation_reaction_diffusion_evolves_correctly():
    solution = TumorpHReactionDiffusion().solve_fenics(mesh={"n_points": 72}, params={"core_ph": 6.5, "boundary_ph": 7.3})
    assert solution.converged
    assert solution.core_ph == pytest.approx(6.5)
    assert solution.boundary_ph == pytest.approx(7.3)
    assert all(delta >= -1.0e-12 for delta in (solution.ph[1:] - solution.ph[:-1]))


def test_protein_viewer_difference_map_for_p53_r175h():
    comparison = P53MutantAnalyzer().analyze_mutation("R175H")
    assert comparison.hbond_loss_count >= 4
    assert comparison.mutant_count < comparison.wild_type_count
    assert "L2/L3" in P53MutantAnalyzer().identify_druggable_cavities("R175H")[0]["site"]


def test_h2_immune_combination_outperforms_monotherapy():
    synergy = H2ImmuneModel().synergy_score()
    assert synergy["combo_final_tumor"] < synergy["h2_final_tumor"]
    assert synergy["combo_final_tumor"] < synergy["anti_pd1_final_tumor"]


def test_all_caveats_are_prominent_in_docs_and_manifest():
    readme = Path("README.md").read_text(encoding="utf-8")
    h2_doc = Path("docs/H2_THERAPY.md").read_text(encoding="utf-8")
    therapy_doc = Path("docs/PROTON_THERAPY.md").read_text(encoding="utf-8")
    hbond_doc = Path("docs/HBOND_DESIGN.md").read_text(encoding="utf-8")
    joined = "\n".join([readme, h2_doc, therapy_doc, hbond_doc, json.dumps(validate_visualization_manifest())])
    for phrase in [
        "H2 selective scavenging is CONTESTED",
        "APR-246 Phase 3 FAILED",
        "Variable RBE",
        "FLASH FAST-01",
        "H2 cancer trials are small",
        "RESEARCH platform",
    ]:
        assert phrase.lower() in joined.lower()


def test_physical_constants_remain_codata_2018_traceable():
    assert m_e == pytest.approx(9.1093837015e-31)
    assert e == pytest.approx(1.602176634e-19)
    assert hbar == pytest.approx(1.054571817e-34)
    assert c == pytest.approx(299792458.0)
    assert a_0 == pytest.approx(5.29177210903e-11, rel=2.0e-9)
    assert alpha == pytest.approx(1 / 137.035999, rel=2.0e-7)


def test_biological_parameters_have_citation_registry_coverage():
    modules = {citation.module for citation in CITATIONS}
    assert {"quantum", "proton_transport", "tumor_ph", "therapy", "hbond", "h2"}.issubset(modules)
    assert len(CITATIONS) >= 12


def test_cli_all_documented_commands_are_functional(capsys):
    commands = [
        ["quantum", "orbital", "--n", "3", "--l", "2", "--m", "0", "--visualize"],
        ["quantum", "spectrum", "--series", "balmer"],
        ["quantum", "h2", "--method", "ccsd", "--basis", "cc-pVTZ", "--scan"],
        ["quantum", "energy", "--n", "2"],
        ["proton", "grotthuss", "--n-waters", "20", "--temperature", "300"],
        ["proton", "cph", "--residue", "HIS", "--ph", "6.7"],
        ["proton", "hbond", "--base-pair", "G-C"],
        ["proton", "qmmm", "--radius", "3.0"],
        ["ph", "warburg", "--glycolysis-rate", "1.5"],
        ["ph", "warburg", "--inhibit", "CAIX", "--dose", "10uM"],
        ["ph", "pinn", "--mesh", "tumor_mri.vtk", "--train", "--epochs", "5000"],
        ["ph", "cest"],
        ["ph", "inhibitor", "--target", "NHE1", "--dose", "10uM"],
        ["therapy", "bragg", "--energy", "150", "--material", "water"],
        ["therapy", "sobp", "--target-depth", "5-10", "--optimize"],
        ["therapy", "dose-predict", "--ct", "patient_ct.nii", "--model", "unet3d"],
        ["therapy", "rl-plan", "--ct", "patient_ct.nii", "--algorithm", "ppo"],
        ["therapy", "flash", "--dose", "8", "--time", "0.1"],
        ["therapy", "rbe", "--let", "9"],
        ["therapy", "dvh", "--dose", "70"],
        ["hbond", "p53", "--mutation", "R175H", "--compare-wt", "--md-ns", "100"],
        ["hbond", "kras", "--mutation", "G12C", "--drug", "sotorasib", "--hbond-analysis"],
        ["hbond", "egfr", "--mutation", "T790M", "--drug", "gefitinib"],
        ["hbond", "gnn", "--protein", "target.pdb", "--ligand", "candidate.sdf"],
        ["h2", "ros", "--h2-conc", "0.78", "--oh-conc", "1.0", "--simulate"],
        ["h2", "immune", "--h2-dose", "3h-day", "--anti-pd1", "--simulate"],
        ["h2", "predict", "--pd1-cd8", "42", "--tumor-ph", "6.7"],
        ["pipeline", "use-case-1"],
        ["pipeline", "use-case-2"],
        ["pipeline", "use-case-3"],
        ["meta", "caveats"],
        ["meta", "citations"],
    ]
    assert len(commands) >= 20
    for command in commands:
        assert main(command) == 0
        payload = json.loads(capsys.readouterr().out)
        assert "command" in payload or "case" in payload
