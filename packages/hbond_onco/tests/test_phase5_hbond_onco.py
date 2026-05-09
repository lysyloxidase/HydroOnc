import numpy as np
import pytest

from hydroonc_hbond.fixtures import MANUAL_COUNTS
from hydroonc_hbond.geometry import automated_hbond_count
from hydroonc_hbond.gnn_drug_design import HBondGNN
from hydroonc_hbond.kras_egfr import KRASMutantAnalyzer
from hydroonc_hbond.p53_mutations import P53MutantAnalyzer


def test_p53_r175h_loses_at_least_four_hbonds_vs_wild_type_2xwr():
    analyzer = P53MutantAnalyzer()
    comparison = analyzer.analyze_mutation("R175H", ph=7.4)
    assert comparison.mutation_type == "structural"
    assert comparison.effect == "zinc_loss"
    assert comparison.hbond_loss_count >= 4
    assert comparison.wild_type_count - comparison.mutant_count >= 4
    assert "ZMC1" in comparison.drug_target


def test_p53_r273h_dna_contact_occupancy_drops_below_ten_percent():
    analyzer = P53MutantAnalyzer()
    comparison = analyzer.analyze_mutation("R273H")
    assert comparison.mutation_type == "contact"
    assert comparison.dna_contact_occupancy < 0.10
    assert "R273_DNA_phosphate" in comparison.lost_hbonds


def test_kras_g12c_his95_hbond_matches_sotorasib_cocrystal_geometry():
    analyzer = KRASMutantAnalyzer()
    result = analyzer.kras_hbond_geometry("G12C", "sotorasib")
    assert result.residue == "His95"
    assert result.distance_A == pytest.approx(2.96, abs=0.15)
    assert result.angle_deg > 120.0
    assert result.matches_cocrystal


def test_egfr_t790m_abolishes_gefitinib_gatekeeper_hbond():
    analyzer = KRASMutantAnalyzer()
    wild_type = analyzer.egfr_gatekeeper_hbond("WT", "gefitinib")
    mutant = analyzer.egfr_gatekeeper_hbond("T790M", "gefitinib")
    assert wild_type.occupancy > 0.50
    assert mutant.abolished
    assert mutant.occupancy < 0.10
    assert mutant.distance_A > 3.5


def test_hbond_gnn_compiles_on_pdbbind_test_complex_and_predicts_affinity():
    model = HBondGNN(hidden_dim=32, message_passing_steps=3)
    graph = model.build_graph("pdbbind_test_protein.pdb", "pdbbind_test_ligand.sdf")
    architecture = model.compile(graph.node_features.shape[1], graph.edge_features.shape[1])
    prediction = model.predict_affinity_detail("pdbbind_test_protein.pdb", "pdbbind_test_ligand.sdf")
    assert architecture["compiled"]
    assert graph.n_hbond_edges >= 1
    assert prediction.graph_nodes == len(graph.nodes)
    assert prediction.graph_edges == len(graph.edges)
    assert prediction.delta_g_kcal_mol < 0.0
    assert prediction.pkd > 0.0


def test_automated_hbond_counts_match_manual_inspection_fixtures():
    p53 = P53MutantAnalyzer()
    kras = KRASMutantAnalyzer()
    assert p53.hbond_count("WT") == MANUAL_COUNTS["p53_wt_2xwr_core"]
    assert p53.hbond_count("R175H") == MANUAL_COUNTS["p53_R175H"]

    kras_count = automated_hbond_count(kras.interaction_fixture("kras_g12c_sotorasib"))
    egfr_wt_count = automated_hbond_count(kras.interaction_fixture("egfr_wt_gefitinib"))
    egfr_mutant_count = automated_hbond_count(kras.interaction_fixture("egfr_t790m_gefitinib"))
    assert kras_count == MANUAL_COUNTS["kras_g12c_sotorasib"]
    assert egfr_wt_count == MANUAL_COUNTS["egfr_wt_gefitinib"]
    assert egfr_mutant_count == MANUAL_COUNTS["egfr_T790M_gefitinib"]


def test_diffdock_style_hbond_reward_prefers_satisfied_pose():
    model = HBondGNN()
    graph = model.build_graph("missing.pdb", "missing.sdf")
    reward = model.hbond_satisfaction_reward(graph)
    pose_score = model.diffusion_docking_pose_score("missing.pdb", "missing.sdf")
    assert reward > 0.0
    assert pose_score["hbond_reward"] == pytest.approx(reward)
    assert np.isfinite(pose_score["pose_score"])
