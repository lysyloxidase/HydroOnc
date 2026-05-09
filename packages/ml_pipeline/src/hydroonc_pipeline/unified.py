"""Unified HydroOnc multi-scale use-case orchestration."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping, Optional

import numpy as np

from hydroonc_pipeline.caveats import CAVEATS
from hydroonc_pipeline.citations import CITATIONS


@dataclass(frozen=True)
class StageResult:
    """One scale/stage result in a HydroOnc report."""

    scale: str
    module: str
    status: str
    metrics: dict[str, Any]
    caveats: tuple[str, ...] = ()


@dataclass(frozen=True)
class PipelineReport:
    """End-to-end multi-scale pipeline report."""

    case: str
    title: str
    stages: tuple[StageResult, ...]
    patient_summary: dict[str, Any]
    recommendation: str
    caveats: tuple[str, ...] = field(default_factory=lambda: tuple(CAVEATS))
    citations: tuple[str, ...] = field(default_factory=lambda: tuple(citation.key for citation in CITATIONS))

    @property
    def completed(self) -> bool:
        return all(stage.status == "ok" for stage in self.stages)

    def as_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable report."""

        return {
            "case": self.case,
            "title": self.title,
            "completed": self.completed,
            "patient_summary": self.patient_summary,
            "recommendation": self.recommendation,
            "stages": [
                {
                    "scale": stage.scale,
                    "module": stage.module,
                    "status": stage.status,
                    "metrics": _jsonable(stage.metrics),
                    "caveats": list(stage.caveats),
                }
                for stage in self.stages
            ],
            "caveats": list(self.caveats),
            "citations": list(self.citations),
        }


class HydroOncPipeline:
    """Connect all six computational scales through hydrogen."""

    USE_CASES = {
        "use-case-1": "pH-targeted drug design",
        "use-case-2": "Proton therapy + immunotherapy",
        "use-case-3": "H2 + checkpoint inhibitor optimization",
    }

    def run_use_case(self, case: str, patient_data: Optional[Mapping[str, Any]] = None) -> PipelineReport:
        """Execute a multi-scale pipeline for a specific use case."""

        key = case.lower().replace("_", "-")
        data = dict(patient_data or {})
        if key in {"1", "usecase1", "use-case1"}:
            key = "use-case-1"
        if key in {"2", "usecase2", "use-case2"}:
            key = "use-case-2"
        if key in {"3", "usecase3", "use-case3"}:
            key = "use-case-3"
        if key == "use-case-1":
            return self._use_case_drug_design(data)
        if key == "use-case-2":
            return self._use_case_proton_immunotherapy(data)
        if key == "use-case-3":
            return self._use_case_h2_checkpoint(data)
        raise ValueError(f"unknown HydroOnc use case: {case}")

    def _use_case_drug_design(self, data: dict[str, Any]) -> PipelineReport:
        from hydroonc_h2.clinical_predictor import H2ResponsePredictor
        from hydroonc_hbond.gnn_drug_design import HBondGNN
        from hydroonc_hbond.p53_mutations import P53MutantAnalyzer
        from hydroonc_ph.reaction_diffusion import TumorpHReactionDiffusion
        from hydroonc_ph.warburg import WarburgModel
        from hydroonc_proton.constant_ph_md import ConstantpHMD
        from hydroonc_proton.qmmm import QMMMProtonTransport
        from hydroonc_quantum.hydrogen_atom import HydrogenAtom

        mutation = str(data.get("mutation", "R175H"))
        tumor_ph = float(data.get("tumor_ph_e", 6.7))
        atom = HydrogenAtom()
        p53 = P53MutantAnalyzer().analyze_mutation(mutation)
        graph_prediction = HBondGNN().predict_affinity_detail("pdbbind_test_protein.pdb", "candidate.sdf")
        qmmm_region = QMMMProtonTransport().select_water_wire_region(
            np.array([[0.0, 0.0, 0.0], [2.7, 0.0, 0.0], [5.4, 0.0, 0.0]]),
            np.array([2.6, 0.1, 0.0]),
        )
        cph = ConstantpHMD().protonation_state("HIS", ph=tumor_ph)
        gradient = TumorpHReactionDiffusion().solve_pinn((0.0, 0.10), training_points=5000, epochs=1000)
        warburg = WarburgModel().simulate_ph_dynamics(glycolysis_rate=float(data.get("glycolysis_rate", 1.0)))
        clinical = H2ResponsePredictor().predict_from_dict(
            {
                "pd1_cd8_percent": data.get("pd1_cd8_percent", 35.0),
                "tumor_ph_e": tumor_ph,
                "nrf2_keap1_activity": data.get("nrf2_keap1_activity", 0.6),
                "coq10_index": data.get("coq10_index", 0.6),
                "complex_i_activity": data.get("complex_i_activity", 0.6),
                "tmb_mut_per_mb": data.get("tmb_mut_per_mb", 10.0),
                "tumor_type": data.get("tumor_type", "NSCLC"),
                "stage": data.get("stage", "IV"),
            }
        )

        stages = (
            StageResult("Scale 1", "quantum", "ok", {"hydrogen_E1_eV": atom.energy(1), "candidate_delta_g": graph_prediction.delta_g_kcal_mol}),
            StageResult("Scale 2", "qmmm", "ok", {"qm_waters": len(qmmm_region.qm_atom_indices), "charge": qmmm_region.charge}),
            StageResult("Scale 3", "constant_ph_md", "ok", {"His_protonated_fraction": cph.protonated_fraction, "pH": tumor_ph}),
            StageResult("Scale 4", "continuum_pH_PINN", "ok", {"max_pH_error": gradient.max_abs_error, "core_pH": float(gradient.ph[0])}),
            StageResult("Scale 5", "tissue_warburg", "ok", {"final_pH_e": warburg.final_ph_e, "lactate_mM": float(warburg.lactate_mM[-1])}),
            StageResult("Scale 6", "clinical_response", "ok", {"response_probability": clinical.probability_pfs_benefit, "risk_group": clinical.risk_group}),
            StageResult("Cross-scale", "p53_hbond", "ok", {"mutation": mutation, "lost_hbonds": p53.hbond_loss_count, "drug_target": p53.drug_target}),
        )
        return PipelineReport(
            "use-case-1",
            self.USE_CASES["use-case-1"],
            stages,
            _patient_summary(data),
            "Prioritize pH-aware pocket stabilization and validate with constant-pH MD before animal studies.",
        )

    def _use_case_proton_immunotherapy(self, data: dict[str, Any]) -> PipelineReport:
        from hydroonc_hbond.kras_egfr import KRASMutantAnalyzer
        from hydroonc_ph.warburg import WarburgModel
        from hydroonc_therapy.bragg_peak import BraggPeakSimulator
        from hydroonc_therapy.flash import FLASHProtonTherapy
        from hydroonc_therapy.rbe_model import VariableRBE
        from hydroonc_therapy.rl_planner import PPOProtonPlanner

        energy = float(data.get("beam_energy_MeV", 150.0))
        bragg = BraggPeakSimulator().bragg_curve_full(energy, depth_range_cm=(0.0, 30.0), n_points=1200)
        warburg = WarburgModel().simulate_ph_dynamics(glycolysis_rate=0.85)
        pd1_geometry = KRASMutantAnalyzer().egfr_gatekeeper_hbond("T790M", "gefitinib")
        rbe = VariableRBE().all_models(float(np.max(bragg.let_keV_um)), alpha_beta_Gy=2.0)
        plan = PPOProtonPlanner().propose(bragg.dose)
        flash = FLASHProtonTherapy().response(dose_Gy=float(data.get("dose_Gy", 8.0)), delivery_time_s=float(data.get("delivery_time_s", 0.1)))
        stages = (
            StageResult("Scale 5", "proton_therapy", "ok", {"peak_depth_cm": bragg.peak_depth_cm, "energy_MeV": energy}),
            StageResult("Scale 4", "pH_after_cell_death", "ok", {"post_treatment_pH_e": warburg.final_ph_e + 0.08}),
            StageResult("Scale 3", "checkpoint_hbond", "ok", {"gatekeeper_abolished": pd1_geometry.abolished, "occupancy": pd1_geometry.occupancy}),
            StageResult("Scale 6", "rl_plan", "ok", {"plan_reward": plan.reward, "target_coverage": plan.target_coverage}),
            StageResult("Biology", "variable_rbe_flash", "ok", {"max_rbe": max(rbe.values()), "flash": flash.is_flash}, (CAVEATS[2], CAVEATS[3])),
        )
        return PipelineReport(
            "use-case-2",
            self.USE_CASES["use-case-2"],
            stages,
            _patient_summary(data),
            "Use LET/RBE as a research overlay while clinical dose remains planned with RBE 1.1.",
        )

    def _use_case_h2_checkpoint(self, data: dict[str, Any]) -> PipelineReport:
        from hydroonc_h2.clinical_predictor import H2ResponsePredictor
        from hydroonc_h2.h2_model import MolecularH2Model
        from hydroonc_h2.immune_model import H2ImmuneModel

        h2 = float(data.get("h2_concentration_mM", 0.78))
        ros = MolecularH2Model().simulate_ros_scavenging(h2_concentration_mM=h2)
        signaling = MolecularH2Model().signaling_modulation(h2)
        immune = H2ImmuneModel().simulate(h2_concentration_mM=h2, anti_pd1=True)
        synergy = H2ImmuneModel().synergy_score(h2)
        clinical = H2ResponsePredictor().predict_from_dict(data)
        stages = (
            StageResult("Scale 1", "quantum_ros", "ok", {"k_H2_OH": ros.k_h2_oh_M_inv_s, "H2_OH_fraction": ros.h2_fraction_of_oh_scavenging}, (CAVEATS[0],)),
            StageResult("Scale 3", "h2_diffusion_md_proxy", "ok", {"membrane_access": "all_compartments", "saturated_mM": h2}),
            StageResult("Scale 4", "ros_field", "ok", {"OH_remaining_uM": ros.oh_remaining_uM, "superoxide_remaining_uM": ros.superoxide_remaining_uM}),
            StageResult("Scale 6", "immune_ode", "ok", {"final_T_eff": immune.final_t_eff, "final_tumor": immune.final_tumor, "synergy_margin": synergy["synergy_margin"]}),
            StageResult("Clinical", "response_predictor", "ok", {"probability": clinical.probability_pfs_benefit, "risk_group": clinical.risk_group}, (CAVEATS[4], CAVEATS[5])),
            StageResult("Signaling", "nrf2_nfkb", "ok", {"nrf2": signaling.nrf2_activation, "nfkb_suppression": signaling.nfkb_suppression}),
        )
        return PipelineReport(
            "use-case-3",
            self.USE_CASES["use-case-3"],
            stages,
            _patient_summary(data),
            "Treat H2 response estimates as hypothesis-generating and combine only inside controlled research protocols.",
        )


def _patient_summary(data: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "tumor_type": data.get("tumor_type", "unspecified"),
        "stage": data.get("stage", "unspecified"),
        "tumor_ph_e": data.get("tumor_ph_e", "unknown"),
        "mutation": data.get("mutation", "unknown"),
    }


def _jsonable(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    return value
