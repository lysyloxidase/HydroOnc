"""HydroOnc command-line interface."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Callable, Optional

import numpy as np

from hydroonc_pipeline.caveats import CAVEATS
from hydroonc_pipeline.citations import CITATIONS
from hydroonc_pipeline.unified import HydroOncPipeline

CommandHandler = Callable[[argparse.Namespace], dict[str, Any]]


def main(argv: Optional[list[str]] = None) -> int:
    """Run the HydroOnc CLI."""

    parser = build_parser()
    args = parser.parse_args(argv)
    if not hasattr(args, "handler"):
        parser.print_help()
        return 0
    result = args.handler(args)
    print(json.dumps(_jsonable(result), indent=2, sort_keys=True))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="hydroonc", description="HydroOnc multi-scale hydrogen-in-cancer research CLI")
    sub = parser.add_subparsers(dest="domain")
    _quantum(sub.add_parser("quantum"))
    _proton(sub.add_parser("proton"))
    _ph(sub.add_parser("ph"))
    _therapy(sub.add_parser("therapy"))
    _hbond(sub.add_parser("hbond"))
    _h2(sub.add_parser("h2"))
    _pipeline(sub.add_parser("pipeline"))
    meta = sub.add_parser("meta")
    meta_sub = meta.add_subparsers(dest="meta_command")
    _leaf(meta_sub, "caveats", _cmd_caveats)
    _leaf(meta_sub, "citations", _cmd_citations)
    return parser


def _quantum(parser: argparse.ArgumentParser) -> None:
    sub = parser.add_subparsers(dest="quantum_command")
    p = _leaf(sub, "orbital", _cmd_quantum_orbital)
    p.add_argument("--n", type=int, default=1)
    p.add_argument("--l", type=int, default=0)
    p.add_argument("--m", type=int, default=0)
    p.add_argument("--visualize", action="store_true")
    p = _leaf(sub, "spectrum", _cmd_quantum_spectrum)
    p.add_argument("--series", default="balmer")
    p = _leaf(sub, "h2", _cmd_quantum_h2)
    p.add_argument("--method", default="ccsd")
    p.add_argument("--basis", default="cc-pVTZ")
    p.add_argument("--scan", action="store_true")
    p = _leaf(sub, "energy", _cmd_quantum_energy)
    p.add_argument("--n", type=int, default=1)


def _proton(parser: argparse.ArgumentParser) -> None:
    sub = parser.add_subparsers(dest="proton_command")
    p = _leaf(sub, "grotthuss", _cmd_proton_grotthuss)
    p.add_argument("--n-waters", type=int, default=20)
    p.add_argument("--temperature", type=float, default=300.0)
    p = _leaf(sub, "cph", _cmd_proton_cph)
    p.add_argument("--residue", default="HIS")
    p.add_argument("--ph", type=float, default=6.7)
    p = _leaf(sub, "hbond", _cmd_proton_hbond)
    p.add_argument("--base-pair", default="A-T")
    p = _leaf(sub, "qmmm", _cmd_proton_qmmm)
    p.add_argument("--radius", type=float, default=3.0)


def _ph(parser: argparse.ArgumentParser) -> None:
    sub = parser.add_subparsers(dest="ph_command")
    p = _leaf(sub, "warburg", _cmd_ph_warburg)
    p.add_argument("--glycolysis-rate", type=float, default=1.0)
    p.add_argument("--inhibit", default="")
    p.add_argument("--dose", default="0uM")
    p = _leaf(sub, "pinn", _cmd_ph_pinn)
    p.add_argument("--mesh", default="")
    p.add_argument("--train", action="store_true")
    p.add_argument("--epochs", type=int, default=5000)
    p = _leaf(sub, "cest", _cmd_ph_cest)
    p.add_argument("--features", type=int, default=64)
    p = _leaf(sub, "inhibitor", _cmd_ph_inhibitor)
    p.add_argument("--target", default="CAIX")
    p.add_argument("--dose", default="10uM")


def _therapy(parser: argparse.ArgumentParser) -> None:
    sub = parser.add_subparsers(dest="therapy_command")
    p = _leaf(sub, "bragg", _cmd_therapy_bragg)
    p.add_argument("--energy", type=float, default=150.0)
    p.add_argument("--material", default="water")
    p = _leaf(sub, "sobp", _cmd_therapy_sobp)
    p.add_argument("--target-depth", default="5-10")
    p.add_argument("--optimize", action="store_true")
    p = _leaf(sub, "dose-predict", _cmd_therapy_dose_predict)
    p.add_argument("--ct", default="")
    p.add_argument("--model", default="unet3d")
    p = _leaf(sub, "rl-plan", _cmd_therapy_rl_plan)
    p.add_argument("--ct", default="")
    p.add_argument("--algorithm", default="ppo")
    p = _leaf(sub, "flash", _cmd_therapy_flash)
    p.add_argument("--dose", type=float, default=8.0)
    p.add_argument("--time", type=float, default=0.1)
    p = _leaf(sub, "rbe", _cmd_therapy_rbe)
    p.add_argument("--let", type=float, default=9.0)
    p = _leaf(sub, "dvh", _cmd_therapy_dvh)
    p.add_argument("--dose", type=float, default=70.0)


def _hbond(parser: argparse.ArgumentParser) -> None:
    sub = parser.add_subparsers(dest="hbond_command")
    p = _leaf(sub, "p53", _cmd_hbond_p53)
    p.add_argument("--mutation", default="R175H")
    p.add_argument("--compare-wt", action="store_true")
    p.add_argument("--md-ns", type=float, default=100.0)
    p = _leaf(sub, "kras", _cmd_hbond_kras)
    p.add_argument("--mutation", default="G12C")
    p.add_argument("--drug", default="sotorasib")
    p.add_argument("--hbond-analysis", action="store_true")
    p = _leaf(sub, "egfr", _cmd_hbond_egfr)
    p.add_argument("--mutation", default="T790M")
    p.add_argument("--drug", default="gefitinib")
    p = _leaf(sub, "gnn", _cmd_hbond_gnn)
    p.add_argument("--protein", default="target.pdb")
    p.add_argument("--ligand", default="candidate.sdf")


def _h2(parser: argparse.ArgumentParser) -> None:
    sub = parser.add_subparsers(dest="h2_command")
    p = _leaf(sub, "ros", _cmd_h2_ros)
    p.add_argument("--h2-conc", type=float, default=0.78)
    p.add_argument("--oh-conc", type=float, default=1.0)
    p.add_argument("--simulate", action="store_true")
    p = _leaf(sub, "immune", _cmd_h2_immune)
    p.add_argument("--h2-dose", default="3h-day")
    p.add_argument("--anti-pd1", action="store_true")
    p.add_argument("--simulate", action="store_true")
    p = _leaf(sub, "predict", _cmd_h2_predict)
    p.add_argument("--pd1-cd8", type=float, default=42.0)
    p.add_argument("--tumor-ph", type=float, default=6.7)


def _pipeline(parser: argparse.ArgumentParser) -> None:
    sub = parser.add_subparsers(dest="pipeline_command")
    for name in ["use-case-1", "use-case-2", "use-case-3"]:
        p = _leaf(sub, name, _cmd_pipeline)
        p.add_argument("--patient-data", default="")
        p.set_defaults(case=name)


def _leaf(sub: Any, name: str, handler: CommandHandler) -> argparse.ArgumentParser:
    parser = sub.add_parser(name)
    parser.set_defaults(handler=handler)
    return parser


def _cmd_quantum_orbital(args: argparse.Namespace) -> dict[str, Any]:
    from hydroonc_quantum.visualization import classify_orbital_shape

    return {"command": "quantum orbital", "n": args.n, "l": args.l, "m": args.m, "shape": classify_orbital_shape(args.n, args.l, args.m), "visualize": args.visualize}


def _cmd_quantum_spectrum(args: argparse.Namespace) -> dict[str, Any]:
    from hydroonc_quantum.spectrum import HydrogenSpectrum

    lines = HydrogenSpectrum().series(args.series.capitalize(), n_max={"lyman": 4, "balmer": 5, "paschen": 6}.get(args.series.lower(), 5))
    return {"command": "quantum spectrum", "series": args.series, "wavelengths_nm": [line.wavelength_nm for line in lines]}


def _cmd_quantum_h2(args: argparse.Namespace) -> dict[str, Any]:
    from hydroonc_quantum.h2_molecule import H2Molecule

    h2 = H2Molecule().compute_energy(method=args.method, basis=args.basis)
    return {"command": "quantum h2", "scan": args.scan, **h2}


def _cmd_quantum_energy(args: argparse.Namespace) -> dict[str, Any]:
    from hydroonc_quantum.hydrogen_atom import HydrogenAtom

    return {"command": "quantum energy", "n": args.n, "energy_eV": HydrogenAtom().energy(args.n)}


def _cmd_proton_grotthuss(args: argparse.Namespace) -> dict[str, Any]:
    from hydroonc_proton.grotthuss import GrotthussSimulator

    sim = GrotthussSimulator(n_waters=args.n_waters)
    traj = sim.simulate(duration_ps=20.0)
    return {"command": "proton grotthuss", "hops": traj.n_hops, "rate_per_ps": traj.hopping_rate_per_ps, "diffusion_cm2_s": sim.diffusion_coefficient_cm2_s(), "temperature_K": args.temperature}


def _cmd_proton_cph(args: argparse.Namespace) -> dict[str, Any]:
    from hydroonc_proton.constant_ph_md import ConstantpHMD

    state = ConstantpHMD().protonation_state(args.residue, args.ph)
    return {"command": "proton cph", **state.__dict__}


def _cmd_proton_hbond(args: argparse.Namespace) -> dict[str, Any]:
    from hydroonc_proton.hbond_analysis import HBondAnalyzer

    bonds = HBondAnalyzer().analyze_base_pair(args.base_pair)
    return {"command": "proton hbond", "base_pair": args.base_pair, "count": len(bonds)}


def _cmd_proton_qmmm(args: argparse.Namespace) -> dict[str, Any]:
    from hydroonc_proton.qmmm import QMMMProtonTransport

    region = QMMMProtonTransport().select_water_wire_region(np.array([[0, 0, 0], [2.7, 0, 0], [5.4, 0, 0]], dtype=float), np.array([2.6, 0.1, 0.0]), radius_A=args.radius)
    return {"command": "proton qmmm", "qm_atom_indices": list(region.qm_atom_indices), "charge": region.charge}


def _cmd_ph_warburg(args: argparse.Namespace) -> dict[str, Any]:
    from hydroonc_ph.warburg import WarburgModel

    model = WarburgModel()
    if args.inhibit:
        return {"command": "ph warburg", **model.inhibitor_effect(args.inhibit, _parse_uM(args.dose))}
    series = model.simulate_ph_dynamics(glycolysis_rate=args.glycolysis_rate)
    return {"command": "ph warburg", "final_ph_e": series.final_ph_e, "final_ph_i": series.final_ph_i}


def _cmd_ph_pinn(args: argparse.Namespace) -> dict[str, Any]:
    from hydroonc_ph.reaction_diffusion import TumorpHReactionDiffusion

    sol = TumorpHReactionDiffusion().solve_pinn((0.0, 0.10), epochs=args.epochs)
    return {"command": "ph pinn", "trained": args.train, "mesh": args.mesh, "max_abs_error": sol.max_abs_error}


def _cmd_ph_cest(args: argparse.Namespace) -> dict[str, Any]:
    from hydroonc_ph.cest_ml import CESTpHMapper

    offsets = np.linspace(-5.0, 5.0, args.features)
    z = 0.90 - 0.08 * np.exp(-0.5 * ((offsets - 3.5) / 0.8) ** 2)
    pred = CESTpHMapper().predict(z)
    return {"command": "ph cest", "predicted_ph": float(pred[0])}


def _cmd_ph_inhibitor(args: argparse.Namespace) -> dict[str, Any]:
    from hydroonc_ph.warburg import WarburgModel

    return {"command": "ph inhibitor", **WarburgModel().inhibitor_effect(args.target, _parse_uM(args.dose))}


def _cmd_therapy_bragg(args: argparse.Namespace) -> dict[str, Any]:
    from hydroonc_therapy.bragg_peak import BraggPeakSimulator

    curve = BraggPeakSimulator().bragg_curve_full(args.energy)
    return {"command": "therapy bragg", "energy_MeV": args.energy, "material": args.material, "peak_depth_cm": curve.peak_depth_cm}


def _cmd_therapy_sobp(args: argparse.Namespace) -> dict[str, Any]:
    from hydroonc_therapy.bragg_peak import BraggPeakSimulator

    target = _parse_range(args.target_depth)
    sobp = BraggPeakSimulator().sobp_full(target)
    return {"command": "therapy sobp", "target_depth_cm": list(target), "uniformity": sobp.target_uniformity(), "optimize": args.optimize}


def _cmd_therapy_dose_predict(args: argparse.Namespace) -> dict[str, Any]:
    from hydroonc_therapy.dose_prediction import DosePredictionUNet

    model = DosePredictionUNet()
    arch = model.compile((32, 32, 32, 3))
    return {"command": "therapy dose-predict", "model": args.model, "ct": args.ct, "compiled": arch.compiled, "output_shape": list(arch.output_shape)}


def _cmd_therapy_rl_plan(args: argparse.Namespace) -> dict[str, Any]:
    from hydroonc_therapy.rl_planner import PPOProtonPlanner

    result = PPOProtonPlanner().propose(np.linspace(0, 1, 64))
    return {"command": "therapy rl-plan", "algorithm": args.algorithm, "reward": result.reward, "weights": result.objective_weights}


def _cmd_therapy_flash(args: argparse.Namespace) -> dict[str, Any]:
    from hydroonc_therapy.flash import FLASHProtonTherapy

    response = FLASHProtonTherapy().response(args.dose, args.time)
    return {"command": "therapy flash", **response.__dict__}


def _cmd_therapy_rbe(args: argparse.Namespace) -> dict[str, Any]:
    from hydroonc_therapy.rbe_model import VariableRBE

    return {"command": "therapy rbe", "models": VariableRBE().all_models(args.let, alpha_beta_Gy=2.0)}


def _cmd_therapy_dvh(args: argparse.Namespace) -> dict[str, Any]:
    from hydroonc_therapy.dvh import compute_dvh

    dose = np.broadcast_to(np.linspace(0, args.dose, 16)[None, None, :], (16, 16, 16))
    dvh = compute_dvh(dose, np.ones_like(dose, dtype=bool))
    return {"command": "therapy dvh", "v35": dvh.v_at_dose(args.dose / 2.0), "d50": dvh.d_at_volume(0.5)}


def _cmd_hbond_p53(args: argparse.Namespace) -> dict[str, Any]:
    from hydroonc_hbond.p53_mutations import P53MutantAnalyzer

    comparison = P53MutantAnalyzer().analyze_mutation(args.mutation)
    return {"command": "hbond p53", "compare_wt": args.compare_wt, "md_ns": args.md_ns, **comparison.__dict__}


def _cmd_hbond_kras(args: argparse.Namespace) -> dict[str, Any]:
    from hydroonc_hbond.kras_egfr import KRASMutantAnalyzer

    result = KRASMutantAnalyzer().kras_hbond_geometry(args.mutation, args.drug)
    return {"command": "hbond kras", "hbond_analysis": args.hbond_analysis, **result.__dict__}


def _cmd_hbond_egfr(args: argparse.Namespace) -> dict[str, Any]:
    from hydroonc_hbond.kras_egfr import KRASMutantAnalyzer

    result = KRASMutantAnalyzer().egfr_gatekeeper_hbond(args.mutation, args.drug)
    return {"command": "hbond egfr", **result.__dict__}


def _cmd_hbond_gnn(args: argparse.Namespace) -> dict[str, Any]:
    from hydroonc_hbond.gnn_drug_design import HBondGNN

    prediction = HBondGNN().predict_affinity_detail(args.protein, args.ligand)
    return {"command": "hbond gnn", **prediction.__dict__}


def _cmd_h2_ros(args: argparse.Namespace) -> dict[str, Any]:
    from hydroonc_h2.h2_model import MolecularH2Model

    result = MolecularH2Model().simulate_ros_scavenging(args.h2_conc, args.oh_conc)
    return {"command": "h2 ros", "simulate": args.simulate, **result.__dict__}


def _cmd_h2_immune(args: argparse.Namespace) -> dict[str, Any]:
    from hydroonc_h2.immune_model import H2ImmuneModel

    h2 = 0.78 if args.h2_dose == "3h-day" else _parse_uM(args.h2_dose) / 1000.0
    trajectory = H2ImmuneModel().simulate(h2_concentration_mM=h2, anti_pd1=args.anti_pd1)
    return {"command": "h2 immune", "simulate": args.simulate, "final_t_eff": trajectory.final_t_eff, "final_tumor": trajectory.final_tumor, "anti_pd1": args.anti_pd1}


def _cmd_h2_predict(args: argparse.Namespace) -> dict[str, Any]:
    from hydroonc_h2.clinical_predictor import H2ResponsePredictor

    prediction = H2ResponsePredictor().predict_from_dict({"pd1_cd8_percent": args.pd1_cd8, "tumor_ph_e": args.tumor_ph})
    return {"command": "h2 predict", **prediction.__dict__}


def _cmd_pipeline(args: argparse.Namespace) -> dict[str, Any]:
    data = _load_patient(args.patient_data)
    return HydroOncPipeline().run_use_case(args.case, data).as_dict()


def _cmd_caveats(args: argparse.Namespace) -> dict[str, Any]:
    del args
    return {"command": "meta caveats", "caveats": CAVEATS}


def _cmd_citations(args: argparse.Namespace) -> dict[str, Any]:
    del args
    return {"command": "meta citations", "citations": [citation.__dict__ for citation in CITATIONS]}


def _parse_uM(value: str) -> float:
    text = str(value).strip().lower().replace("µ", "u")
    if text.endswith("um"):
        return float(text[:-2])
    if text.endswith("mm"):
        return float(text[:-2]) * 1000.0
    return float(text)


def _parse_range(value: str) -> tuple[float, float]:
    left, right = value.replace(",", "-").split("-", 1)
    return float(left), float(right)


def _load_patient(path: str) -> dict[str, Any]:
    if not path:
        return {}
    file_path = Path(path)
    if not file_path.exists():
        return {}
    return json.loads(file_path.read_text(encoding="utf-8"))


def _jsonable(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, tuple):
        return [_jsonable(item) for item in value]
    if isinstance(value, list):
        return [_jsonable(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    return value


if __name__ == "__main__":
    raise SystemExit(main())
