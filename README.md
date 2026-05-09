# HydroOnc

**From Schrödinger to the clinic — hydrogen at every scale of cancer.**

HydroOnc is a multi-scale hydrogen-in-cancer research platform. It connects
quantum hydrogen physics, proton transport, tumor pH, proton therapy, mutant
oncoprotein H-bonds, molecular hydrogen therapy, and a unified clinical-facing
pipeline. It is not clinical decision support.

## Repository Layout

```text
hydroonc/
├── packages/quantum          # Phase 1: Schrodinger equation, orbitals, H2
├── packages/proton_transport # Phase 2: Grotthuss, QM/MM, constant-pH MD
├── packages/tumor_ph         # Phase 3: Warburg, pH PDE, CEST-ML
├── packages/proton_therapy   # Phase 4: Bragg peak, RBE, FLASH, dose ML
├── packages/hbond_onco       # Phase 5: mutant protein H-bonds, GNN design
├── packages/h2_therapy       # Phase 6: H2 ROS, immune ODE, response ML
├── packages/ml_pipeline      # Phase 7: unified pipeline + CLI
├── apps/api                  # FastAPI
├── apps/web                  # Next.js + R3F visualization
├── data
└── docs
```

## Six Scales

```text
Scale 1 Quantum      H atom, H2, orbitals, spectra
Scale 2 QM/MM        Grotthuss, water wires, proton transfer setup
Scale 3 MD           constant-pH MD, H-bond occupancy, immune ODEs
Scale 4 Continuum    tumor pH reaction-diffusion and PINN surrogates
Scale 5 Tissue       Bragg peak, SOBP, LET/RBE, FLASH, DVH
Scale 6 Clinical     response prediction, RL planning, multi-scale reports
```

## Phase 1 Quick Start

```bash
python -m pip install -e "packages/quantum[dev]"
pytest packages/quantum
```

## Phase 2 Quick Start

```bash
python -m pip install -e "packages/proton_transport[dev]"
pytest packages/proton_transport
```

## Phase 3 Quick Start

```bash
python -m pip install -e "packages/tumor_ph[dev]"
pytest packages/tumor_ph
```

## Phase 4 Quick Start

```bash
python -m pip install -e "packages/proton_therapy[dev]"
pytest packages/proton_therapy
```

## Phase 5 Quick Start

```bash
python -m pip install -e "packages/hbond_onco[dev]"
pytest packages/hbond_onco
```

## Phase 6 Quick Start

```bash
python -m pip install -e "packages/h2_therapy[dev]"
pytest packages/h2_therapy
```

## Phase 7 Quick Start

```bash
python -m pip install -e "packages/ml_pipeline[dev]"
hydroonc pipeline use-case-1
pytest packages/ml_pipeline
```

Open the dependency-free visualization shell:

```text
apps/web/static/index.html
```

Optional molecular and visualization backends:

```bash
python -m pip install -e "packages/quantum[dev,h2,viz]"
```

## Use Cases

**pH-targeted drug design**

```text
Quantum H-bond energetics -> QM/MM proton transfer -> constant-pH MD ->
tumor pH PINN -> response prediction from CEST pH and molecular features
```

**Proton therapy + immunotherapy**

```text
Bragg peak and LET/RBE -> pH shift after cell death -> checkpoint H-bonds ->
PPO plan score and anti-PD-1 timing hypothesis
```

**H2 + checkpoint inhibitor optimization**

```text
H2 ROS kinetics -> tumor diffusion proxy -> ROS field -> immune ODE ->
research-grade response probability
```

## Key Results

- Quantum: H atom `E_n` matches Bohr values to 4 decimals.
- Tumor pH: Warburg model reproduces tumor `pH_e` around `6.4-7.0`.
- Proton therapy: Bragg/stopping references match PSTAR-style tables within 2%.
- H-bond oncology: p53 R175H loses 4+ H-bonds versus wild-type 2XWR.
- H2 immune: ODE predicts H2 plus anti-PD-1 outperforms monotherapy.
- GNN: H-bond-aware affinity prediction runs on a PDBbind-style test complex.

## Physical Constants

| Constant | Value | Source |
| --- | ---: | --- |
| Electron mass | `9.1093837015e-31 kg` | CODATA 2018 |
| Elementary charge | `1.602176634e-19 C` | CODATA 2018 |
| Reduced Planck constant | `1.054571817e-34 J s` | CODATA 2018 |
| Speed of light | `299792458 m/s` | CODATA 2018 |
| Bohr radius | `5.29177e-11 m` | derived CODATA 2018 |
| Fine-structure constant | `~1/137.036` | derived CODATA 2018 |

## Caveats

- H2 selective scavenging is CONTESTED: physiological H2 contributes less than 1% of hydroxyl quenching with GSH present.
- APR-246 Phase 3 FAILED despite Phase 2 promise.
- Variable RBE is not yet endorsed by ICRU; RBE = 1.1 remains standard.
- FLASH FAST-01 was n=10 feasibility, not efficacy.
- H2 cancer trials are small, n<200 total, and NOT Phase 3.
- HydroOnc is a RESEARCH platform, not clinical decision support.

## Visualization Pages

- Quantum Explorer: orbitals, spectra, H2 PEC, TDSE.
- Proton Transport: Grotthuss hopping, H-bond networks, pH gradient.
- Tumor pH Dashboard: Warburg, reaction-diffusion, CEST overlay.
- Proton Therapy: Bragg, SOBP, dose, RBE, FLASH.
- Oncoprotein H-Bonds: p53/KRAS/EGFR difference maps and GNN affinity.
- H2 Therapy: ROS dynamics, immune ODE, clinical response predictor.
