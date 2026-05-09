# HydroOnc

HydroOnc is a multi-scale hydrogen-in-cancer platform. Phase 1 builds the
quantum-mechanical foundation for hydrogen atom physics, H2 electronic
structure, time-dependent wave mechanics, perturbations, spectra, and orbital
visualization. Phase 2 adds proton transport across biological water wires and
tumor-relevant pH contexts. Phase 3 models Warburg-driven tumor acidosis,
pH regulators, reaction-diffusion gradients, and CEST-MRI pH mapping. Phase 4
adds proton therapy physics, FLASH, variable RBE, DVH, and ML dose prediction.
Phase 5 analyzes mutant oncoprotein H-bond disruption and GNN-assisted drug
design. Phase 6 models molecular hydrogen therapy through ROS caveats, immune
restoration, and response prediction.

## Repository Layout

```text
hydroonc/
├── packages/quantum          # Phase 1: Schrodinger equation, orbitals, H2
├── packages/proton_transport # Phase 2: Grotthuss, QM/MM, constant-pH MD
├── packages/tumor_ph         # Phase 3: Warburg, pH PDE, CEST-ML
├── packages/proton_therapy   # Phase 4: Bragg peak, RBE, FLASH, dose ML
├── packages/hbond_onco       # Phase 5: mutant protein H-bonds, GNN design
├── packages/h2_therapy       # Phase 6: H2 ROS, immune ODE, response ML
├── packages/ml_pipeline      # Phase 7
├── apps/api                  # FastAPI
├── apps/web                  # Next.js + R3F visualization
├── data
└── docs
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

Optional molecular and visualization backends:

```bash
python -m pip install -e "packages/quantum[dev,h2,viz]"
```
