# HydroOnc

HydroOnc is a multi-scale hydrogen-in-cancer platform. Phase 1 builds the
quantum-mechanical foundation for hydrogen atom physics, H2 electronic
structure, time-dependent wave mechanics, perturbations, spectra, and orbital
visualization. Phase 2 adds proton transport across biological water wires and
tumor-relevant pH contexts.

## Repository Layout

```text
hydroonc/
├── packages/quantum          # Phase 1: Schrodinger equation, orbitals, H2
├── packages/proton_transport # Phase 2: Grotthuss, QM/MM, constant-pH MD
├── packages/tumor_ph         # Phase 3
├── packages/proton_therapy   # Phase 4
├── packages/hbond_onco       # Phase 5
├── packages/h2_therapy       # Phase 6
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

Optional molecular and visualization backends:

```bash
python -m pip install -e "packages/quantum[dev,h2,viz]"
```
