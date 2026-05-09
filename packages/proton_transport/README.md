# HydroOnc Proton Transport

Phase 2 models proton transport in biological water networks and tumor
microenvironments.

Implemented components:

- simplified MS-EVB-style Grotthuss hopping on water wires
- diffusion calibration to the experimental proton diffusion coefficient in water
- CP2K input generation for AIMD validation
- QM/MM region and electrostatic-embedding configuration helpers
- constant-pH MD setup with an OpenMM adapter and deterministic mock backend
- residue protonation-state sampling for tumor and normal pH ranges
- hydrogen-bond geometry analysis for DNA and protein motifs
- PROPKA3 adapter with lysozyme benchmark fallback values

Install for development:

```bash
pip install -e "packages/proton_transport[dev]"
```

Optional MD and pKa integrations:

```bash
pip install -e "packages/proton_transport[dev,md,propka]"
```
