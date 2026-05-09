# HydroOnc Proton Therapy

Phase 4 models clinical proton therapy physics and planning:

- Bethe-Bloch mass stopping power with PSTAR-style water reference values
- Bragg-Kleeman range-energy mapping and pristine Bragg curves
- spread-out Bragg peak (SOBP) target coverage
- LET and variable-RBE models for distal-edge risk
- FLASH dose-rate and oxygen-depletion response
- lightweight 3D U-Net-compatible dose prediction surface
- DVH computation and PPO-style planning helpers

Install for development:

```bash
pip install -e "packages/proton_therapy[dev]"
```

Optional ML backend:

```bash
pip install -e "packages/proton_therapy[dev,ml]"
```
