# HydroOnc H-Bond Oncology

Phase 5 analyzes hydrogen-bond disruption in mutant cancer proteins and exposes
a GNN-style drug-design surface.

Implemented components:

- p53 hotspot mutation H-bond occupancy comparisons
- curated R175H, R273H, R248W, G245S, R249S, and R282W impact summaries
- KRAS G12C/D/V switch-II pocket H-bond checks
- EGFR T790M/C797S gatekeeper-resistance analysis
- geometric H-bond counting utilities
- protein-ligand interaction graph construction
- lightweight H-bond-aware affinity predictor and docking reward

Install for development:

```bash
pip install -e "packages/hbond_onco[dev]"
```

Optional ML chemistry stack:

```bash
pip install -e "packages/hbond_onco[dev,ml]"
```
