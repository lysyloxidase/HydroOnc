# HydroOnc Quantum

Phase 1 of HydroOnc implements a compact quantum-mechanical foundation for
hydrogen physics:

- analytical hydrogen atom wavefunctions and energies
- radial probability distributions and orbital grids
- numerical radial solving with a Numerov shooting method
- emission spectra for the Lyman, Balmer, Paschen, Brackett, and Pfund series
- first-order perturbation helpers
- variational hydrogen ground-state estimates
- split-operator time-dependent Schrodinger evolution
- H2 electronic-structure wrappers for PySCF, with explicit reference fallback
- 3D orbital isosurface sampling for visualization

Install for development:

```bash
pip install -e "packages/quantum[dev]"
```

Install optional H2 and visualization backends:

```bash
pip install -e "packages/quantum[dev,h2,viz]"
```
