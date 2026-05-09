# Research Report

Living report for HydroOnc assumptions, sources, benchmarks, and validation
results.

## Phase 2 Proton Transport

Reference proton diffusion in room-temperature water is
`D_H+ ~= 9.31e-5 cm^2/s`, roughly four times water self-diffusion. HydroOnc
models this with a structural diffusion fraction of 0.90, one hop per
picosecond, 2.65 Angstrom effective water-wire spacing, and a correlation factor
chosen to reproduce the experimental value within the Phase 2 quality gate.

Tumor extracellular pH is represented by `pH = 6.7`; normal extracellular
comparisons use `pH = 7.4`. Histidine protonation therefore changes strongly
across the simulated tumor-normal pH contrast.
