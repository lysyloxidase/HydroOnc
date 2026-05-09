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

## Phase 3 Tumor pH

The Warburg module calibrates aerobic glycolysis and pH-regulator export so
extracellular pH drops from `7.4` toward `6.7` over one simulated hour, while
intracellular pH remains near `7.0-7.2`. NHE1 inhibition lowers intracellular pH
by about `0.2-0.3` units, and CAIX inhibition raises extracellular pH toward the
normal tissue range.

The reaction-diffusion module represents an acidic tumor core around `pH 6.5`
and a vascular boundary around `pH 7.3`. A finite-difference fallback preserves
the FEniCS-facing API for local testing and CI.
