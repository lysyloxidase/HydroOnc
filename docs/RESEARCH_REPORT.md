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

## Phase 4 Proton Therapy

Clinical proton energies span roughly `70-250 MeV`, corresponding to water
ranges from about `4.0 cm` to `37.7 cm`. HydroOnc calibrates the Bragg peak
module to place a `150 MeV` pristine peak near `15.6 cm` in water and uses
PSTAR-style mass stopping power references for quality gates.

The RBE module keeps the ICRU clinical convention of `1.1` as a floor while
showing distal-edge increases for McNamara, Wedenberg, McMahon, and
Carabe-Fernandez-style models. FLASH delivery is represented by dose rates above
`40 Gy/s`, oxygen depletion, and a normal-tissue sparing factor.

## Phase 5 H-Bond Oncology

p53 is mutated in roughly half of human cancers, with hotspot residues including
R175, G245, R248, R249, R273, and R282. HydroOnc models structural mutations as
loss or rearrangement of stabilizing H-bond networks and contact mutations as
loss of direct protein-DNA H-bonds.

KRAS G12C inhibitors depend on the switch-II pocket network, including a His95
H-bond that can be disrupted by resistance mutations such as Y96D/S and H95Q/R.
EGFR T790M is represented as simultaneous steric blocking and loss of the
gefitinib/erlotinib gatekeeper H-bond, with osimertinib bypassing via C797
covalency until C797S resistance.

## Phase 6 Molecular H2 Therapy

Saturated molecular hydrogen in water is represented as `1.57 mg/L`, about
`0.78 mM`. The ROS module uses `k(H2 + OH) = 4.2e7 M^-1 s^-1` and
`k(H2 + ONOO-) = 3.6e4 M^-1 s^-1`, while assigning negligible direct reaction
with superoxide and hydrogen peroxide.

The critical caveat is encoded directly: with physiological GSH around `5 mM`
and `k(GSH + OH) ~= 1e10 M^-1 s^-1`, GSH dominates hydroxyl radical quenching by
more than three orders of magnitude. H2 oncology effects are therefore modeled
primarily through Nrf2/NF-kB/mitochondrial signaling and immune restoration of
exhausted CD8+ T cells, including synergy with anti-PD-1 therapy.
