# Proton Transport

Phase 2 models why proton mobility in water is anomalously fast. In addition to
vehicular hydronium diffusion, excess charge moves structurally through
Grotthuss hopping along hydrogen-bonded water wires.

The package implements:

- a simplified MS-EVB-style water-wire model calibrated to
  `D_H+ = 9.31e-5 cm^2/s`
- Eigen/Zundel carrier-state hopping trajectories near `1 ps^-1`
- CP2K AIMD input generation for BLYP-D3 or RPBE-D3, 0.5 fs, NVT 300 K
- QM/MM region selection around excess-proton positions
- tumor and normal constant-pH MD setup surfaces
- Henderson-Hasselbalch protonation-state estimates for His, Asp, and Glu
- Watson-Crick and protein hydrogen-bond geometry analysis
- PROPKA3 lysozyme benchmark fallbacks for CI
