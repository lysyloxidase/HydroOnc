# Proton Therapy

Phase 4 models clinical proton therapy from beam physics through planning:

- proton ranges for `70-250 MeV` beams in water
- Bethe-Bloch mass stopping power with bundled PSTAR-style references
- pristine Bragg curves and spread-out Bragg peaks
- LET rise at the distal edge and variable-RBE models
- FLASH dose-rate response above `40 Gy/s`
- DVH metrics from 3D dose grids
- 3D U-Net-compatible dose prediction and PPO-style planning helpers

The reference gate for a `150 MeV` proton places the Bragg peak near
`15.6 cm` depth in water. The clinical RBE convention is `1.1`, but implemented
variable-RBE models increase above this value at high distal-edge LET.
