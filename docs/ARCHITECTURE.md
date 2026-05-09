# HydroOnc Architecture

HydroOnc is organized as a set of focused packages connected by later API, web,
data, and machine-learning layers. Phase 1 lives in `packages/quantum` and
establishes trusted hydrogen physics primitives used by downstream modules.

Phase 2 lives in `packages/proton_transport`. It translates the quantum hydrogen
foundation into biological proton mobility: Grotthuss hopping on water wires,
QM/MM and CP2K validation setup, constant-pH MD configuration, hydrogen-bond
geometry analysis, and pKa benchmark checks.

Phase 3 lives in `packages/tumor_ph`. It models tumor acidosis from Warburg
glycolysis and pH regulators, solves tumor core-to-boundary pH gradients, and
provides CEST-MRI pH mapping surfaces for later API and web visualization work.

Phase 4 lives in `packages/proton_therapy`. It covers Bragg peak physics,
spread-out Bragg peaks, LET/RBE modeling, FLASH dose-rate response, DVH metrics,
ML dose prediction, and PPO-style automated planning interfaces.
