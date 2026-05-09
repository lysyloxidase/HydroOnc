# HydroOnc Architecture

HydroOnc is organized as a set of focused packages connected by later API, web,
data, and machine-learning layers. Phase 1 lives in `packages/quantum` and
establishes trusted hydrogen physics primitives used by downstream modules.

Phase 2 lives in `packages/proton_transport`. It translates the quantum hydrogen
foundation into biological proton mobility: Grotthuss hopping on water wires,
QM/MM and CP2K validation setup, constant-pH MD configuration, hydrogen-bond
geometry analysis, and pKa benchmark checks.
