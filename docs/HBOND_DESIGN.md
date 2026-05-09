# Hydrogen-Bond Design

Phase 5 focuses on hydrogen-bond networks in mutant cancer proteins:

- p53 hotspot mutations: R175H, G245S, R248W, R249S, R273H, and R282W
- KRAS G12C switch-II pocket H-bonds used by sotorasib and adagrasib
- EGFR T790M/C797S gatekeeper and covalent-resistance logic
- H-bond-aware protein-ligand graph construction
- lightweight GNN-style binding-affinity prediction

The p53 R175H gate models loss of the L2/L3 zinc-site H-bond network, while
R273H models loss of direct DNA phosphate contact occupancy. These curated
fixtures are deterministic CI controls for later MD/FoldX/AlphaFold-derived
structures.
