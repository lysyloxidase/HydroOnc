# H2 Therapy

Phase 6 models molecular hydrogen as a research-grade anti-cancer adjuvant.

Implemented layers:

- selective ROS chemistry for hydroxyl radical and peroxynitrite
- explicit competing-scavenger caveat showing GSH dominates hydroxyl quenching
- Nrf2/NF-kB/mitochondrial signaling summaries
- CD8+ effector and exhausted T-cell ODEs in acidic tumor microenvironment
- H2 and anti-PD-1 combination/synergy simulation
- clinical response probability estimates from multi-omic features

The model uses saturated aqueous H2 near `0.78 mM`. It deliberately encodes the
central caveat from the research report: at physiological concentration, H2
contributes less than 1% of total hydroxyl radical scavenging when GSH is
included, so the clinically relevant mechanism is modeled as signaling and
immune restoration rather than bulk antioxidant stoichiometry.
