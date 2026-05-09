"""Curated oncology H-bond fixtures used by Phase 5 deterministic gates."""

from hydroonc_hbond.geometry import HBondInteraction


P53_WT_2XWR_HBONDS = (
    HBondInteraction("R175_bb_L2_1", "R175:N", "D184:O", 2.91, 158.0, 0.82, "backbone", True),
    HBondInteraction("R175_bb_L2_2", "R175:NH1", "P177:O", 2.88, 151.0, 0.79, "backbone", True),
    HBondInteraction("R175_bb_L3_1", "R175:NH2", "C176:O", 2.95, 149.0, 0.76, "backbone", True),
    HBondInteraction("R175_bb_L3_2", "R175:N", "H179:O", 3.02, 144.0, 0.68, "backbone", True),
    HBondInteraction("R175_bb_zinc_1", "R175:NH1", "C238:O", 2.84, 153.0, 0.74, "backbone", True),
    HBondInteraction("R175_bb_zinc_2", "R175:NH2", "C242:O", 2.93, 148.0, 0.70, "backbone", True),
    HBondInteraction("R175_bb_core_1", "R175:N", "G245:O", 3.10, 139.0, 0.61, "backbone", False),
    HBondInteraction("R175_bb_core_2", "R175:NH1", "R249:O", 3.05, 142.0, 0.58, "backbone", False),
    HBondInteraction("R175_sc_1", "R175:NH1", "D184:OD1", 2.78, 164.0, 0.86, "sidechain", True),
    HBondInteraction("R175_sc_2", "R175:NH2", "D184:OD2", 2.82, 161.0, 0.83, "sidechain", True),
    HBondInteraction("R175_sc_3", "R175:NE", "C176:O", 2.96, 150.0, 0.66, "sidechain", True),
    HBondInteraction("R175_sc_4", "R175:NH1", "H179:ND1", 3.08, 136.0, 0.54, "sidechain", False),
    HBondInteraction("R273_DNA_phosphate", "R273:NH1", "DNA:OP1", 2.72, 166.0, 0.88, "dna_contact", True),
    HBondInteraction("R248_minor_groove", "R248:NH2", "DNA:N3", 2.86, 159.0, 0.81, "dna_contact", True),
    HBondInteraction("G245_turn_water", "G245:N", "HOH451:O", 2.95, 145.0, 0.52, "water", False, True),
    HBondInteraction("R282_H2_helix", "R282:NH1", "E286:OE1", 2.89, 157.0, 0.77, "sidechain", True),
)


P53_MUTANT_HBONDS = {
    "R175H": (
        HBondInteraction("R273_DNA_phosphate", "R273:NH1", "DNA:OP1", 2.73, 164.0, 0.72, "dna_contact", True),
        HBondInteraction("R248_minor_groove", "R248:NH2", "DNA:N3", 2.92, 154.0, 0.64, "dna_contact", True),
        HBondInteraction("H175_new_weak", "H175:NE2", "D184:OD1", 3.21, 124.0, 0.18, "sidechain", False),
    ),
    "R273H": tuple(
        interaction
        if interaction.identifier != "R273_DNA_phosphate"
        else HBondInteraction("R273_DNA_phosphate", "H273:NE2", "DNA:OP1", 3.68, 101.0, 0.06, "dna_contact", True)
        for interaction in P53_WT_2XWR_HBONDS
    ),
    "R248W": tuple(
        interaction
        if interaction.identifier != "R248_minor_groove"
        else HBondInteraction("R248_minor_groove", "W248:NE1", "DNA:N3", 3.91, 88.0, 0.03, "dna_contact", False)
        for interaction in P53_WT_2XWR_HBONDS
    ),
    "G245S": tuple(P53_WT_2XWR_HBONDS)
    + (HBondInteraction("S245_inappropriate", "S245:OG", "R249:NH1", 2.79, 151.0, 0.71, "mispaired", False),),
    "R249S": tuple(
        interaction
        if interaction.identifier != "R175_bb_core_2"
        else HBondInteraction("R175_bb_core_2", "R175:NH1", "S249:OG", 3.52, 116.0, 0.12, "backbone", False)
        for interaction in P53_WT_2XWR_HBONDS
    ),
    "R282W": tuple(
        interaction
        if interaction.identifier != "R282_H2_helix"
        else HBondInteraction("R282_H2_helix", "W282:NE1", "E286:OE1", 3.74, 104.0, 0.08, "sidechain", True)
        for interaction in P53_WT_2XWR_HBONDS
    ),
}


MANUAL_COUNTS = {
    "p53_wt_2xwr_core": 16,
    "p53_R175H": 2,
    "kras_g12c_sotorasib": 1,
    "egfr_wt_gefitinib": 1,
    "egfr_T790M_gefitinib": 0,
}
