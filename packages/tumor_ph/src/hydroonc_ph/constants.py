"""Reference values for tumor microenvironment pH modeling."""

NORMAL_EXTRACELLULAR_PH = 7.4
NORMAL_EXTRACELLULAR_PH_RANGE = (7.2, 7.5)
TUMOR_EXTRACELLULAR_PH = 6.7
TUMOR_INTRACELLULAR_PH_RANGE = (7.0, 7.2)
NORMAL_INTRACELLULAR_PH_RANGE = (6.8, 7.2)

HCC_EXTRACELLULAR_PH_MEAN = 6.66
HCC_EXTRACELLULAR_PH_SD = 0.19
BREAST_EXTRACELLULAR_PH_RANGE = (6.3, 6.9)
MELANOMA_EXTRACELLULAR_PH = 6.96

D_H_TISSUE_CM2_S = 9.0e-5
D_LACTATE_TISSUE_CM2_S = 5.0e-6
D_O2_TISSUE_CM2_S = 2.0e-5

DEFAULT_REGULATORS = {
    "v_atpase": 1.0,
    "nhe1": 1.0,
    "mct1": 1.0,
    "mct4": 1.0,
    "caix": 1.0,
    "nbce1": 1.0,
    "nbcn1": 1.0,
    "ae2": 1.0,
}

INHIBITOR_IC50_UM = {
    "v_atpase": 0.01,
    "nhe1": 1.0,
    "mct1": 0.10,
    "mct4": 2.0,
    "caix": 0.25,
}

INHIBITOR_ALIASES = {
    "bafilomycin": "v_atpase",
    "bafilomycin a1": "v_atpase",
    "v-atpase": "v_atpase",
    "vatpase": "v_atpase",
    "cariporide": "nhe1",
    "eipa": "nhe1",
    "slc9a1": "nhe1",
    "nhe1": "nhe1",
    "azd3965": "mct1",
    "mct1": "mct1",
    "mct4": "mct4",
    "slc-0111": "caix",
    "slc0111": "caix",
    "caix": "caix",
    "ca9": "caix",
}
