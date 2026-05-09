"""Reference constants for HydroOnc proton therapy."""

PROTON_MASS_MEV = 938.2720813
ELECTRON_MASS_MEV = 0.51099895
BETHE_K_MEV_CM2_G = 0.307075

WATER_Z_OVER_A = 0.55509
WATER_MEAN_EXCITATION_MEV = 75.0e-6
WATER_DENSITY_G_CM3 = 1.0

CLINICAL_ENERGY_RANGE_MEV = (70.0, 250.0)
CLINICAL_RANGE_TABLE_WATER_CM = {
    70.0: 4.0,
    100.0: 7.7,
    150.0: 15.6,
    200.0: 26.0,
    250.0: 37.7,
}

# Compact PSTAR-like total mass stopping power references for protons in water.
PSTAR_WATER_MASS_STOPPING_POWER = {
    70.0: 7.28,
    100.0: 6.10,
    150.0: 5.22,
    200.0: 4.90,
    250.0: 4.72,
}

FLASH_DOSE_RATE_THRESHOLD_GY_S = 40.0
CLINICAL_RBE = 1.1
