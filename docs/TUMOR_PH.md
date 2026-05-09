# Tumor pH

Phase 3 models tumor microenvironment acidity across three layers:

- Warburg-effect acid production and regulator-driven pH dynamics
- reaction-diffusion gradients from acidic tumor core to vascular boundary
- CEST-MRI pH map reconstruction with a lightweight ML architecture

Reference values used by the package:

- HCC extracellular pH: `6.66 +/- 0.19`
- breast tumor extracellular pH: `6.3-6.9`
- melanoma extracellular pH: about `6.96`
- normal extracellular pH: `7.2-7.5`
- tumor intracellular pH: `7.0-7.2`

The model includes V-ATPase, NHE1, MCT1/MCT4, CAIX/XII, NBCe1/NBCn1, and AE2.
Inhibitor response helpers cover bafilomycin A1, cariporide/EIPA, AZD3965, and
SLC-0111.
