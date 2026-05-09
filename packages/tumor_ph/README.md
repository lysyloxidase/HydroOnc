# HydroOnc Tumor pH

Phase 3 models tumor acidosis and pH mapping:

- Warburg-effect acid production and regulator-driven pH dynamics
- inhibitor response estimates for V-ATPase, NHE1, MCT1/MCT4, and CAIX
- reaction-diffusion pH gradients with finite-difference fallback for FEniCS
- PINN-style surrogate solutions for pH fields
- CEST-MRI pH mapping with a lightweight NumPy architecture

Install for development:

```bash
pip install -e "packages/tumor_ph[dev]"
```

Optional numerical backends:

```bash
pip install -e "packages/tumor_ph[dev,fenics,ml]"
```
