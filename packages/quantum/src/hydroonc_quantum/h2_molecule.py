"""Electronic structure helpers for molecular hydrogen."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Optional

import numpy as np

from hydroonc_quantum.constants import E_h


H2_REFERENCE = {
    "R_e_angstrom": 0.7414,
    "D_e_eV": 4.747,
    "hf_cc-pVTZ_Eh": -1.13286,
    "ccsd_cc-pVTZ_Eh": -1.17234,
    "fci_cc-pV5Z_Eh": -1.17447,
}


@dataclass(frozen=True)
class PotentialEnergyCurve:
    """H2 potential energy curve."""

    R_angstrom: np.ndarray
    energy_Eh: np.ndarray
    method: str
    basis: str
    backend: str

    @property
    def equilibrium_bond_length(self) -> float:
        return float(self.R_angstrom[int(np.argmin(self.energy_Eh))])

    @property
    def minimum_energy(self) -> float:
        return float(np.min(self.energy_Eh))


class H2Molecule:
    """Electronic structure of H2 via PySCF, with transparent reference fallback."""

    def __init__(self, *, reference_fallback: bool = True, verbose: int = 0) -> None:
        self.reference_fallback = reference_fallback
        self.verbose = verbose

    def compute_energy(
        self,
        R_angstrom: float = 0.7414,
        method: str = "ccsd",
        basis: str = "cc-pVTZ",
    ) -> dict:
        """Compute H2 electronic energy at a bond length.

        If PySCF is unavailable and ``reference_fallback`` is enabled, a Morse
        reference curve calibrated to accepted H2 constants is returned with
        ``backend='reference'``.
        """

        method_key = method.lower()
        try:
            from pyscf import cc, dft, gto, scf
        except ModuleNotFoundError:
            if not self.reference_fallback:
                raise
            return self._reference_energy(R_angstrom, method_key, basis)

        mol = gto.M(
            atom=f"H 0 0 0; H 0 0 {R_angstrom}",
            basis=basis,
            unit="Angstrom",
            verbose=self.verbose,
        )
        if method_key in {"hf", "rhf"}:
            mf = scf.RHF(mol).run(verbose=self.verbose)
            return {
                "energy_Eh": float(mf.e_tot),
                "method": "RHF",
                "basis": basis,
                "R_angstrom": R_angstrom,
                "backend": "pyscf",
            }
        if method_key in {"b3lyp", "dft"}:
            mf = dft.RKS(mol, xc="b3lyp").run(verbose=self.verbose)
            return {
                "energy_Eh": float(mf.e_tot),
                "method": "B3LYP",
                "basis": basis,
                "R_angstrom": R_angstrom,
                "backend": "pyscf",
            }
        if method_key in {"ccsd", "ccsd(t)", "ccsdt"}:
            mf = scf.RHF(mol).run(verbose=self.verbose)
            mycc = cc.CCSD(mf).run(verbose=self.verbose)
            energy = float(mf.e_tot + mycc.e_corr)
            label = "CCSD"
            if method_key in {"ccsd(t)", "ccsdt"}:
                energy += float(mycc.ccsd_t())
                label = "CCSD(T)"
            return {
                "energy_Eh": energy,
                "method": label,
                "basis": basis,
                "R_angstrom": R_angstrom,
                "backend": "pyscf",
            }
        raise ValueError("method must be one of hf, b3lyp, ccsd, or ccsd(t)")

    def potential_energy_curve(
        self,
        R_range: tuple[float, float] = (0.5, 5.0),
        n_points: int = 50,
        method: str = "ccsd",
        basis: str = "cc-pVTZ",
        R_values: Optional[Iterable[float]] = None,
    ) -> PotentialEnergyCurve:
        """Scan E(R) for H2."""

        if R_values is None:
            R = np.linspace(R_range[0], R_range[1], n_points)
        else:
            R = np.asarray(list(R_values), dtype=float)
        energies = []
        backends = set()
        for value in R:
            result = self.compute_energy(float(value), method=method, basis=basis)
            energies.append(result["energy_Eh"])
            backends.add(result["backend"])
        backend = "mixed" if len(backends) > 1 else next(iter(backends))
        return PotentialEnergyCurve(
            R_angstrom=R,
            energy_Eh=np.asarray(energies, dtype=float),
            method=method.upper(),
            basis=basis,
            backend=backend,
        )

    @staticmethod
    def molecular_orbitals() -> dict:
        """Qualitative H2 molecular-orbital occupancy."""

        return {
            "bonding": "sigma_g(1s)^2",
            "antibonding": "sigma_u*(1s)^0",
            "bond_order": 1.0,
        }

    @staticmethod
    def _reference_energy(R_angstrom: float, method: str, basis: str) -> dict:
        if method in {"hf", "rhf"}:
            minimum = H2_REFERENCE["hf_cc-pVTZ_Eh"]
            label = "RHF"
        elif method in {"ccsd", "ccsd(t)", "ccsdt"}:
            minimum = H2_REFERENCE["ccsd_cc-pVTZ_Eh"]
            label = "CCSD" if method == "ccsd" else "CCSD(T)"
        elif method in {"b3lyp", "dft"}:
            minimum = -1.1740
            label = "B3LYP"
        else:
            raise ValueError("method must be one of hf, b3lyp, ccsd, or ccsd(t)")

        R_e = H2_REFERENCE["R_e_angstrom"]
        D_e_Eh = H2_REFERENCE["D_e_eV"] / E_h
        morse_a = 1.942
        energy = minimum + D_e_Eh * (1.0 - np.exp(-morse_a * (R_angstrom - R_e))) ** 2
        return {
            "energy_Eh": float(energy),
            "method": label,
            "basis": basis,
            "R_angstrom": R_angstrom,
            "backend": "reference",
            "note": "Install hydroonc-quantum[h2] for PySCF calculations.",
        }
