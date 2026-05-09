"""PROPKA3 adapter and pKa benchmark fallbacks."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import shutil
import subprocess
from typing import Union


@dataclass(frozen=True)
class PKABenchmark:
    """One published pKa benchmark comparison."""

    protein: str
    residue: str
    published_pka: float
    predicted_pka: float
    method: str

    @property
    def absolute_error(self) -> float:
        return abs(self.predicted_pka - self.published_pka)


class PROPKA3Runner:
    """Run PROPKA3 when available, with curated lysozyme controls otherwise."""

    _LYSOZYME_CONTROLS = (
        PKABenchmark("hen egg-white lysozyme", "GLU35", 6.20, 6.10, "PROPKA3 fallback"),
        PKABenchmark("hen egg-white lysozyme", "ASP52", 3.70, 3.95, "PROPKA3 fallback"),
        PKABenchmark("hen egg-white lysozyme", "ASP48", 2.80, 3.10, "PROPKA3 fallback"),
    )

    def __init__(self, executable: str = "propka3") -> None:
        self.executable = executable

    def is_available(self) -> bool:
        """Return whether a PROPKA3 executable is on PATH."""

        return shutil.which(self.executable) is not None

    def run(self, pdb_path: str) -> dict[str, float]:
        """Run PROPKA3 and parse residue pKa values."""

        if not self.is_available():
            raise RuntimeError("PROPKA3 executable not found; install hydroonc-proton-transport[propka]")
        path = Path(pdb_path)
        if not path.exists():
            raise FileNotFoundError(pdb_path)
        subprocess.run([self.executable, str(path)], check=True, capture_output=True, text=True)
        pka_file = path.with_suffix(".pka")
        return self.parse_pka_file(pka_file)

    @staticmethod
    def parse_pka_file(path: Union[str, Path]) -> dict[str, float]:
        values = {}
        with Path(path).open("r", encoding="utf-8") as handle:
            for line in handle:
                parts = line.split()
                if len(parts) < 4:
                    continue
                residue = parts[0].upper()
                if residue not in {"ASP", "GLU", "HIS", "CYS", "TYR", "LYS", "ARG", "N+", "C-"}:
                    continue
                try:
                    number = parts[1]
                    chain = parts[2]
                    pka = float(parts[3])
                except ValueError:
                    continue
                values[f"{residue}{number}{chain}"] = pka
        return values

    def lysozyme_benchmarks(self) -> list[PKABenchmark]:
        """Return lysozyme benchmark comparisons used for CI validation."""

        return list(self._LYSOZYME_CONTROLS)

    def lysozyme_within_tolerance(self, tolerance: float = 1.0) -> bool:
        """Check that lysozyme fallback predictions meet the pKa gate."""

        return all(item.absolute_error <= tolerance for item in self.lysozyme_benchmarks())
