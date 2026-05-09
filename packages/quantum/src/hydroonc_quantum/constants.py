"""CODATA 2018 constants for hydrogen physics.

All primary constants are SI values. Energy conveniences with no suffix are in
electronvolts because the public HydroOnc quantum API reports spectroscopy
energies in eV; SI variants carry a ``_J`` suffix.
"""

from math import pi

# Fundamental constants, CODATA 2018 exact or recommended values.
m_e = 9.1093837015e-31  # Electron mass (kg)
m_p = 1.67262192369e-27  # Proton mass (kg)
mu_H = m_e * m_p / (m_e + m_p)  # Reduced mass of hydrogen (kg)
e = 1.602176634e-19  # Elementary charge (C)
hbar = 1.054571817e-34  # Reduced Planck constant (J s)
h = 2.0 * pi * hbar  # Planck constant (J s)
epsilon_0 = 8.8541878128e-12  # Vacuum permittivity (F/m)
k_B = 1.380649e-23  # Boltzmann constant (J/K)
c = 299792458.0  # Speed of light (m/s)

# Derived constants.
a_0 = 4.0 * pi * epsilon_0 * hbar**2 / (m_e * e**2)  # Bohr radius (m)
a_H = a_0 * (m_e / mu_H)  # Hydrogen reduced-mass Bohr radius (m)
E_h_J = e**2 / (4.0 * pi * epsilon_0 * a_0)  # Hartree energy (J)
E_h = E_h_J / e  # Hartree energy (eV)
Ry_J = E_h_J / 2.0  # Rydberg energy (J)
Ry = E_h / 2.0  # Rydberg energy (eV)
Ry_H_J = Ry_J * (mu_H / m_e)  # Hydrogen reduced-mass Rydberg energy (J)
Ry_H = Ry * (mu_H / m_e)  # Hydrogen reduced-mass Rydberg energy (eV)
R_inf = Ry_J / (h * c)  # Rydberg constant for infinite nuclear mass (1/m)
R_H = R_inf * (mu_H / m_e)  # Hydrogen Rydberg constant (1/m)
alpha = e**2 / (4.0 * pi * epsilon_0 * hbar * c)  # Fine-structure constant
mu_B = e * hbar / (2.0 * m_e)  # Bohr magneton (J/T)


def E_n(n: int) -> float:
    """Hydrogenic Bohr level in eV using the infinite-mass Rydberg."""

    if n < 1:
        raise ValueError("n must be a positive principal quantum number")
    return -Ry / n**2
