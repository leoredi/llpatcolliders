"""
Physical constants and particle data for LLP decay calculations.

All masses are in GeV, following PDG conventions.
"""

import numpy as np

# Fundamental constants
ALPHA_EM = 1.0 / 137.036
G_FERMI = 1.1663788e-5      # GeV^-2
HBAR_C_GEV_M = 1.97327e-16  # GeV·m
HBAR_GEV_S = 6.582119e-25   # GeV·s

# Meson decay constants (GeV)
F_PI = 0.1304   # Pion decay constant
F_K = 0.1556    # Kaon decay constant

# Particle masses in GeV
MASS = {
    # Leptons
    'electron': 0.000510999,
    'muon': 0.105658,
    'tau': 1.77686,

    # Neutrinos (treated as massless for decay kinematics)
    'nu_e': 0.0,
    'nu_mu': 0.0,
    'nu_tau': 0.0,

    # Mesons
    'pi_charged': 0.13957,
    'pi_neutral': 0.13498,
    'K_charged': 0.49368,
    'K_neutral': 0.49761,
    'K_short': 0.49761,
    'K_long': 0.49761,

    # Vector mesons
    'rho_charged': 0.77526,
    'rho_neutral': 0.77526,
    'omega': 0.78265,
    'eta': 0.54786,
    'eta_prime': 0.95778,

    # Gauge bosons
    'W': 80.377,
    'Z': 91.1876,
    'photon': 0.0,
}

# PDG particle IDs
PDG = {
    # Leptons (negative charge)
    'electron': 11,
    'muon': 13,
    'tau': 15,

    # Neutrinos
    'nu_e': 12,
    'nu_mu': 14,
    'nu_tau': 16,

    # Mesons
    'pi_charged': 211,    # π+
    'pi_neutral': 111,    # π0
    'K_charged': 321,     # K+
    'K_neutral': 311,     # K0
    'K_short': 310,       # K_S
    'K_long': 130,        # K_L

    # Vector mesons
    'rho_charged': 213,   # ρ+
    'rho_neutral': 113,   # ρ0
    'omega': 223,         # ω
    'eta': 221,           # η
    'eta_prime': 331,     # η'

    # Gauge bosons
    'W': 24,
    'Z': 23,
    'photon': 22,
}

# CKM matrix elements (magnitudes, PDG 2024)
V_UD = 0.97435
V_US = 0.2243
V_UB = 0.00382
V_CD = 0.221
V_CS = 0.975
V_CB = 0.0408
V_TD = 0.0086
V_TS = 0.0415
V_TB = 1.014


def get_mass(particle_name: str) -> float:
    """
    Get particle mass by name.

    Parameters
    ----------
    particle_name : str
        Name of the particle (e.g., 'electron', 'pi_charged')

    Returns
    -------
    float
        Mass in GeV

    Raises
    ------
    KeyError
        If particle not found
    """
    return MASS[particle_name]


def get_pdg(particle_name: str) -> int:
    """
    Get PDG ID by particle name.

    Parameters
    ----------
    particle_name : str
        Name of the particle

    Returns
    -------
    int
        PDG particle ID

    Raises
    ------
    KeyError
        If particle not found
    """
    return PDG[particle_name]


def is_charged(particle_name: str) -> bool:
    """
    Check if a particle is electrically charged.

    Parameters
    ----------
    particle_name : str
        Name of the particle

    Returns
    -------
    bool
        True if particle is charged
    """
    neutral_particles = {
        'nu_e', 'nu_mu', 'nu_tau',
        'pi_neutral', 'K_neutral', 'K_short', 'K_long',
        'rho_neutral', 'omega', 'eta', 'eta_prime',
        'photon', 'Z'
    }
    return particle_name not in neutral_particles
