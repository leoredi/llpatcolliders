"""
Core functionality for LLP decay physics.

This module provides:
- Physical constants and particle data
- Relativistic kinematics (Lorentz boosts, phase space)
- Random sampling utilities
"""

from .particle import (
    MASS, PDG,
    ALPHA_EM, G_FERMI, HBAR_C_GEV_M, HBAR_GEV_S,
    F_PI, F_K,
    V_UD, V_US, V_UB, V_CD, V_CS, V_CB, V_TD, V_TS, V_TB,
    get_mass, get_pdg, is_charged
)

from .kinematics import (
    kallen,
    four_vector_mass,
    boost_to_lab,
    two_body_decay_momenta,
    sample_two_body_decay,
    invariant_mass_from_daughters,
    transverse_momentum,
    rapidity,
    pseudorapidity
)

from .sampling import (
    create_rng,
    sample_discrete,
    rejection_sampling,
    weighted_sample,
    uniform_sphere
)

__all__ = [
    # Constants
    'MASS', 'PDG',
    'ALPHA_EM', 'G_FERMI', 'HBAR_C_GEV_M', 'HBAR_GEV_S',
    'F_PI', 'F_K',
    'V_UD', 'V_US', 'V_UB', 'V_CD', 'V_CS', 'V_CB', 'V_TD', 'V_TS', 'V_TB',

    # Particle functions
    'get_mass', 'get_pdg', 'is_charged',

    # Kinematics
    'kallen',
    'four_vector_mass',
    'boost_to_lab',
    'two_body_decay_momenta',
    'sample_two_body_decay',
    'invariant_mass_from_daughters',
    'transverse_momentum',
    'rapidity',
    'pseudorapidity',

    # Sampling
    'create_rng',
    'sample_discrete',
    'rejection_sampling',
    'weighted_sample',
    'uniform_sphere',
]
