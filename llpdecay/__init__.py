"""
LLPDecay: Long-Lived Particle Decay Generator

A Python package for generating LLP decay products with proper kinematics
and matrix elements. Designed to plug into existing HNL sensitivity pipelines
that use MadGraph+Pythia for production.

Key Features
------------
- Production-independent decay sampling
- Proper Lorentz boosts and phase space
- Matrix element weighted branching ratios
- Support for polarized decays
- Modular architecture for multiple LLP types

Quick Start
-----------
>>> from llpdecay import HNL
>>> import numpy as np

>>> # Create HNL model
>>> hnl = HNL(mass=2.0, Umu=1e-6)  # 2 GeV, muon mixing

>>> # Sample decay
>>> parent_4vec = np.array([10.0, 3.0, 0.5, 9.3])  # [E, px, py, pz] in GeV
>>> daughters, channel = hnl.sample_decay(parent_4vec, return_channel=True)
>>> print(f"Decayed via {channel}")
>>> print(f"Daughter 4-vectors:\\n{daughters}")

>>> # Get branching ratios
>>> brs = hnl.branching_ratios()
>>> for ch, br in sorted(brs.items(), key=lambda x: -x[1])[:5]:
...     print(f"{ch}: {br:.1%}")

Package Structure
-----------------
llpdecay/
├── core/          # Physical constants, kinematics, sampling utilities
├── models/        # LLP model implementations (HNL, ALP, ...)
├── decays/        # Decay channel definitions
└── tests/         # Validation tests

References
----------
- Bondarenko et al., JHEP 1901 (2019) 127, arXiv:1805.08567
- Gorbunov & Shaposhnikov, Phys.Rev.D75:083010, arXiv:hep-ph/9911364
- Atre et al., JHEP 0905:030, arXiv:0901.3589
"""

__version__ = '0.1.0'
__author__ = 'LLPatColliders Team'

# Core functionality
from .core import (
    # Constants
    MASS, PDG,
    ALPHA_EM, G_FERMI, HBAR_C_GEV_M, HBAR_GEV_S,
    F_PI, F_K,
    get_mass, get_pdg, is_charged,

    # Kinematics
    kallen,
    four_vector_mass,
    boost_to_lab,
    two_body_decay_momenta,
    sample_two_body_decay,
    invariant_mass_from_daughters,
    transverse_momentum,
    rapidity,
    pseudorapidity,

    # Sampling
    create_rng,
    sample_discrete,
    rejection_sampling,
    weighted_sample,
    uniform_sphere,
)

# Decay channels
from .decays import (
    DecayChannel,
    HNL_CHANNELS,
    get_channel,
    list_channels,
    group_by_flavor,
)

# Models
from .models import (
    LLPModel,
    HNL,
    ALP,
)

# Three-body decays (Phase 2)
from .decays import (
    ThreeBodyPhaseSpace,
    sample_three_body_decay,
    hnl_three_body_leptonic_me,
)

__all__ = [
    # Version
    '__version__',
    '__author__',

    # Constants
    'MASS', 'PDG',
    'ALPHA_EM', 'G_FERMI', 'HBAR_C_GEV_M', 'HBAR_GEV_S',
    'F_PI', 'F_K',
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

    # Decay channels
    'DecayChannel',
    'HNL_CHANNELS',
    'get_channel',
    'list_channels',
    'group_by_flavor',

    # Three-body decays
    'ThreeBodyPhaseSpace',
    'sample_three_body_decay',
    'hnl_three_body_leptonic_me',

    # Models
    'LLPModel',
    'HNL',
    'ALP',
]
