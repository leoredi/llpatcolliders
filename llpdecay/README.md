# LLPDecay: Long-Lived Particle Decay Generator

A Python package for generating LLP decay products with proper kinematics and matrix elements.

## Overview

**LLPDecay** plugs into existing HNL sensitivity pipelines that use MadGraph+Pythia for production. While production tools give us LLP 4-vectors, this package handles:

- ✅ Proper decay kinematics (Lorentz boosts, phase space)
- ✅ Matrix element weighted branching ratios
- ✅ Support for polarized decays
- ✅ Modular architecture for multiple LLP types (HNL, ALP, ...)

## Installation

```bash
# From the repository root
cd llpatcolliders
python -m pip install -e llpdecay/
```

## Quick Start

```python
from llpdecay import HNL
import numpy as np

# Create HNL model
hnl = HNL(mass=2.0, Umu=1e-6)  # 2 GeV mass, muon mixing

# Get branching ratios
brs = hnl.branching_ratios()
for channel, br in sorted(brs.items(), key=lambda x: -x[1])[:5]:
    print(f"{channel}: {br:.1%}")

# Sample decay
parent_4vec = np.array([10.0, 3.0, 0.5, 9.3])  # [E, px, py, pz] in GeV
daughters, channel = hnl.sample_decay(parent_4vec, return_channel=True)

print(f"\nDecayed via: {channel}")
print(f"Daughter 1: {daughters[0, 0]}")  # [E, px, py, pz]
print(f"Daughter 2: {daughters[0, 1]}")

# Get decay properties
print(f"\nLifetime: {hnl.lifetime():.3e} s")
print(f"Decay length: {hnl.ctau():.3e} m")
```

## Features

### Implemented (Phase 1)
- ✅ Core kinematics (Lorentz boosts, 2-body phase space)
- ✅ HNL model with proper matrix elements
- ✅ Two-body decay channels (ℓπ, ℓK, ℓρ, νπ⁰)
- ✅ Branching ratio calculations
- ✅ Polarization effects
- ✅ Majorana vs Dirac modes
- ✅ Lifetime and decay probability calculations

### Planned (Phase 2)
- ⏳ Three-body decay channels (ν ℓ ℓ)
- ⏳ Integration with HNLCalc for BR validation
- ⏳ ALP model implementation
- ⏳ Radiative corrections

## Physics References

- Bondarenko et al., [arXiv:1805.08567](https://arxiv.org/abs/1805.08567) - HNL phenomenology
- Gorbunov & Shaposhnikov, [arXiv:hep-ph/9911364](https://arxiv.org/abs/hep-ph/9911364) - Sterile neutrino decays
- Atre et al., [arXiv:0901.3589](https://arxiv.org/abs/0901.3589) - HNL production and decay

## Package Structure

```
llpdecay/
├── core/              # Physical constants, kinematics, sampling
│   ├── particle.py    # Particle data (masses, PDG IDs)
│   ├── kinematics.py  # Lorentz boosts, phase space
│   └── sampling.py    # Random sampling utilities
├── models/            # LLP model implementations
│   ├── base.py        # Abstract base class
│   ├── hnl.py         # Heavy Neutral Lepton
│   └── alp.py         # Axion-Like Particle (Phase 2)
├── decays/            # Decay channel definitions
│   └── channels.py    # Channel database
└── tests/             # Validation tests
    ├── test_kinematics.py
    ├── test_hnl_decays.py
    └── test_validation.py
```

## Running Tests

```bash
cd llpdecay
pytest tests/ -v
```

## Integration with Existing Pipeline

The package is designed to integrate with `analysis_pbc/limits/expected_signal.py`:

```python
from llpdecay import HNL

# In your signal calculation
hnl_model = HNL(mass=m_N, Ue=Ue, Umu=Umu, Utau=Utau)

for event in hnl_events:
    # event has: E, px, py, pz from MadGraph+Pythia
    parent_4vec = np.array([event.E, event.px, event.py, event.pz])

    # Sample decay
    daughters, channel = hnl_model.sample_decay(
        parent_4vec,
        polarization=event.polarization,
        return_channel=True
    )

    # Apply reconstruction cuts
    n_charged = hnl_model.get_charged_count(channel)
    if n_charged >= 2:
        # Apply pT cuts, geometric acceptance, etc.
        ...
```

## API Documentation

### HNL Model

```python
HNL(mass, Ue=0, Umu=0, Utau=0, is_majorana=True)
```

**Parameters:**
- `mass` (float): HNL mass in GeV
- `Ue`, `Umu`, `Utau` (float): Mixing angles squared |U_α|²
- `is_majorana` (bool): Majorana (True) vs Dirac (False)

**Methods:**
- `available_channels()`: List kinematically accessible channels
- `branching_ratios()`: Get {channel: BR} dictionary
- `sample_decay(parent_4vec, ...)`: Generate decay products
- `total_width()`: Total decay width (GeV)
- `lifetime()`: Mean lifetime (seconds)
- `ctau()`: Mean decay length (meters)
- `decay_probability(parent_4vec, x_min, x_max)`: P(decay in interval)

### Decay Channels

Supported HNL decay modes:

**Two-body charged current:**
- `e_pi`, `mu_pi`, `tau_pi`: N → ℓ⁻ π⁺
- `e_K`, `mu_K`, `tau_K`: N → ℓ⁻ K⁺
- `e_rho`, `mu_rho`, `tau_rho`: N → ℓ⁻ ρ⁺

**Two-body neutral current:**
- `nu_e_pi0`, `nu_mu_pi0`, `nu_tau_pi0`: N → ν π⁰ (invisible)

**Three-body (Phase 2):**
- `nu_e_e`, `nu_mu_mu`, `nu_tau_tau`: N → ν ℓ⁺ ℓ⁻
- `nu_e_mu`, `nu_mu_e`: N → ν ℓ₁⁺ ℓ₂⁻

## Contributing

Please add new features in separate modules following the existing structure. All physics formulas should include literature references.

## License

Part of the LLPatColliders project.
