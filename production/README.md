# HNL Production Simulation for Far-Detector Studies

## Overview

This package provides publication-quality Heavy Neutral Lepton (HNL) production 
simulation in PYTHIA8, following the methodology established by:
- **MATHUSLA** (David Curtin et al., GitHub: MATHUSLA_LLPfiles_RHN_*)
- **ANUBIS** (arXiv:1909.13022, arXiv:2101.xxxxx)
- **Physics Beyond Colliders** benchmarks (arXiv:1901.09966)
- **HNL phenomenology review** (arXiv:1805.08567, Bondarenko et al.)

## Key Physics Principles

### What This Package Does (Stage 1: Production Kinematics)
- Generates HNL 4-vectors from various parent particles
- Records production vertex, parent ID, boost factor
- Keeps HNL **stable** in Pythia (no decay simulated)
- Outputs weighted events for later analysis

### What Stage 2 Handles (NOT done here)
- HNL lifetime: cτ₀(m_N, |U|²) from HNLCalc
- Decay probability: P_decay = exp(-d₁/λ) × (1 - exp(-L/λ))
- Where λ = βγ × cτ₀ (proper decay length × boost)
- Specific decay channels and kinematics

## Physics Beyond Colliders Benchmarks

We follow the PBC benchmark scenarios (arXiv:1901.09966):

| Benchmark | Coupling | Description |
|-----------|----------|-------------|
| **BC6** | |U_e|² | Single HNL, electron dominance |
| **BC7** | |U_μ|² | Single HNL, muon dominance |
| **BC8** | |U_τ|² | Single HNL, tau dominance |

Each benchmark assumes one Majorana HNL with mixing to a single SM neutrino flavor.

## Production Modes by Mass Regime

### Low Mass (m_N < 0.5 GeV): Kaon-dominated
- **K⁺ → ℓ⁺ N** (2-body leptonic)
- **K_L → π ℓ N** (semileptonic, subdominant)
- Kinematic limit: m_N < m_K - m_ℓ ≈ 0.39 GeV (for muon)

### Intermediate Mass (0.5 < m_N < 2 GeV): Charm-dominated  
- **D⁺ → ℓ⁺ N** (2-body leptonic)
- **D_s⁺ → ℓ⁺ N** (2-body leptonic)
- **D → K ℓ N** (semileptonic 3-body)
- Kinematic limits: m_N < m_D - m_ℓ ≈ 1.76 GeV

### High Mass (2 < m_N < 5 GeV): Beauty-dominated
- **B⁺ → ℓ⁺ N** (2-body leptonic)
- **B → D ℓ N** (semileptonic 3-body)
- **B → D* ℓ N** (semileptonic 3-body)
- **B_c⁺ → ℓ⁺ N** (2-body, small cross section)
- Kinematic limits: m_N < m_B - m_ℓ ≈ 5.17 GeV

### Electroweak Regime (m_N > 5 GeV)
- **W → ℓ N** (Drell-Yan production)
- **Z → ν N** (neutral current)
- **t → b W → b ℓ N** (top decays)
- Kinematic limit: m_N < m_W ≈ 80 GeV

## Methodology Comparison

### Option A: FONLL + External Decay (Gold Standard)
Used by MATHUSLA, most accurate for meson production:
1. Use FONLL for B/D meson dσ/dpT/dy (NLO+NNLL)
2. Sample meson 4-vectors from FONLL distributions
3. Decay externally using proper form factors
4. **Pros**: Most accurate meson kinematics
5. **Cons**: Requires external tools, complex

### Option B: MadGraph + Pythia (For EW production)
Used for W/Z production:
1. Use MadGraph5 with SM_HeavyN_CKM_AllMasses_LO model
2. Generate `p p > w+, w+ > e+ n1` at parton level
3. Shower with Pythia8
4. **Pros**: Correct matrix elements for EW
5. **Cons**: Requires MadGraph installation

### Option C: Pythia-Only (This Package)
Validated approach (arXiv:2103.11494):
1. Generate inclusive heavy flavor or EW processes in Pythia
2. Let Pythia hadronize → realistic meson kinematics
3. Force parent → ℓ N decays
4. Extract HNL 4-vectors
5. **Pros**: Self-contained, fast, good for sensitivity estimates
6. **Cons**: Less accurate than FONLL for meson spectra

**This package uses Option C** for simplicity and self-containment.

## Physical vs Unphysical Decay Channels

### ✅ PHYSICAL channels (used here)
```
Charged pseudoscalar mesons (2-body leptonic):
  K⁺ → ℓ⁺ N    (PDG: 321)
  D⁺ → ℓ⁺ N    (PDG: 411)  
  D_s⁺ → ℓ⁺ N  (PDG: 431)
  B⁺ → ℓ⁺ N    (PDG: 521)
  B_c⁺ → ℓ⁺ N  (PDG: 541)

Semileptonic (3-body):
  D⁰ → K⁻ ℓ⁺ N    (PDG: 421)
  B⁰ → D⁻ ℓ⁺ N    (PDG: 511)
  B_s⁰ → D_s⁻ ℓ⁺ N (PDG: 531)

Electroweak (2-body):
  W⁺ → ℓ⁺ N       (PDG: 24)
  Z → ν N         (PDG: 23)

Tau decays (for BC8):
  τ⁻ → π⁻ N       (2-body)
  τ⁻ → ρ⁻ N       (2-body) 
  τ⁻ → ℓ⁻ ν̄ N    (3-body)
```

### ❌ UNPHYSICAL channels (DO NOT USE)
```
K⁰_L → ν N     ← WRONG: No 2-body leptonic decay exists
K⁰_S → ν N     ← WRONG: Same reason
D⁰ → ν N       ← WRONG: Neutral D has no 2-body leptonic
B⁰ → ν N       ← WRONG: Neutral B has no 2-body leptonic
```

## Meson Fractions (from Pythia8 hadronization)

Following MATHUSLA methodology:

**B mesons** (from bb̄ hadronization):
- B⁰ / B̄⁰: 44.8%
- B⁺ / B⁻: 44.8%  
- B_s⁰ / B̄_s⁰: 10.3%
- B_c⁺ / B_c⁻: 0.018%

**D mesons** (from cc̄ hadronization):
- D⁰ / D̄⁰: ~60%
- D⁺ / D⁻: ~25%
- D_s⁺ / D_s⁻: ~15%

## Cross Sections at 14 TeV (for normalization)

| Process | σ (pb) | Source |
|---------|--------|--------|
| pp → bb̄ X | ~500,000 | FONLL |
| pp → cc̄ X | ~10,000,000 | FONLL |
| pp → W X | ~200,000 | NLO |
| pp → Z X | ~60,000 | NLO |
| pp → tt̄ X | ~900 | NNLO |

## Output Format

CSV file with columns:
```
event, weight, hnl_id, parent_pdg,
pt, eta, phi, p, E, mass,
prod_x_mm, prod_y_mm, prod_z_mm,
boost_gamma
```

Where:
- `event`: Event number
- `weight`: Relative MC weight (pythia.info.weight())
- `hnl_id`: HNL PDG code (±9900015)
- `parent_pdg`: Parent particle PDG code (521=B⁺, 511=B⁰, etc.)
- `pt, eta, phi`: Transverse momentum [GeV], pseudorapidity, azimuthal angle
- `p, E, mass`: 3-momentum magnitude [GeV], energy [GeV], mass [GeV]
- `prod_x_mm, prod_y_mm, prod_z_mm`: Production vertex in mm
- `boost_gamma`: Lorentz boost factor γ = E/m

## Usage

```bash
# Compile
make

# Run single mass point (100k events default)
./main_hnl_production 2.0 muon

# Override number of events
./main_hnl_production 10.0 electron 500000

# More examples:
./main_hnl_production 0.3 muon 50000    # 300 MeV muon-coupled (Kaon regime)
./main_hnl_production 1.5 electron      # 1.5 GeV electron-coupled (Charm regime)
./main_hnl_production 3.0 muon          # 3.0 GeV muon-coupled (Beauty regime)
./main_hnl_production 15.0 tau 200000   # 15 GeV tau-coupled (EW regime)

# Run all benchmarks
./run_benchmarks.sh
```

**Note:** The program automatically chooses the production regime (Kaon / Charm / Beauty / EW) based on the HNL mass:
- m < 0.5 GeV: Kaon regime (SoftQCD)
- 0.5 ≤ m < 2 GeV: Charm regime (cc̄ production)
- 2 ≤ m < 5 GeV: Beauty regime (bb̄ production)
- m ≥ 5 GeV: Electroweak regime (W/Z production)

## References

1. **HNL Phenomenology**: Bondarenko et al., JHEP 11 (2018) 032, arXiv:1805.08567
2. **PBC Benchmarks**: Beacham et al., J.Phys.G 47 (2020) 010501, arXiv:1901.09966  
3. **Pythia Validation**: Gorkavenko et al., J.Phys.G 48 (2021) 105001, arXiv:2103.11494
4. **MATHUSLA**: Curtin et al., arXiv:1806.07396
5. **ANUBIS**: Bauer et al., arXiv:1909.13022
6. **Heavy Neutral Leptons at ANUBIS**: Hirsch et al., PRD 101 (2020) 055034

## Author
Generated for transverse far-detector studies, following PBC methodology.
