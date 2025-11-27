# Physics Documentation: HNL Production for Far-Detector Studies

## 1. Introduction

Heavy Neutral Leptons (HNLs), also known as sterile neutrinos or right-handed neutrinos, 
are hypothetical particles that appear in many extensions of the Standard Model (SM). 
They are prime targets for far-detector experiments at the LHC such as MATHUSLA, ANUBIS, 
CODEX-b, and FASER.

This document explains the physics behind our HNL production simulation and how it 
interfaces with downstream analysis (HNLCalc, geometry code).

## 2. HNL Basics

### 2.1 The Neutrino Portal

HNLs interact with the SM through the **neutrino portal**:

```
L = F_α (L̄_α H̃) N + h.c.
```

where `L_α` is the SM lepton doublet, `H` is the Higgs doublet, and `F_α` is the 
Yukawa coupling. After electroweak symmetry breaking, this generates active-sterile 
mixing:

```
|U_α|² = |F_α v / M_N|²
```

where `v = 246 GeV` is the Higgs VEV and `M_N` is the HNL mass.

### 2.2 PBC Benchmark Points

The Physics Beyond Colliders (PBC) initiative defines simplified benchmarks:

| Benchmark | Dominant Mixing | Parameter Space |
|-----------|-----------------|-----------------|
| BC6 | |U_e|² | (m_N, |U_e|²) |
| BC7 | |U_μ|² | (m_N, |U_μ|²) |
| BC8 | |U_τ|² | (m_N, |U_τ|²) |

Each benchmark assumes a **single Majorana HNL** with dominant mixing to one 
SM flavor. This is a simplification but covers the most important phenomenology.

## 3. Production Mechanisms

### 3.1 Two-Body Leptonic Decays

For **charged pseudoscalar mesons** M⁺:

```
M⁺ → ℓ⁺ N
```

The branching ratio scales as:

```
BR(M → ℓ N) ∝ |U_ℓ|² × (kinematic factor)
```

**Physical channels:**
- K⁺ → ℓ⁺ N (m_N < m_K - m_ℓ ≈ 0.39 GeV for μ)
- D⁺ → ℓ⁺ N (m_N < m_D - m_ℓ ≈ 1.76 GeV for μ)
- D_s⁺ → ℓ⁺ N (m_N < m_Ds - m_ℓ ≈ 1.86 GeV for μ)
- B⁺ → ℓ⁺ N (m_N < m_B - m_ℓ ≈ 5.17 GeV for μ)
- B_c⁺ → ℓ⁺ N (m_N < m_Bc - m_ℓ ≈ 6.16 GeV for μ)

**NOT physical (DO NOT USE):**
- K⁰ → ν N (neutral mesons don't have 2-body leptonic decays to HNL!)
- D⁰ → ν N (same reason)
- B⁰ → ν N (same reason)

### 3.2 Three-Body Semileptonic Decays

For both charged and neutral mesons:

```
M → M' ℓ N
```

**Physical channels:**
- D⁰ → K⁻ ℓ⁺ N
- D⁺ → K⁰ ℓ⁺ N
- B⁰ → D⁻ ℓ⁺ N
- B⁺ → D⁰ ℓ⁺ N
- B_s → D_s⁻ ℓ⁺ N
- Λ_b → Λ_c ℓ⁻ N

These have wider kinematic reach but smaller branching ratios.

### 3.3 Electroweak Production (Drell-Yan)

For m_N > few GeV:

```
pp → W* → ℓ N       (charged current)
pp → Z* → ν N       (neutral current)
pp → tt̄ → Wb Wb → ℓNb + X
```

The cross section is:

```
σ(pp → ℓ N) ≈ σ(pp → W) × BR(W → ℓ ν) × |U_ℓ|²
```

### 3.4 Tau Decays (BC8)

For tau-coupled HNL (BC8), tau leptons from W/Z decay can produce HNLs:

```
τ → π N       (m_N < m_τ - m_π ≈ 1.64 GeV)
τ → ρ N       (m_N < m_τ - m_ρ ≈ 1.00 GeV)
τ → ℓ ν N     (3-body)
```

## 4. Kinematic Considerations

### 4.1 Boost Distribution

For far-detector studies, the HNL boost factor γ = E/m is crucial:

```
λ = βγ × cτ₀
```

where `cτ₀` is the proper decay length and `λ` is the lab-frame decay length.

Meson decays produce **softer** HNLs (lower γ) than EW production, which affects 
detector acceptance.

### 4.2 Angular Distribution

- **Meson production**: Mostly forward (low pT) from soft QCD
- **B mesons from hard QCD**: More central with higher pT
- **EW production**: Follows W/Z rapidity distribution

For **transverse detectors** (MATHUSLA, ANUBIS), central production is important.
For **forward detectors** (FASER), forward production matters.

### 4.3 Production Vertex

HNLs from meson decay inherit the meson's decay vertex:
- **Light mesons (K)**: Decay tens of cm from IP
- **D mesons**: Decay ~mm from IP
- **B mesons**: Decay ~mm from IP

This is handled by Pythia's hadronization and stored in our output.

## 5. Why Pythia Works for Kinematics

The paper arXiv:2103.11494 validated Pythia's approach for HNL production:

> "We conclude that the PYTHIA approach is quite suitable for the estimation 
> of the sensitivity region for the intensity frontier experiments."

Key reasons:
1. Pythia correctly models meson pT and rapidity spectra
2. 2-body decay kinematics are exact
3. 3-body phase space is adequate for acceptance studies
4. The main uncertainty is in absolute cross sections (handled by reweighting)

### 5.1 Limitations

1. **Meson spectra**: FONLL provides more accurate pT distributions
2. **Form factors**: Pythia uses phase space, not full matrix elements
3. **Cross sections**: Should use FONLL or measured cross sections for normalization

For **sensitivity estimates**, Pythia is sufficient. For **final publications**, 
consider supplementing with FONLL-based samples.

## 6. Interface with Stage 2 (HNLCalc)

### 6.1 What This Code Produces

- HNL 4-momentum (pT, η, φ, E)
- Production vertex (x, y, z) in mm
- Parent particle ID
- Event weight
- Boost factor γ

### 6.2 What HNLCalc Provides

- Proper lifetime cτ₀(m_N, |U|²) from analytical formulas
- Decay branching ratios
- Phase space for decay products

### 6.3 Geometry Code

Given an HNL with momentum p⃗ and production point x⃗₀:

1. **Propagation**: The HNL travels in direction p̂ = p⃗/|p⃗|
2. **Decay probability**: 
   ```
   P_decay = exp(-d₁/λ) × (1 - exp(-L/λ))
   ```
   where:
   - `d₁` = distance from IP to detector entry
   - `L` = path length through detector
   - `λ = βγ × cτ₀` = lab-frame decay length

3. **Signal yield**: 
   ```
   N_signal = Σᵢ wᵢ × P_decay,i × ε_i
   ```
   where `ε_i` is the reconstruction efficiency.

## 7. Cross Section Normalization

For correct absolute yields, use these cross sections (14 TeV):

| Process | σ [μb] | Notes |
|---------|--------|-------|
| pp → bb̄ X | 500 | FONLL, inclusive |
| pp → cc̄ X | 10,000 | FONLL, inclusive |
| B⁺ production | ~100 | After fragmentation |
| D⁺ production | ~500 | After fragmentation |
| pp → W X | 200 | NNLO |
| pp → Z X | 60 | NNLO |
| pp → tt̄ X | 0.9 | NNLO |

### 7.1 B-Meson Fractions (from Pythia hadronization)

- B⁰ / B̄⁰: 44.8%
- B⁺ / B⁻: 44.8%
- B_s / B̄_s: 10.3%
- B_c / B_c: 0.018%

Reference: LHCb measurement arXiv:1910.09934

## 8. Validation Checklist

Before using results for a publication:

- [ ] Kinematic distributions (pT, η) match FONLL/MadGraph where available
- [ ] Parent fractions match expected meson ratios
- [ ] Cross-check a few mass points against MATHUSLA/ANUBIS published results
- [ ] Verify kinematic limits are correctly applied
- [ ] Ensure HNL is truly stable in Pythia (no decays)

## 9. References

### Primary Physics References

1. **HNL Phenomenology** (THE reference for production/decay):
   - Bondarenko et al., JHEP 11 (2018) 032
   - arXiv:1805.08567

2. **PBC Benchmarks**:
   - Beacham et al., J.Phys.G 47 (2020) 010501
   - arXiv:1901.09966

3. **Updated Benchmarks**:
   - Abdullahi et al., arXiv:2207.02742

### Experimental/Simulation References

4. **MATHUSLA**:
   - Curtin et al., arXiv:1806.07396
   - GitHub: davidrcurtin/MATHUSLA_LLPfiles_RHN_*

5. **ANUBIS**:
   - Bauer et al., arXiv:1909.13022
   - Hirsch et al., PRD 101 (2020) 055034

6. **Pythia Validation**:
   - Gorkavenko et al., J.Phys.G 48 (2021) 105001
   - arXiv:2103.11494

7. **FONLL for B/D production**:
   - Cacciari et al., JHEP 05 (2012) 007
   - http://www.lpthe.jussieu.fr/~cacciari/fonll/fonllform.html

---

*Document prepared following MATHUSLA methodology for far-detector HNL studies.*
