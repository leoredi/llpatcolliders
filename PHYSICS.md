# Physics Documentation

## 1. Overview

This project computes exclusion sensitivity for Heavy Neutral Leptons (HNLs, also called sterile neutrinos or right-handed neutrinos) at a proposed long-lived particle detector near CMS at the HL-LHC. The detector is a drainage gallery — a tube-shaped cavity at O(100 m) from the interaction point.

An HNL with mass m_N and mixing |U_ℓ|² to a lepton flavour ℓ ∈ {e, μ, τ} is produced in meson/baryon decays or W/Z decays, travels macroscopic distances (because its coupling is small), and decays inside the gallery into visible charged particles. The analysis finds, for each mass, the range of |U_ℓ|² where the expected signal exceeds the 95% CL Poisson threshold (N = 2.996) at 3000 fb⁻¹ — the "exclusion island".

The lower bound on |U|² comes from the HNL being too long-lived to decay inside the gallery (exponential decay suppression). The upper bound comes from the HNL being too short-lived and decaying before reaching the gallery.

## 2. HNL production channels

### 2.1 Meson production (Pythia, low mass)

HNLs below ~5 GeV are dominantly produced from meson and baryon decays at the LHC interaction point. The production is simulated with Pythia 8.315 in three mass regimes:

**Kaon regime (m_N < 0.5 GeV):**
- K± → ℓ N (2-body, charged current)
- K_L → π ℓ N (3-body, semileptonic)
- Uses card file `hnl_Kaon.cmnd`

**Charm regime (0.5 < m_N < 2.0 GeV):**
- D± → ℓ N (2-body)
- D_s± → ℓ N (2-body)
- D⁰ → K ℓ N (3-body, semileptonic)
- Λ_c → Λ ℓ N (3-body, baryonic)
- Uses card file `hnl_Dmeson.cmnd`

**Beauty regime (m_N > 2.0 GeV):**
- B± → ℓ N (2-body)
- B_c± → ℓ N (2-body)
- B⁰ → D ℓ N, B_s → D_s ℓ N (3-body)
- Λ_b → Λ_c ℓ N (3-body, baryonic)
- Uses card file `hnl_Bmeson.cmnd`

Each decay channel is configured with 100% branching ratio in Pythia (exclusive mode). The physical branching ratio is computed analytically later by HNLCalc and applied as a weight.

### 2.2 Tau chain production (fromTau mode)

For tau-coupled HNLs with m_N < 1.77 GeV (the threshold used in the production script, slightly below m_τ = 1.777 GeV), a second production path exists:

```
parent meson → τ ν_τ    (SM decay)
                τ → N X  (BSM decay, X = π, ρ, 3π, ℓν)
```

This "fromTau" mode is implemented as a separate Pythia run. The tau decay channels are configured with adaptive branching ratios that depend on the available phase space:
- τ → ρ N (dominant when kinematically allowed)
- τ → 3π N
- τ → π N
- τ → μ ν N, τ → e ν N (leptonic)

The parent mesons that produce taus are D_s± → τ ν (BR = 5.3%) and B mesons → τ ν X (BR ≈ 2.3% each for B⁰, B⁺, B_s⁰). These BR values are hardcoded in `get_parent_tau_br()` in `production_xsecs.py`.

### 2.3 Electroweak production (MadGraph, high mass)

For m_N above a few GeV, electroweak production pp → W/Z → ℓ N becomes relevant:

- pp → W± → ℓ± N (dominant)
- pp → Z → ν N (subdominant)

Generated with MadGraph5 using the `SM_HeavyN_CKM_AllMasses_LO` model at √s = 14 TeV. Events are generated with |U|² = 1 and rescaled analytically. Note: an NLO K-factor of 1.3 is defined in `run_hnl_scan.py` but is not currently applied in the analysis pipeline.

In the limit calculation, EW normalization uses fixed σ(W) and σ(Z) values from `production_xsecs.py`; MadGraph event weights are used for kinematic averaging (efficiency), not as the absolute cross-section normalisation.

The LHE output is converted to CSV with the same column format as Pythia production for pipeline compatibility.

### 2.4 Production cross-sections and fragmentation

The analysis weights each event by the physical production rate. Key cross-sections at 14 TeV:

| Process | σ | In code (pb) |
|---|---|---|
| pp → cc̄ | 24 mb | 2.4×10¹⁰ |
| pp → bb̄ | 500 μb | 5.0×10⁸ |
| pp → K± (inclusive) | 50 mb | 5.0×10¹⁰ |
| pp → K_L (inclusive) | 25 mb | 2.5×10¹⁰ |
| pp → W | — | 2.0×10⁸ |
| pp → Z | — | 6.0×10⁷ |

Charm fragmentation fractions: D⁰ 59%, D⁺ 24%, D_s⁺ 10%, Λ_c⁺ 6%.

Beauty fragmentation fractions: B⁰ 40%, B⁺ 40%, B_s⁰ 10%, Λ_b⁰ 10%.

Additional rare b-hadrons: B_c⁺ (0.1%), Ξ_b (3%), Ω_b (1%) — these are included in the code as small fractions of σ(bb̄).

These are defined in `analysis_pbc/config/production_xsecs.py`.

## 3. HNL decay physics (HNLCalc)

The HNLCalc module (`analysis_pbc/HNLCalc/HNLCalc.py`) computes the total decay width Γ_tot and proper lifetime cτ₀ for an HNL of given mass and couplings (U_e, U_μ, U_τ).

### 3.1 Decay channels

**2-body decays (N → hadron + lepton):**
- N → ℓ + pseudoscalar (π, K, D, D_s, B, B_c)
- N → ℓ + vector (ρ, K*, D*, ω, φ)
- N → ν + pseudoscalar (π⁰, η, η', K⁰)
- N → ν + vector (ρ⁰, ω, φ)

**3-body leptonic decays:**
- N → ℓ_α ℓ_β ν (charged current)
- N → ν ℓ⁺ℓ⁻ (neutral current)
- N → ν ν ν (invisible)

**3-body hadronic decays (m_N > 1 GeV):**
- N → ℓ u d̄ (charged current, with QCD corrections)
- N → ν q q̄ (neutral current, with QCD corrections)

### 3.2 Form factors

Semileptonic meson decays use form factors with pole parameterisation:

**Pseudoscalar daughters:** two form factors f₊(q²), f₀(q²) with vector and scalar pole masses. Example: D → K has f₀₀ = 0.747, M_V = M(D*_s) = 1.968 GeV.

**Vector daughters:** four form factors V, A₀, A₁, A₂ with monopole or dipole poles. Example: B → D* has V₀ = 0.76, A₁₀ = 0.66.

**Baryon daughters:** six form factors f₁–f₃, g₁–g₃ with vector and axial pole masses. Example: Λ_c → Λ has f₁₀ = 0.29, g₁₀ = 0.38.

### 3.3 CKM matrix elements

CKM elements used in production and decay:
- V_ud = 0.97373, V_us = 0.2243, V_ub = 3.82×10⁻³
- V_cd = 0.221, V_cs = 0.975, V_cb = 40.8×10⁻³

These are hardcoded in HNLCalc and mapped to specific parent→daughter meson transitions.

### 3.4 Scaling properties

The key scaling laws that enable fast parameter scans:
- **Lifetime:** cτ₀(|U|²) = cτ₀_ref × (|U|²_ref / |U|²)
- **Production BR:** BR(|U|²) = BR_ref × (|U|² / |U|²_ref)

Both follow from the fact that all production and decay amplitudes are proportional to |U|². The scaling check script enforces a default relative-error tolerance of 5×10⁻⁴ (0.05%) in `check_hnlcalc_scaling.py`.

## 4. Detector geometry

The detector is a drainage gallery near CMS, modelled as a tube with:
- 47 path vertices defining the centreline in the x-y plane at z = 22 m
- Tube radius: 1.4 m × 1.1 safety margin = 1.54 m effective radius
- Distance from IP: O(100 m) along most of the path

The geometry is implemented as a 3D triangle mesh (trimesh) in `per_parent_efficiency.py`. For each HNL direction (η, φ), a ray is cast from the origin and intersected with the mesh to find:
- **entry_distance**: distance from IP to where the ray enters the tube (metres)
- **path_length**: distance the ray travels inside the tube (metres)
- **hits_tube**: boolean flag

These are cached per production CSV to avoid recomputation.

## 5. Signal calculation

The expected number of signal events for a given (mass, |U|², flavour) is computed in `expected_signal.py`:

```
N_sig = Σ_parents  L × σ_parent × BR(parent → HNL) × ε_parent
```

where ε_parent is the weighted-average detector efficiency for that parent species:

```
ε_parent = Σ_i (w_i × P_decay_i × P_sep_i) / Σ_i w_i
```

and:
- L = 3000 fb⁻¹ (HL-LHC integrated luminosity)
- σ_parent = production cross-section for the parent particle (in pb; a factor 10³ converts L[fb⁻¹] × σ[pb] to events)
- BR(parent → HNL) = production branching ratio from HNLCalc
- w_i = event weight from Pythia/MadGraph (normalises the MC sample)
- P_decay_i = exp(-d_entry / λ) × [1 - exp(-d_path / λ)] with λ = βγ × cτ₀
- P_sep_i = 1 if ≥ 2 charged decay products separated by ≥ 1 mm at the detector, else 0

For tau-chain production, the "parent" in the sum is the grandfather meson (D_s, B, etc.), and the signal formula becomes:

```
N_sig_tau = Σ_grandparents  L × σ_grandparent × BR(grandparent → τν) × BR(τ → N X) × ε
```

where σ_grandparent is the production cross-section of the meson that decays to τ (not the tau itself), BR(grandparent → τν) is the SM branching ratio from `get_parent_tau_br()`, and BR(τ → N X) comes from HNLCalc (keyed as parent PDG 15).

**W/Z production BRs:** The W and Z branching ratios to HNL are computed in `hnl_model_hnlcalc.py` using the **total** mixing (|U_e|² + |U_μ|² + |U_τ|²), not a flavour-specific coupling. For pure-coupling benchmarks (e.g., 100, 010, 001) this is equivalent, but for mixed-coupling scenarios the distinction matters.

Dirac HNLs get a factor of 2 compared to Majorana (both N and N̄ contribute).

### 5.1 Exclusion threshold

The exclusion limit is set at N_sig = 2.996 events, corresponding to the 95% CL Poisson upper limit for zero observed events (−ln 0.05 = 2.996).

### 5.2 eps2 scan

For each mass point, |U|² is scanned over [10⁻¹², 10⁻²] on a log-uniform grid of 100 points. By default (`limits/run.py`), HNLCalc is evaluated once at a reference eps2 and then scaled analytically across the scan; `--hnlcalc-per-eps2` disables this and recomputes HNLCalc at every point. The lower and upper crossings of N_sig = 2.996 define eps2_min and eps2_max — the boundaries of the exclusion island.

Linear interpolation in log-space refines the crossing points to sub-grid precision.

## 6. Decay acceptance

The decay acceptance calculation in `decay_detector.py` determines whether an HNL decay produces detectable tracks.

### 6.1 Decay product kinematics

Pre-computed HNL decay events from the MATHUSLA RHN libraries are loaded and boosted from the HNL rest frame to the lab frame using the HNL's βγ and direction. (A separate script can generate custom decay samples with MadGraph+Pythia, but this is not the default runtime path.)

### 6.2 Track separation cut

A minimum separation of 1 mm between charged decay products at the detector surface rejects backgrounds from neutrino interactions. The separation is computed by projecting charged particle directions from the decay vertex onto the detector mesh.

### 6.3 Decay file selection

The MATHUSLA RHN decay libraries provide pre-computed decay events at discrete mass points. File selection is flavour-dependent:
- Below the low-mass threshold (0.42 GeV for e/τ, 0.53 GeV for μ): prefer analytical files, with nearest-file fallback
- Above threshold: filter to the flavour's allowed category set, then choose the nearest mass within that filtered set
- Warn if |Δm| > 0.5 GeV

## 7. Output and plots

The final output is `output/csv/analysis/HNL_U2_limits_summary.csv` with columns:
- mass_GeV, flavour, benchmark
- eps2_min, eps2_max (boundaries of the exclusion island, NaN if no sensitivity)
- peak_events (maximum N_sig over the scan)
- separation_mm (track-separation cut used for that run)

If `--timing` is enabled, additional timing/count columns are also written.

The money plot (`plot_money_island.py`) shows the exclusion island in the (m_N, |U_ℓ|²) plane for each flavour. The shaded region is excluded. The lower boundary ("too long-lived") is eps2_min. The upper boundary ("too prompt") is eps2_max. A tip-point interpolation closes the island at the high-mass sensitivity edge.
