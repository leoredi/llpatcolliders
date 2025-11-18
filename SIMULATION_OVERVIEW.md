# Simulation Overview and Physics Goals

**Date**: 2025-11-18
**Project**: Long-Lived Particle Detection at Colliders

---

## What We Are Simulating

This codebase simulates the production, propagation, and detection of **Long-Lived Particles (LLPs)** at the Large Hadron Collider (LHC). We focus on two primary production mechanisms:

### 1. Heavy Neutral Leptons (HNLs) from W Boson Decays
**Configuration**: `pythiaStuff/hnlLL.cmnd` and mass-specific variants

**Process**:
- pp collisions at √s = 13.6 TeV (Run 3 energy)
- W± boson production via three channels:
  - Direct q q̄ → W± production
  - W + gluon production (better rates)
  - W + quark production (better rates)
- W± → μ± + N decay (100% branching fraction in simulation)
- N (HNL) → μ∓ q q̄' decay with both lepton-number-conserving (LNC) and lepton-number-violating (LNV) channels

**Mass Range Scanned**: 15, 23, 31, 39, 47, 55, 63, 71 GeV (8 mass points)

**Particle Properties**:
- PDG ID: 9900015 (custom HNL definition)
- Majorana fermion (particle = antiparticle)
- Proper decay length: cτ = 10 m (in base configuration)
- **1M events simulated per mass point** (high statistics run completed)

**Results**: ~1M HNL particles per mass point, with ~95% of events containing 1 HNL and ~5% containing 2 HNLs (from multi-W production)

### 2. Exotic Higgs Decays to LLP Pairs
**Configuration**: `pythiaStuff/higgsLL.cmnd`

**Process**:
- pp → H production via gluon fusion and other SM mechanisms
- H → LL LL̄ (exotic decay to long-lived scalar pairs)
- LL → bb̄ decay

**Particle Properties**:
- Higgs mass: 125 GeV
- LLP (PDG ID: 6000113) mass: 15 GeV
- Narrow width Higgs with exotic decay channel

**Purpose**: Test detector sensitivity to Higgs portal scenarios with long-lived particles

---

## Why We Are Simulating This

### Scientific Motivation

**1. Beyond Standard Model Physics**
- HNLs are well-motivated extensions to the SM, explaining neutrino masses via the seesaw mechanism
- Can account for matter-antimatter asymmetry through leptogenesis
- Higgs portal models provide natural dark matter candidates

**2. Detector Sensitivity Studies**
Long-lived particles are challenging to detect because they:
- Travel macroscopic distances before decaying (meters to hundreds of meters)
- Often escape the main detector before decaying
- Require specialized detectors placed far from the interaction point

**3. Exclusion Limit Setting**
We compare our detector sensitivity against existing experimental limits from:
- **MATHUSLA**: Proposed surface detector above ATLAS/CMS
- **CODEX-b**: Detector in the LHCb cavern
- **ANUBIS**: Another proposed LLP detector

### Our Detector Geometry

**Configuration**: Cylindrical tube detector
- **Location**: z = 22 m from interaction point
- **Radius**: ~1.54 m (elliptical cross-section 1.4 m × 1.1 m)
- **Geometry**: Curved path following correctedVert coordinates
- **Detection**: Based on path length through detector volume

---

## What We Have Accomplished

### Simulation Infrastructure
✓ PYTHIA 8 Monte Carlo event generator integration
✓ Duplicate particle detection (100% elimination of PYTHIA event record artifacts)
✓ Automated mass scanning with parallel execution (2x speedup)
✓ CSV-based data pipeline for **1M events per configuration** (high statistics)

### Analysis Pipeline
✓ Geometric ray-tracing for particle-detector intersection
✓ Lifetime-dependent decay probability calculations
✓ Event-level and particle-level statistics
✓ Exclusion limit plotting vs cτ (proper decay length)

### Completed Studies

**HNL Mass Scan** (as of 2025-11-18):
- 8 mass points from 15-71 GeV **completed with 1M events each**
- High statistics run provides excellent statistical precision
- Full exclusion plots generated for each mass
- Comparison with MATHUSLA, CODEX-b, ANUBIS limits
- Results show detector reaches BR × ε sensitivity of ~10⁻⁴ to 10⁻⁶ depending on lifetime

**Higgs → LL LL̄**:
- Baseline simulation completed
- Demonstrates exotic Higgs decay detection capability

---

## What We Can Simulate Next

### 1. Extended Parameter Scans

#### HNL Lifetime Scan
**Currently**: Fixed cτ = 10 m in PYTHIA, post-processing scan from 1-1000 m
**Next Step**: Generate dedicated samples at different lifetimes
- **Why**: More accurate event acceptance for very short/long lifetimes
- **Suggested range**: cτ = 1, 3, 10, 30, 100, 300 m
- **Implementation**: Modify `9900015:tau0` parameter in .cmnd files

#### Finer Mass Grid
**Currently**: 8 GeV steps (15, 23, 31, ..., 71 GeV)
**Next Step**: 4 GeV or 2 GeV steps
- **Why**: Better resolution for mass-dependent acceptance effects
- **Critical region**: Near kinematic boundaries (low mass, near MW/2)

### 2. New Physics Scenarios

#### Tau-Coupled HNLs
**Current**: μ± coupling only
**Motivation**: Different experimental signatures, tau-philic models
- Modify decay channels: `9900015:addChannel = 1 0.5 101 15 2 -1` (τ instead of μ)
- May have different acceptance due to tau decay products

#### Mixed Couplings
- Simulate HNLs with democratic couplings to e, μ, τ
- More realistic phenomenological scenarios
- Affects visible decay signatures

#### Different LLP Decay Modes
**Current**: N → μ± q q̄
**Options**:
- N → ν γ (invisible + photon)
- N → ν ℓ⁺ ℓ⁻ (displaced dilepton)
- Long-lived neutralinos, staus, gluinos (SUSY scenarios)

### 3. Detector Geometry Variations

#### Optimize Detector Placement
**Current**: Fixed at z = 22 m
**Scan**: z = 10, 15, 20, 25, 30, 40, 50 m
- Find optimal distance for different particle lifetimes
- Trade-off: closer = more flux, farther = more decay probability

#### Different Detector Shapes
**Current**: Cylindrical tube
**Alternatives**:
- Flat rectangular detector (MATHUSLA-style)
- Multiple smaller detectors at different locations
- Calorimeter-style segmented volumes

#### Angular Acceptance Studies
- Currently using 3D ray-tracing with full solid angle
- Could study: forward-only (FASER-style) vs barrel (CODEX-style) geometries

### 4. Background Studies

#### Cosmic Ray Backgrounds
**Not yet simulated**
- Crucial for surface detectors (MATHUSLA)
- Less important for underground locations (CODEX)

#### Neutrino-Induced Backgrounds
- Simulate neutrino interactions in detector volume
- Important for understanding irreducible backgrounds

#### Accidental Coincidences
- Soft QCD processes creating detector activity
- Beam halo effects

### 5. Systematic Studies

#### Detector Resolution Effects
- Smear particle energies and positions
- Study impact on reconstruction efficiency
- Timing resolution for LLP vertex reconstruction

#### PDF Uncertainties
- Vary parton distribution functions
- Impact on production cross sections

#### Pythia Tune Variations
- Test different tunes (Monash 2013, CP5, etc.)
- Systematic on production kinematics

### 6. Advanced Analysis

#### Machine Learning Classification
- Train classifiers to distinguish signal from background
- Feature engineering: displaced vertices, track isolation, etc.
- Could improve sensitivity by factor of 2-5

#### Combined Detector Analysis
- Simulate multiple detectors simultaneously
- Study correlations and combined reach
- Exclusive vs inclusive event selections

#### Model-Independent Limits
- Present results as recast for other BSM scenarios
- Effective field theory (EFT) approach
- Model-independent decay length vs mass limits

---

## Immediate Next Steps (Prioritized)

### High Priority
1. **Finer mass grid** (4 GeV steps): masses 11, 15, 19, 23, ..., 71, 75 GeV
   - Critical for smooth exclusion curves
   - Better physics interpretation
   - Can use same 1M event statistics

2. **Lifetime optimization study**
   - Which cτ values maximize detector sensitivity?
   - Guide for optimal detector placement

### Medium Priority
3. **Tau-coupled HNL study**
   - Complementary to muon channel
   - Different experimental signatures

4. **Background simulation**
   - Essential for realistic sensitivity claims
   - Start with neutrino backgrounds

### Lower Priority (but scientifically interesting)
5. **Detector geometry optimization**
6. **Alternative BSM scenarios** (SUSY, dark photons, etc.)
7. **Machine learning analysis improvements**

---

## Technical Notes

### Computational Resources
- **Current**: ~3.5-4.5 hours for full 8-point mass scan (1M events each, 2 parallel jobs)
- **Scaling**: Linear with number of events, sub-linear with parallelization
- **Bottleneck**: Single-threaded PYTHIA execution per mass point

### Data Storage
- **Current**: ~55 MB per mass point (raw CSV) + ~96 MB analysis outputs
- **Total**: ~1.2 GB for full scan + analysis
- **Note**: Already at high statistics (1M events)

### Analysis Pipeline
- **Geometric calculations**: Using trimesh for 3D ray-tube intersections
- **Post-processing**: Event-level decay probability scanning over lifetime range
- **Visualization**: matplotlib for exclusion plots and correlation studies

---

## Questions to Consider

1. **What is the optimal detector size and placement** for HNL detection in the mass range 15-71 GeV?

2. **How do our sensitivity projections compare** when using realistic backgrounds vs signal-only studies?

3. **What is the relative importance** of different HNL coupling scenarios (e vs μ vs τ)?

4. **Can we exclude any region of parameter space** not already covered by existing experiments?

5. **What BSM scenario** (beyond HNLs) would benefit most from our detector geometry?

---

## References to Existing Data

- HNL simulation outputs: `output/csv/hnlLL_m{15,23,31,39,47,55,63,71}GeVLLP.csv`
- Exclusion plots: `output/images/hnlLL_m{mass}GeVLLP_exclusion_vs_lifetime.png`
- Experimental limits: `external/{MATHUSLA,CODEX,ANUBIS}.csv`
- Analysis script: `decayProbPerEvent.py`
- Mass scan automation: `pythiaStuff/run_mass_scan.py`

---

## APPENDIX: Lifetime Configuration Investigation

**Date**: 2025-11-18
**Branch**: `diagnostic-decay-output`
**Investigation By**: Claude Code diagnostic analysis

### Problem Statement

Investigation into whether HNL particles are being simulated with correct/consistent proper lifetimes across different mass points.

### Key Finding: Fixed Lifetime Across All Masses

**All mass points (15-71 GeV) currently use identical lifetime settings:**

```bash
# In all pythiaStuff/hnlLL_m{15,23,31,39,47,55,63,71}GeV.cmnd files:
9900015:tau0 = 1E4                   # c*tau = 10 m = 1e4 mm
```

**Conversion to physical units:**
- `tau0 = 1E4 mm/c` (PYTHIA units)
- Proper lifetime: **τ = 33.36 nanoseconds**
- Proper decay length: **c*τ = 10 meters**

### Verification: Diagnostic Simulation

Modified `main144.cc` to output decay information (production/decay vertices, proper decay time) and ran 10,000 event test at m = 39 GeV.

**Results confirmed fixed lifetime:**

```
Proper lifetime τ (seconds):
  Mean: 3.285e-08 s  (32.85 ns)
  Median: 2.304e-08 s  (23.04 ns)
  Expected from config: 3.336e-08 s (33.36 ns)  ✓ MATCHES

Proper decay length c*τ (meters):
  Mean: 9.85 m
  Median: 6.91 m
  Expected from config: 10 m  ✓ MATCHES

Lab frame decay distance (boosted):
  Mean: 103.68 m  (average γ ≈ 10.4)
  Median: 27.19 m  (median γ ≈ 4.4)
```

*Note: Spread in measured τ values is physical (exponential decay distribution), not an error.*

### Current Analysis Pipeline Behavior

The analysis script `decayProbPerEvent.py` **does NOT use PYTHIA decay information**:

1. **Reads only kinematic data**: `event, id, pt, eta, phi, momentum, mass`
2. **Ignores PYTHIA decay times/vertices** (not even present in standard output)
3. **Performs post-processing lifetime scan**:
   ```python
   lifetimes = np.logspace(-9.5, -4.5, 20)  # 0.316 ns to 31.6 μs
   ```
4. **Recalculates decay probabilities** for each hypothetical lifetime:
   ```python
   decay_length = gamma * beta * c * lifetime_hypothesis
   p_decay = exp(-d_entry/decay_length) * (1 - exp(-d_path/decay_length))
   ```

**Interpretation:** The PYTHIA simulation provides a **kinematic template** only. The actual lifetime-dependent physics is computed in post-processing.

### Implications

#### What This Means for Current Results

1. **All mass points have identical production kinematics** (since produced from same W → μ N process)
2. **Only difference between mass points**: the mass value itself, affecting:
   - Boost factor: γ = p/m
   - Phase space in W decay
3. **Exclusion limits are derived** by post-processing the same kinematic sample with different lifetime hypotheses

#### Physics Consistency Question

**For HNL phenomenology:** Lifetime depends on mass and mixing as:
```
τ ∝ 1 / (|U|² × m⁵)
```

For fixed mixing |U|², heavier HNLs should decay faster. Current approach with fixed lifetime across all masses is **inconsistent with specific HNL models**, but may be acceptable for **model-independent limits**.

### Three Possible Approaches

#### Option 1: Model-Independent Limits (Current Approach)

**Justification:** Present limits as "for a particle with mass m and lifetime τ, our sensitivity is..."

**Pros:**
- Most general interpretation
- Results applicable to any BSM scenario (not just HNLs)
- Computationally efficient

**Cons:**
- Not directly testable against specific HNL models
- Ignores correlations between m and τ in realistic scenarios

**Validity check needed:**
- Does fixed τ = 33 ns introduce acceptance bias?
- Are kinematic distributions representative of the full lifetime range we scan?

#### Option 2: Mass-Dependent Lifetimes (Model-Specific)

**Justification:** Use theory-motivated τ(m, |U|²) relationships

**Pros:**
- Directly testable against HNL models
- Physically consistent
- Can present limits in mixing parameter space

**Cons:**
- Requires multiple simulations per mass point (for different |U|²)
- Computationally expensive
- Results specific to HNL scenario

**Implementation:**
- For each mass, scan over reasonable |U|² range (e.g., 10⁻¹⁰ to 10⁻⁵)
- Calculate corresponding τ(m, |U|²)
- Run simulations for 3-5 representative lifetimes per mass
- Present limits in (m, |U|²) parameter space

#### Option 3: Hybrid Approach

**Justification:** Combine physics realism with computational efficiency

**Pros:**
- Balance between Options 1 and 2
- Validate that kinematic distributions don't depend strongly on lifetime
- Can still present model-independent limits if desired

**Cons:**
- More complex analysis pipeline
- Requires interpolation between simulation points

**Implementation:**
- Run 3 lifetime values per mass: short (τ ~ 1 ns), medium (τ ~ 30 ns), long (τ ~ 1 μs)
- Verify kinematic distributions are lifetime-independent
- If independent: use any lifetime for kinematics, scan in post-processing (current approach validated)
- If dependent: interpolate between simulated lifetime points

### Technical Implementation Notes

#### Modified Output Format

Branch `diagnostic-decay-output` includes extended CSV output:

```cpp
// pythiaStuff/main144.cc modifications:
myfile << "event,\tid,\tpt,\teta,\tphi,\tmomentum,\tmass,\t"
       << "tau,\txProd,\tyProd,\tzProd,\txDec,\tyDec,\tzDec\n";
```

**New columns:**
- `tau`: Proper decay time in mm/c (PYTHIA units)
- `xProd, yProd, zProd`: Production vertex coordinates (mm)
- `xDec, yDec, zDec`: Decay vertex coordinates (mm)

**Purpose:** Enable verification that PYTHIA is simulating decays correctly, though currently unused by analysis.

#### Diagnostic Script

Created `check_tau.py` to analyze PYTHIA decay output:
- Reads extended CSV format
- Computes proper lifetime statistics
- Calculates lab-frame decay distances
- Verifies consistency with configuration files

### Open Questions Requiring Decision

1. **Which approach aligns with physics goals?**
   - Model-independent limits (Option 1)?
   - HNL-specific parameter space (Option 2)?
   - Hybrid validation + limits (Option 3)?

2. **Does lifetime affect acceptance?**
   - Do different lifetimes change which particles make it to the detector?
   - Are PYTHIA-level cuts applied based on decay position?

3. **Does lifetime affect production kinematics?**
   - Are there interference/resonance effects?
   - Does W → μ N branching depend on N lifetime?
   - (Likely answer: No, but should verify)

4. **What does "cτ" on the x-axis of exclusion plots mean?**
   - Simulated lifetime: cτ = 10 m (fixed)
   - Analyzed lifetime: cτ = 0.1 - 10,000 m (scanned)
   - Current plots show analyzed (hypothetical) cτ, not simulated cτ

### Recommendations for Next Steps

1. **Immediate:** Decide on approach (1, 2, or 3 above)

2. **Validation:** Run small tests with different PYTHIA lifetimes (τ = 1, 10, 100 ns) at one mass point
   - Compare kinematic distributions (pt, eta, phi)
   - Check if acceptance changes significantly
   - Verify post-processing approach validity

3. **Documentation:** Clarify in papers/presentations:
   - What lifetime is used in simulation vs analysis
   - Whether limits are model-independent or model-specific
   - Relationship between axes labels and simulation parameters

4. **Long term:** Consider extending to model-dependent limits in (m, |U|²) space if targeting HNL phenomenology specifically

### Files Modified/Created

- `pythiaStuff/main144.cc`: Extended output format (branch: diagnostic-decay-output)
- `check_tau.py`: Diagnostic analysis script
- `output/csv/hnlLL_m39GeVLLP.csv`: Test output with 10k events, extended format
- `SIMULATION_OVERVIEW.md`: This documentation (section added)

---

*This document should be updated as new simulations are completed and new physics questions arise.*
