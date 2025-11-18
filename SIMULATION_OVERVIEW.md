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

*This document should be updated as new simulations are completed and new physics questions arise.*
