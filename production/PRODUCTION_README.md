# HNL Production Pipeline - Complete Guide

## Quick Start

```bash
# Single mass point
./main_hnl_production <mass_GeV> <flavor> [nEvents] [mode]

# Examples
./main_hnl_production 2.6 muon 200000              # Muon coupling, beauty regime
./main_hnl_production 3.0 tau 200000 direct        # Tau MODE A (direct production)
./main_hnl_production 1.0 tau 200000 fromTau       # Tau MODE B (from tau decay)
```

## Output Files

CSV files created in current directory with format:
```
HNL_{mass}GeV_{flavor}_{regime}[_{mode}].csv
```

**Examples:**
- `HNL_2p60GeV_muon_beauty.csv` - 2.6 GeV muon, B-meson production
- `HNL_3p00GeV_tau_beauty_direct.csv` - 3.0 GeV tau direct mode
- `HNL_1p00GeV_tau_charm_fromTau.csv` - 1.0 GeV tau fromTau mode

**Regime naming:**
- `kaon`: m < 0.5 GeV (K-meson dominated)
- `charm`: 0.5 ≤ m < 2 GeV (D-meson dominated)
- `beauty`: 2 ≤ m < 5 GeV (B-meson dominated)
- `ew`: m ≥ 5 GeV (W/Z boson dominated)

---

## Enhanced Mass Grid (202 Total Files)

Production uses **enhanced mass grids** with increased sampling near kinematic thresholds for publication-quality exclusion curves.

### Electron Coupling (56 mass points)

```bash
0.2 0.22 0.25 0.28 0.3 0.32 0.35 0.38 0.40 0.42 0.45 0.48 0.50 0.52 0.55
0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.75 1.8 1.82 1.85 1.9
2.0 2.3 2.6 3.0 3.4 3.8 4.2 4.6 4.8 5.0 5.2
6.0 7.0 8.0 9.0 10.0 12.0 15.0 20.0 30.0 40.0 50.0 60.0 80.0 GeV
```

### Muon Coupling (56 mass points)

```bash
0.2 0.22 0.25 0.28 0.3 0.32 0.35 0.37 0.38 0.39 0.40 0.42 0.45 0.48 0.50 0.55
0.6 0.7 0.8 0.9 1.0 1.2 1.4 1.6 1.65 1.70 1.75 1.76 1.78 1.80 1.85 1.90
2.0 2.3 2.6 3.0 3.4 3.8 4.2 4.6 4.8 5.0 5.2
6.0 8.0 10.0 12.0 15.0 20.0 30.0 40.0 50.0 60.0 80.0 GeV
```

### Tau Coupling (45 mass points, dual-mode)

```bash
# Direct mode (45 points):
0.50 0.55 0.60 0.65 0.70 0.80 0.90 1.00 1.10 1.20 1.30 1.40 1.45 1.50 1.55
1.60 1.62 1.64 1.66 1.70 1.74 1.78 1.80 1.85 1.90
2.0 2.4 2.8 3.0 3.2 3.6 4.0 4.5
6.0 7.0 8.0 10.0 15.0 20.0 30.0 40.0 50.0 60.0 80.0 GeV

# FromTau mode (19 points, m < 1.6 GeV only):
0.50 0.55 0.60 0.65 0.70 0.80 0.90 1.00 1.10 1.20 1.30 1.40 1.45 1.50 1.55
1.60 1.62 1.64 1.66 GeV
```

**Total Production Files:**
- Electron: 56 files
- Muon: 56 files
- Tau: 45 direct + 19 fromTau = 64 files
- **Grand Total: 176 files (standard grid) + enhanced sampling = 202 planned**

### Enhanced Threshold Sampling

**1. K Threshold (m_K - m_ℓ ≈ 0.39-0.49 GeV)** - Kaon → D-meson transition
- Enhanced muon sampling: 0.37, 0.38, 0.39, 0.40, 0.42, 0.45, 0.48, 0.50 GeV
- Captures BR(K→μN) → BR(D→μN) crossover

**2. D Threshold (m_D - m_ℓ ≈ 1.65-1.87 GeV)** - Charm → beauty transition
- Enhanced muon sampling: 1.65, 1.70, 1.75, 1.76, 1.78, 1.80, 1.85, 1.90 GeV
- Captures BR(D→μN) → BR(B→μN) crossover

**3. B Threshold (m_B - m_ℓ ≈ 4.8-5.3 GeV)** - Beauty → EW transition
- Enhanced sampling: 4.6, 4.8, 5.0, 5.2, 6.0 GeV
- Captures meson → W/Z transition

**4. Tau Threshold (m_τ = 1.777 GeV)** - Dual-mode boundary
- Enhanced tau sampling: 1.40, 1.45, 1.50, 1.55, 1.60, 1.62, 1.64, 1.66, 1.70, 1.74, 1.78, 1.80 GeV
- Precisely maps fromTau kinematic cutoff

**Physics Benefits:**
- Captures sharp features in exclusion limit curves
- Maps geometric acceptance × BR transitions
- Identifies sensitivity "island" boundaries
- Enables publication-quality interpolation

---

## Dual-Mode Tau Production

For tau coupling, **TWO production modes** are run to maximize physics reach:

### MODE A: Direct Production (direct)
```bash
./main_hnl_production 2.0 tau 200000 direct
```

**Physics:** O(Uτ²) mixing at production vertex
- B/Ds/Bc → τ N (2-body leptonic)
- B → D τ N (3-body semileptonic)
- W/Z → τ N (electroweak, m ≥ 5 GeV)
- **Parent PDG:** 511 (B⁰), 521 (B±), 531 (Bs), 411 (D±), 431 (Ds), 24 (W±), 23 (Z)

### MODE B: From Tau Decay (fromTau)
```bash
./main_hnl_production 1.0 tau 200000 fromTau
```

**Physics:** O(Uτ²) mixing at tau decay vertex
- B/W → τ ν (Standard Model production of tau)
- τ → π N (forced decay to HNL)
- **Parent PDG:** 15 (τ) - NOT meson!
- **Kinematic limit:** m_HNL + m_π < m_τ → m_HNL < 1.637 GeV

**Important:** Modes are kept separate to avoid O(Uτ⁴) contamination. Combine in analysis:
```
N_total(tau coupling) = N_direct + N_fromTau
```

### Mode Selection by Mass

| Mass Range | Direct Mode | FromTau Mode | Notes |
|------------|-------------|--------------|-------|
| 0.5-1.6 GeV | B/D mesons | τ → π N available | Both modes contribute |
| 1.6-5.0 GeV | B/Bc mesons | Kinematically forbidden | Direct only |
| ≥ 5.0 GeV | W/Z bosons | Kinematically forbidden | Direct only |

---

## Parallel Production

**DO NOT run sequentially** - use parallel execution to reduce total runtime from ~150 hours to ~1-2 hours.

### Example: 10 Parallel Workers

```bash
#!/bin/bash
NEVENTS=200000
MAX_PARALLEL=10

wait_for_slot() {
    while [ $(jobs -r | wc -l) -ge $MAX_PARALLEL ]; do
        sleep 2
    done
}

# Electron masses
for mass in 0.2 0.25 0.3 0.35 0.4 0.45 0.5 ...; do
    wait_for_slot
    ./main_hnl_production $mass electron $NEVENTS > "electron_${mass}.log" 2>&1 &
done

# Muon masses
for mass in 0.2 0.25 0.3 0.35 0.4 0.45 0.5 ...; do
    wait_for_slot
    ./main_hnl_production $mass muon $NEVENTS > "muon_${mass}.log" 2>&1 &
done

# Tau direct mode
for mass in 0.5 0.7 1.0 1.3 1.6 2.0 2.4 ...; do
    wait_for_slot
    ./main_hnl_production $mass tau $NEVENTS direct > "tau_direct_${mass}.log" 2>&1 &
done

# Tau fromTau mode (m < 1.6 GeV only)
for mass in 0.5 0.7 1.0 1.3 1.6; do
    wait_for_slot
    ./main_hnl_production $mass tau $NEVENTS fromTau > "tau_fromTau_${mass}.log" 2>&1 &
done

wait  # Wait for all jobs to complete
```

**Resource Requirements:**
- **CPU:** 10 cores for 10 parallel jobs
- **Memory:** ~500 MB per job → 5 GB total
- **Runtime:** ~1-2 hours for full grid (vs 120-150 hours sequential)
- **Disk:** ~1.1 GB final output (170 files × 5-10 MB each)

---

## Production Status (Current)

**As of Nov 28, 2025:**

```
✅ Completed: 170 files (84% of 202 planned)
❌ Failed:     32 files (16% - all high mass ≥5 GeV)
```

### Working Mass Ranges

**Complete coverage (0.2 - 4.6 GeV):**
- Kaon regime (m < 0.5 GeV): ✅ 100%
- Charm regime (0.5-2 GeV): ✅ 100%
- Beauty regime (2-5 GeV): ✅ 100%

### Failed Mass Ranges (Z Decay Bug)

**High-mass EW regime (m ≥ 5 GeV):** ❌ 0% success

Missing files by flavor:
- Electron: 5.0, 5.2, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 15.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0 GeV (15 files)
- Muon: 5.0, 5.2, 6.0, 8.0, 10.0, 12.0, 15.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0 GeV (13 files)
- Tau: 6.0, 7.0, 8.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0 GeV (11 files)
- **Total missing: 39 files**

---

## Benchmark Statistics

**2.6 GeV Muon (BC7 - Beauty Regime):**

From actual production run (200k events):
```
Events generated:     200,000
HNLs found:           399,956  (200% efficiency - normal!)
File size:            ~39 MB
Runtime:              ~70 minutes (single core)

Parent composition:
  B⁰/B̄⁰:  43.1%  (PDG ±511)
  B±:      43.1%  (PDG ±521)
  Bs/B̄s:   9.7%  (PDG ±531)
  Λb/Λ̄b:   3.9%  (PDG ±5122)
  Bc±:     0.1%  (PDG ±541)
```

**Notes:**
- Efficiency > 100% is expected: Multiple HNLs per pp collision
- Typical: 2-4 HNLs per event in B-meson regime
- Particle/antiparticle balance: ±0.5% (expected from CP-even pp collisions)

---

## Critical Physics Implementation

### HNL Particle Definition

```cpp
// In main_hnl_production.cc
pythia.readString("9900015:all = ν₄ nu4 2 0 0 " + mass + " 0.0 0.0 0.0 0.0");
//                                    ↑
//                               Spin type = 2 (fermion)
```

**Parameters:**
- **PDG code:** 9900015 (user-defined fermion)
- **Spin type:** 2 (spin-1/2 fermion, NOT scalar)
- **Mass:** User-specified (0.2-80 GeV)
- **Lifetime:** Stable in Pythia (`mayDecay = off`)
  - Decays handled in Stage 2 (Python geometry + HNLCalc)
- **Charge:** 0 (neutral)

### Production Normalization (Critical!)

**What Pythia Does:**
- Generates pp collision kinematics at √s = 14 TeV
- Forces parent meson decays: BR(M→ℓN) = 1.0 for efficient sampling
- Outputs HNL 4-momenta, production vertex, parent PDG

**What Pythia Does NOT Do:**
- Apply real branching ratios (would be ~10⁻⁹, too rare)
- Apply production cross-sections (internal σ differs from NLO theory)
- Simulate HNL decays (lifetime calculated in analysis)

**Stage 2 Analysis Applies:**
- External cross-sections (σ_K, σ_D, σ_B, σ_W from PBC literature)
- Real branching ratios from HNLCalc: BR(M→ℓN) ∝ |U_ℓ|²
- Geometric acceptance via ray-tracing
- Decay probability from boost factors and HNLCalc lifetimes

**Result:** Pythia generates kinematics ONLY. Physics normalization in analysis.

---

## CSV File Format

**Header:**
```csv
event,weight,hnl_id,parent_pdg,pt,eta,phi,p,E,mass,prod_x_mm,prod_y_mm,prod_z_mm,boost_gamma
```

**Columns:**
- `event`: Pythia event number (0-indexed, repeats for multiple HNLs per event)
- `weight`: Relative MC weight (all = 1.0, unweighted generation)
- `hnl_id`: HNL PDG code (9900015 or -9900015)
- `parent_pdg`: Parent particle PDG code **[CRITICAL FOR ANALYSIS]**
  - Examples: 511 (B⁰), 521 (B±), 411 (D⁺), 431 (Ds), 321 (K±), 24 (W±), 23 (Z)
- `pt`, `eta`, `phi`: Transverse momentum [GeV], pseudorapidity, azimuthal angle
- `p`, `E`, `mass`: 3-momentum magnitude [GeV], energy [GeV], mass [GeV]
- `prod_x_mm`, `prod_y_mm`, `prod_z_mm`: Production vertex [mm from IP]
- `boost_gamma`: Boost factor γ = E/m (for lifetime calculations)

**Example line:**
```csv
44,1.0,9900015,511,1.75,3.07,-2.75,18.78,18.96,2.6,0.19,0.023,2.3,7.29
```
Event 44 produced HNL with mass 2.6 GeV from B⁰ (PDG 511), momentum 18.78 GeV, boost γ=7.29.

---

## Known Issues

### CRITICAL: Z Boson Decay Failure (UNFIXED)

**Status:** ❌ **BUG NOT FIXED** - All high-mass production fails

**Symptom:**
```
PYTHIA Error in ResonanceDecays::next: no open decay channel for id = 23
Events generated: 0
HNLs found: 0
```

**Affected Masses:**
- ALL masses ≥ 5 GeV (EW regime)
- 100% failure rate: 0/39 high-mass files succeeded
- Low/mid masses (< 5 GeV) unaffected

**Root Cause:**

In `main_hnl_production.cc` lines 422-425 (`configureEWDecays()` function):

```cpp
// Z -> ν N (for neutral current production)
if (mHNL < mZ) {
    pythia.readString("23:onMode = off");  // ← Turns OFF all Z decays
    pythia.readString("23:addChannel = 1 0.5 0 " + nu + " " + hnl);
    pythia.readString("23:addChannel = 1 0.5 0 " + nuBar + " " + hnl);
}
```

**The Problem:**
1. Code disables ALL Z boson decay channels
2. Adds ONLY: Z → ν N and Z → ν̄ N (custom decays)
3. At m ≥ 5 GeV, EW card enables Z production: `WeakSingleBoson:ffbar2gmZ = on`
4. When Z → ν N decay fails, Z has ZERO open channels
5. Pythia aborts: "no open decay channel"

**Why Z → ν N Fails:**
- Pythia's resonance decay handler may not support user-defined fermions in Z decays
- Phase space calculation issues
- Numerical instabilities in matrix elements
- BR=1.0 forcing may conflict with resonance kinematics

**Impact:**
- Missing 39/202 files (19% of planned production)
- No coverage of EW regime (m ≥ 5 GeV)
- Physics-critical low/mid-mass range (0.2-5 GeV) **completely unaffected**

**Workaround Options:**
1. Keep some standard Z decay channels open as fallback
2. Disable Z production entirely, rely only on W→ℓN
3. Fix Pythia's handling of custom fermion Z decays
4. Accept lower coverage at high masses (detector sensitivity low there anyway)

**Files for Bug Fix:**
To fix this issue, share with another LLM:
1. `main_hnl_production.cc` (lines 417-450, `configureEWDecays()` function)
2. `cards/hnl_EW.cmnd` (shows Z production enabled)
3. Example error log from `../output/logs/simulation_new/HNL_10p00GeV_electron_ew.log`
4. This README section

---

### Issue 2: Low Efficiency Warnings

**Symptom:**
```
HNLs found: 234 (0.1% efficiency)
```

**Causes:**
1. Mass too high for regime (e.g., m=5 GeV with kaon card)
2. Kinematically forbidden (e.g., tau at m=0.3 GeV where m < m_τ - m_π)
3. FromTau mode at m > 1.6 GeV (tau decay forbidden)

**Solution:**
- Check mass is appropriate for regime (code auto-selects, but verify in logs)
- Tau fromTau mode only valid for m < 1.6 GeV

---

### Issue 3: File Naming

**Current format (Nov 2025):**
```
HNL_{mass}GeV_{flavor}_{regime}[_{mode}].csv
```

Examples:
- `HNL_2p60GeV_muon_beauty.csv` (2.6 GeV → "2p60")
- `HNL_10p00GeV_electron_ew.csv` (10.0 GeV → "10p00")

**Old format (deprecated):**
```
HNL_mass_{mass}_{flavor}_{regime}.csv
```

If analysis code expects old format, update filename parsing.

---

## Validation Checklist

After production, verify:

```bash
# 1. File count
ls ../output/csv/simulation_new/*.csv | wc -l
# Expected: 170 (current), 202 (after Z bug fix)

# 2. Check for empty files
find ../output/csv/simulation_new/ -name "*.csv" -size 0

# 3. Check HNL counts
for f in ../output/csv/simulation_new/*.csv; do
    echo "$f: $(tail -n +2 $f | wc -l) HNLs"
done

# 4. Verify parent PDG codes make sense
# Example: 2.6 GeV should have B mesons (511, 521, 531, 5122)
awk -F',' 'NR>1 {print $4}' ../output/csv/simulation_new/HNL_2p60GeV_muon_beauty.csv | sort | uniq -c

# 5. Check for NaN/Inf
grep -E "nan|inf|NaN|Inf" ../output/csv/simulation_new/*.csv
```

**Expected Results:**
- 170 files with non-zero size
- No NaN/Inf values
- Parent PDG codes match regime (511/521/531 for beauty, 411/431 for charm, etc.)
- Boost factors γ > 1 (HNLs are boosted)
- Production vertices near IP (|x|, |y|, |z| < ~1 cm for most events)

---

## Next Steps After Production

### 1. Move Files to Analysis Directory
```bash
mv HNL_*.csv ../output/csv/simulation_new/
```

### 2. Run Analysis Pipeline
```bash
cd ../analysis_pbc_test

# Geometry preprocessing (computes ray-tracing, cached)
conda run -n llpatcolliders python geometry/per_parent_efficiency.py

# Calculate exclusion limits
conda run -n llpatcolliders python limits/u2_limit_calculator.py

# Generate plots
conda run -n llpatcolliders python ../money_plot/plot_money_island.py
```

### 3. Verify Output
```bash
# Check analysis results
cat ../output/csv/analysis/HNL_U2_limits_summary.csv

# View final plot
open ../output/images/HNL_moneyplot_island.png
```

---

## Technical Details

### Code Version
- **Date:** November 2025
- **Pythia version:** 8.315 (bundled in `pythia/pythia8315/`)
- **Compiler:** g++ with C++11
- **Platform:** macOS (Darwin 25.1.0)

### Features
- ✅ Dual-mode tau production (direct + fromTau)
- ✅ Enhanced mass grids with threshold sampling
- ✅ Correct HNL spin type (fermion, not scalar)
- ✅ Parent tracking algorithm
- ✅ Zero compiler warnings
- ✅ Parallel production capability
- ❌ Z decay bug at high mass (UNFIXED)

### Regime Boundaries
- **Kaon:** m < 0.5 GeV (SoftQCD processes)
- **Charm:** 0.5 ≤ m < 2.0 GeV (`HardQCD:hardccbar`)
- **Beauty:** 2.0 ≤ m < 5.0 GeV (`HardQCD:hardbbbar`)
- **EW:** m ≥ 5.0 GeV (`WeakSingleBoson`, W/Z production)

Boundaries chosen to match dominant production mechanism at each mass.

---

## References

**Code Files:**
- `main_hnl_production.cc` - Main simulation code
- `run_production_benchmark.sh` - Validation script (2.6 GeV muon)
- `cards/hnl_*.cmnd` - Pythia configuration files

**Analysis Pipeline:**
- `../analysis_pbc_test/limits/u2_limit_calculator.py` - Limit calculation
- `../analysis_pbc_test/models/hnl_model_hnlcalc.py` - HNLCalc wrapper
- `../money_plot/plot_money_island.py` - Exclusion plot generator

**Documentation:**
- `../CLAUDE.md` - Complete project-level guide for LLMs
- `../llm_validation_guide/production/README_PRODUCTION.md` - Validation procedures
