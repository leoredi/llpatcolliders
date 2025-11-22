# PBC-Grade HNL Analysis Pipeline

Physics Beyond Colliders (PBC) compatible analysis for the CMS drainage gallery long-lived particle detector.

## Quick Start

### 1. Installation

```bash
# Ensure conda environment exists
conda activate llpatcolliders

# Install HNLCalc (if not already installed)
cd analysis_pbc_test
git clone https://github.com/laroccod/HNLCalc.git
conda run -n llpatcolliders pip install sympy mpmath particle numba
conda run -n llpatcolliders pip install 'scikit-hep==0.4.0'
```

### 2. Run Pipeline Tests

```bash
# From llpatcolliders/analysis_pbc_test/
conda run -n llpatcolliders python tests/test_pipeline.py
```

**Expected output:**
```
✓ TEST 1 passed (HNL model wrapper)
✓ TEST 2 passed (Geometry preprocessing)
✓ TEST 3 passed (Expected signal events)
```

### 3. Calculate |U|² Limits

```bash
# For all mass points (muon coupling)
conda run -n llpatcolliders python limits/u2_limit_calculator.py

# Or run benchmark test (2.6 GeV muon)
conda run -n llpatcolliders python tests/test_26gev_muon.py
```

---

## Pipeline Overview

The analysis follows a **three-stage** workflow matching PBC methodology (MATHUSLA/CODEX-b/ANUBIS):

### Stage 1: Production (Pythia C++)
**Location:** `../production/main_hnl_single.cc`

- Pythia 8.315 generates pp collisions at √s = 14 TeV
- Two regimes: Mesons (K/D/B, m < 5 GeV) and EW (W/Z, m ≥ 5 GeV)
- Forced decays: BR(M→ℓN) = 1.0 for efficient sampling
- HNL is **stable** in Pythia (`mayDecay = off`)
- Output: CSV with kinematics and `parent_id` (PDG code)

**Key:** Real BRs come from HNLCalc in Stage 3, not Pythia.

### Stage 2: Geometry (Python)
**Location:** `geometry/per_parent_efficiency.py`

- Ray-traces each HNL trajectory from IP (0,0,0) through detector mesh
- Detector: Tube at z = 22m above CMS (drainage gallery)
- Computes: entry distance, path length, hits detector (boolean)
- Calculates boost factor: β γ = p/m for proper lifetime
- Output: Preprocessed CSV with geometry columns (cached)

**Key:** Each HNL processed individually (per-parent counting).

### Stage 3: Limits (Python + HNLCalc)
**Location:** `limits/u2_limit_calculator.py`, `models/hnl_model_hnlcalc.py`

- HNLCalc provides BR(parent→ℓN, |U|²) and cτ₀(mass, |U|²)
- Per-parent counting: `N_sig = Σ_parents [L × σ_parent × BR × ε_geom]`
- Scans 100 |U|² values (10⁻¹² to 10⁻²) to find N_sig = 3 crossings
- Output: Exclusion range [|U|²_min, |U|²_max] at 90% CL

**Key:** Cross-sections from `config/production_xsecs.py` (PBC standard).

---

## Validated Methodology

✅ **Weight handling:** Relative MC weights (not absolute σ)
✅ **Per-parent counting:** Independent σ_D, σ_B, σ_K per species
✅ **HNLCalc integration:** Real BRs and lifetimes
✅ **Geometry:** Ray-tracing with relativistic boosts
✅ **No double-counting:** External cross-sections from literature

**Validation test:** 2.6 GeV muon coupling → |U_mu|² ∈ [6.9×10⁻⁹, 2.4×10⁻⁵]

See [`VALIDATION.md`](VALIDATION.md) for full report.

---

## Directory Structure

```
analysis_pbc_test/
├── README.md                        # This file
├── VALIDATION.md                    # Full validation report
├── HNLCalc/                         # Cloned repository (arXiv:2405.07330)
│   └── HNLCalc.py                  # 150+ production modes, 100+ decays
├── config/
│   └── production_xsecs.py         # LHC cross-sections (σ_K, σ_D, σ_B, σ_W, σ_Z)
├── models/
│   └── hnl_model_hnlcalc.py        # Wrapper: HNLModel(mass, Ue2, Umu2, Utau2)
├── geometry/
│   └── per_parent_efficiency.py    # Ray-tracing and boost calculations
├── limits/
│   ├── u2_limit_calculator.py      # Main analysis driver
│   ├── MULTI_HNL_METHODOLOGY.md    # Per-parent counting explanation
│   └── ROBUSTNESS_FIXES.md         # Defensive programming guide
└── tests/
    ├── test_pipeline.py            # Smoke tests (1.0 GeV muon)
    └── test_26gev_muon.py          # Benchmark (2.6 GeV muon)
```

---

## Example: 2.6 GeV Muon Benchmark

```bash
conda run -n llpatcolliders python tests/test_26gev_muon.py
```

**Output:**
```
Mass:      2.6 GeV
Flavour:   muon (benchmark 010)
Lumi:      3000 fb⁻¹
Peak N_sig: 2.88e+03 events

90% CL Exclusion Range:
  |U_mu|²_min = 6.893e-09
  |U_mu|²_max = 2.364e-05
  Island width: 3.54 decades
```

**Interpretation:**
- **Parent composition:** 86% B-mesons (B⁰/B⁺/Bs), 14% Λb (heavy flavor production)
- **Geometric acceptance:** 1.41% of HNLs reach z = 22m detector
- **Island structure:** Lifetime sweet spot at cτ ~ 10-100m
  - Too short-lived (|U|² > 10⁻⁵): Decays before detector
  - Too long-lived (|U|² < 10⁻⁹): Flies through without decaying

---

## Key Concepts

### Per-Parent Counting

**Why not per-event counting?**

A single pp collision can produce multiple HNLs:
```
Event #44: [B0→N, Bs→N, D+→N, Ds→N]
```

- **Per-event logic:** Count as 1 event with P = 1 - ∏(1-P_i) → **Loses ~50% sensitivity**
- **Per-parent logic:** Count as 4 independent channels → **Physically correct**

Each parent has different production cross-section:
- σ(pp → D⁰) ≈ 2.8 × 10¹⁰ pb
- σ(pp → B⁰) ≈ 4.0 × 10⁸ pb

Cannot assign single σ to multi-parent events!

**See:** `limits/MULTI_HNL_METHODOLOGY.md` for detailed explanation.

### Weight Semantics

Pythia CSV contains `weight` column used for **relative MC reweighting** (not absolute cross-sections).

**Current:** All weights = 1.0 (unweighted generation)

**Validated at:** `production/main_hnl_single.cc:231`
```cpp
double weight = pythia.info.weight();  // ✅ Relative weight
// NOT pythia.info.sigmaGen()          // ❌ Would double-count σ
```

Runtime checks prevent accidental use of absolute σ (w > 10⁶ → error).

---

## Troubleshooting

### Import errors (HNLCalc not found)
```bash
cd analysis_pbc_test
git clone https://github.com/laroccod/HNLCalc.git
conda run -n llpatcolliders pip install sympy mpmath particle numba
```

### Geometry cache errors
Delete corrupted cache files:
```bash
rm ../output/csv/geometry/*.csv
```
Geometry will be recomputed on next run.

### PDG coverage warnings
```
[WARN] Mass 2.60 GeV: 1 parent PDG(s) have no HNLCalc BR: [310]
```

This is expected (KS0 not modeled in HNLCalc). Run diagnostic:
```bash
conda run -n llpatcolliders python limits/diagnostic_pdg_coverage.py
```

---

## References

**Methodology:**
- MATHUSLA: arXiv:1811.00927
- ANUBIS: arXiv:1909.13022
- CODEX-b: arXiv:1911.00481
- PBC Report: CERN-PBC-REPORT-2018-007

**HNL Physics:**
- HNLCalc: arXiv:2405.07330 (Feng, Hewitt, Kling, La Rocco)
- Repository: https://github.com/laroccod/HNLCalc

**Our Implementation:**
- Main documentation: `../CLAUDE.md`
- Validation: `VALIDATION.md`
- Defensive programming: `limits/ROBUSTNESS_FIXES.md`
- Multi-HNL handling: `limits/MULTI_HNL_METHODOLOGY.md`
