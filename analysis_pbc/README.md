# PBC-Grade HNL Analysis Pipeline

Physics Beyond Colliders (PBC) compatible analysis for the CMS drainage gallery long-lived particle detector.

## Quick Start

### 1. Installation

```bash
# Ensure conda environment exists
conda activate llpatcolliders

# Install HNLCalc (if not already installed)
cd analysis_pbc
git clone https://github.com/laroccod/HNLCalc.git
conda run -n llpatcolliders pip install sympy mpmath particle numba
conda run -n llpatcolliders pip install 'scikit-hep==0.4.0'
```

### 2. Run Pipeline Tests

```bash
# From llpatcolliders/analysis_pbc/
conda run -n llpatcolliders python tests/closure_anubis/test_expected_signal_events_kernel.py
conda run -n llpatcolliders python tests/test_26gev_muon.py
```

### 3. Combine Production Channels

**IMPORTANT:** Before running analysis, you must combine MadGraph (EW) and Pythia (meson) production files to avoid double-counting:

```bash
# Combine all overlapping production channels
conda run -n llpatcolliders python limits/combine_production_channels.py

# Expected output:
# - Creates combined files for masses with multiple production regimes
# - Deletes original separate files (data preserved in combined files)
# - Saves ~2 GB of disk space by removing duplicates
```

**What this does:**
- At overlapping masses (e.g., 2-6 GeV), HNLs are produced from BOTH:
  - Meson decays (B/D/K → ℓN) via Pythia
  - EW decays (W/Z → ℓN) via MadGraph
- These are different production mechanisms that must be ADDED (not double-counted)
- The script combines CSV files at each mass into a single unified file

**Why this is critical:**
- Without combining: Analysis would only use one production channel per mass
- After combining: Analysis correctly includes all production mechanisms
- Disk space: Original files are deleted after combining (data preserved)

### 4. Calculate |U|² Limits

```bash
# For all mass points, all flavors (parallel processing)
conda run -n llpatcolliders python limits/run_serial.py --parallel

# Or serial processing (slower)
conda run -n llpatcolliders python limits/run_serial.py

# Or run benchmark test (2.6 GeV muon)
conda run -n llpatcolliders python tests/test_26gev_muon.py
```

---

## Pipeline Overview

The analysis follows a **three-stage** workflow matching PBC methodology (MATHUSLA/CODEX-b/ANUBIS):

### Stage 1: Production (Pythia C++)
**Location:** `../production/pythia_production/main_hnl_production.cc`

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
**Location:** `limits/run_serial.py`, `limits/expected_signal.py`, `models/hnl_model_hnlcalc.py`

- HNLCalc provides BR(parent→ℓN, |U|²) and cτ₀(mass, |U|²)
- Per-parent counting: `N_sig = Σ_parents [L × σ_parent × BR × ε_geom]`
- Scans 100 |U|² values (10⁻¹² to 10⁻²) to find N_sig = 3 crossings
- Output: Exclusion range [|U|²_min, |U|²_max] at 95% CL

**Key:** Cross-sections from `config/production_xsecs.py` (PBC standard).

---

## Validated Methodology

✅ **Weight handling:** Relative MC weights (not absolute σ)
✅ **Per-parent counting:** Independent σ_D, σ_B, σ_K per species
✅ **HNLCalc integration:** Real BRs and lifetimes
✅ **Geometry:** Ray-tracing with relativistic boosts
✅ **No double-counting:** External cross-sections from literature

**Validation test (repo sample):** 2.6 GeV muon coupling → |U_mu|² ∈ [5.5×10⁻⁹, 9.5×10⁻⁵]

See [`VALIDATION.md`](VALIDATION.md) for full report.

---

## Directory Structure

```
analysis_pbc/
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
│   ├── combine_production_channels.py  # Combine MadGraph + Pythia (run FIRST!)
│   ├── run_serial.py                   # Main analysis driver (parallel/serial)
│   ├── expected_signal.py              # Signal-yield kernel (expected_signal_events)
│   ├── MULTI_HNL_METHODOLOGY.md        # Per-parent counting explanation
│   └── ROBUSTNESS_FIXES.md             # Defensive programming guide
└── tests/
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
Peak N_sig: 2.81e+05 events

95% CL Exclusion Range:
  |U_mu|²_min = 5.462e-09
  |U_mu|²_max = 9.545e-05
  Island width: 4.24 decades
```

**Interpretation:**
- **Parent composition:** mixed heavy-flavor + EW (`B` mesons plus `W/Z` in the combined sample)
- **Geometric acceptance:** ~1–2% of HNLs hit the detector (sample-dependent)
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

**Validated at:** `production/pythia_production/main_hnl_production.cc`
```cpp
double weight = pythia.info.weight();  // ✅ Relative weight
// NOT pythia.info.sigmaGen()          // ❌ Would double-count σ
```

If you ever see extremely large `weight` values, it likely means absolute σ was written into the CSV (which would double-count σ in the analysis).

---

## Troubleshooting

### Import errors (HNLCalc not found)
```bash
cd analysis_pbc
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

This is expected (KS0 not modeled in HNLCalc). These events are dropped; impact is typically negligible.

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
