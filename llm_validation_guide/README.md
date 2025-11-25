# LLM Validation Guide: HNL Pipeline

**Purpose:** Validate the Heavy Neutral Lepton (HNL) analysis pipeline for the CMS Drainage Gallery detector

**Created:** 2025-11-25

---

## 📊 Pipeline Overview

The HNL analysis has **two independent stages**:

```
┌─────────────────────────────────────────────────────────────┐
│  STAGE 1: PRODUCTION (C++ / Pythia 8.315)                   │
│  Location: production/                                       │
│  Output: CSV files with HNL kinematics                       │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│  STAGE 2: ANALYSIS (Python / HNLCalc)                       │
│  Location: analysis_pbc_test/                               │
│  Output: Exclusion limits |U|²_min, |U|²_max                │
└─────────────────────────────────────────────────────────────┘
```

---

## 🗂️ Validation Structure

This guide is organized into two independent validation sections:

### 📁 `production/`
**Validates:** Pythia event generation (Stage 1)

**Key questions:**
- Are weights handled correctly (relative vs absolute)?
- Is the two-regime production model correct (meson vs EW)?
- Are HNLs kept stable in Pythia?
- Is CSV output format correct?

**Files:** 3 critical files
**Time to review:** ~10 minutes

→ **See:** `production/README_PRODUCTION.md`

---

### 📁 `analysis/`
**Validates:** Limit calculation pipeline (Stage 2)

**Key questions:**
- Is per-parent counting implemented correctly?
- Are cross-sections from external sources (not double-counted)?
- Is geometry/lifetime calculation correct?
- Do tests confirm expected results?

**Files:** 7 critical files
**Time to review:** ~20-30 minutes

→ **See:** `analysis/README_ANALYSIS.md`

---

## 🚀 Quick Start

### For a Complete Review:
```bash
cd llm_validation_guide

# 1. Read overview
cat CLAUDE.md  # Full physics context

# 2. Validate production
cd production
cat README_PRODUCTION.md

# 3. Validate analysis
cd ../analysis
cat README_ANALYSIS.md
```

---

### For Express Validation (5 minutes):

**Three most critical files to check:**

1. **`production/main_hnl_single.cc`** (line ~231)
   - ✅ Must use `.weight()` not `.sigmaGen()`
   - ❌ If wrong: 10⁶× overcounting

2. **`analysis/MULTI_HNL_METHODOLOGY.md`**
   - ✅ Must use per-parent counting
   - ❌ If wrong: ~50% sensitivity loss

3. **`analysis/production_xsecs.py`**
   - ✅ Cross-sections from literature
   - ❌ If wrong: Double-counting

---

## 📖 Files in This Directory

```
llm_validation_guide/
├── README.md                    ← You are here
├── CLAUDE.md                    ← Physics overview (symlink)
│
├── production/                  ← Stage 1: Event generation
│   ├── README_PRODUCTION.md
│   ├── main_hnl_single.cc
│   ├── make.sh
│   └── hnl_*_Template.cmnd
│
└── analysis/                    ← Stage 2: Limit calculation
    ├── README_ANALYSIS.md
    ├── MULTI_HNL_METHODOLOGY.md
    ├── VALIDATION.md
    ├── u2_limit_calculator.py
    ├── per_parent_efficiency.py
    ├── hnl_model_hnlcalc.py
    ├── production_xsecs.py
    ├── test_26gev_muon.py
    └── test_expected_signal_events_kernel.py
```

---

## 🎯 Validation Workflow

### Step 1: Production Validation
Go to `production/` folder and verify:
- [ ] Pythia configuration is correct
- [ ] Weight handling is relative (not absolute)
- [ ] Two-regime model (meson/EW) is implemented
- [ ] CSV output format is correct

### Step 2: Analysis Validation
Go to `analysis/` folder and verify:
- [ ] Per-parent counting methodology
- [ ] Cross-section normalization
- [ ] Geometry and boost calculations
- [ ] Tests pass with expected results

### Step 3: Integration Check
- [ ] Production CSV format matches analysis expectations
- [ ] Parent PDG codes are consistent
- [ ] Weight column semantics are documented

---

## 🔍 Critical Decision Points

### Production (Stage 1)
**Decision:** How to handle event weights?

**Correct:**
```cpp
double weight = pythia.info.weight();  // Relative MC weight
```

**Wrong:**
```cpp
double weight = pythia.info.sigmaGen();  // Absolute σ → double-counts!
```

**Impact if wrong:** Results off by 10⁶×

---

### Analysis (Stage 2)
**Decision:** Per-parent or per-event counting?

**Correct:**
```python
# Per-parent: Each HNL contributes to its parent's cross-section
N_sig = Σ_parents [L × σ_parent × BR(parent→ℓN) × ε_geom(parent)]
```

**Wrong:**
```python
# Per-event: Combines multiple HNLs from one collision
N_sig = L × σ_??? × P_event  # Cannot assign single σ!
```

**Impact if wrong:** Lose ~50% sensitivity

---

## 📚 Documentation Structure

### Top Level (This Directory)
- **`README.md`**: This file - navigation guide
- **`CLAUDE.md`**: Complete physics overview and pipeline description

### Production Subfolder
- Focused on Pythia event generation
- C++ code validation
- Configuration file review

### Analysis Subfolder
- Focused on Python limit calculation
- Methodology validation (per-parent counting)
- Test suite verification

---

## 🧪 Running Tests

### Production Tests
```bash
cd ../production

# Test single mass point
./main_hnl_single 2.6 muon

# Check output
ls -lh ../output/csv/simulation/HNL_mass_2p6_muon*.csv
head ../output/csv/simulation/HNL_mass_2p6_muon_Meson.csv
```

### Analysis Tests
```bash
cd ../analysis_pbc_test

# Closure test (pure math)
conda run -n llpatcolliders python tests/closure_anubis/test_expected_signal_events_kernel.py

# Benchmark test (end-to-end)
conda run -n llpatcolliders python tests/test_26gev_muon.py
```

---

## ✅ Validation Checklist

### Production Stage
- [ ] Weight handling verified (`production/main_hnl_single.cc:~231`)
- [ ] Two-regime logic correct (m<5 GeV: meson, m≥5 GeV: EW)
- [ ] HNL is stable in Pythia (`mayDecay = off`)
- [ ] CSV format includes: event, weight, parent_id, kinematics

### Analysis Stage
- [ ] Per-parent counting confirmed (`analysis/MULTI_HNL_METHODOLOGY.md`)
- [ ] External cross-sections used (`analysis/production_xsecs.py`)
- [ ] HNLCalc integration correct (`analysis/hnl_model_hnlcalc.py`)
- [ ] Geometry/boost calculation correct (`analysis/per_parent_efficiency.py`)
- [ ] Tests pass (`analysis/test_*.py`)

### Integration
- [ ] Production CSV → Analysis input works seamlessly
- [ ] Parent PDG codes consistent between stages
- [ ] Weight semantics documented and correct

---

## 🚨 Red Flags

### Production
- ⛔ Weights > 10⁶ in CSV (indicates absolute σ usage)
- ⛔ HNL decays in Pythia (should be stable)
- ⛔ Missing parent_id column in CSV

### Analysis
- ⛔ Per-event probability formula `1 - ∏(1-P_i)`
- ⛔ Cross-sections from `pythia.info.sigmaGen()`
- ⛔ Tests fail or produce unrealistic results

---

## 📖 References

**LLP Detector Methodology:**
- MATHUSLA: arXiv:1811.00927
- ANUBIS: arXiv:1909.13022
- CODEX-b: arXiv:1911.00481

**HNL Physics:**
- HNLCalc: arXiv:2405.07330

**This Project:**
- Complete guide: `CLAUDE.md`
- Production: `production/README_PRODUCTION.md`
- Analysis: `analysis/README_ANALYSIS.md`

---

## 💡 For LLM Reviewers

**Recommended reading order:**

1. **This file** - Understand the two-stage structure
2. **`CLAUDE.md`** - Physics context (optional but helpful)
3. **`production/README_PRODUCTION.md`** - Validate Stage 1
4. **`analysis/README_ANALYSIS.md`** - Validate Stage 2

**Time investment:**
- Express review: 5 minutes (3 critical files)
- Production only: 10 minutes (3 files)
- Analysis only: 20 minutes (7 files)
- Complete review: 30 minutes (all files)

---

**Created by:** Claude Code (Anthropic)
**Date:** 2025-11-25
**Repository:** llpatcolliders/llpatcolliders
