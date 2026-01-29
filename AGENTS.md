# HNL LLP Detector Pipeline - Quick Reference

**Project:** Heavy Neutral Lepton (HNL) search at CMS Drainage Gallery Detector
**Physics:** Sterile neutrino displaced decays at LHC (√s = 14 TeV)
**Mass Range:** 0.2 - 80 GeV | **Luminosity:** 3000 fb⁻¹ (HL-LHC)

Policy: From tag hnl-pre-decay-final onward, all HNL analyses require explicit decay simulation; geometry-only acceptance is deprecated and removed.

---

## Quick Start (4 Commands)

```bash
# 1. Generate events (Pythia meson + MadGraph EW)
cd production/pythia_production && ./run_full_production.sh          # conda env required
cd production/madgraph_production && python3 scripts/run_hnl_scan.py # uses MG5/madgraph env (NOT conda)

# 2. Combine production channels ⚠️ REQUIRED!
cd analysis_pbc && conda run -n llpatcolliders python limits/combine_production_channels.py

# 3. Calculate limits
cd analysis_pbc && conda run -n llpatcolliders python limits/run.py --parallel

# 4. Generate plot
mkdir -p output/images && cd money_plot && conda run -n llpatcolliders python plot_money_island.py
```

**Output:** `output/images/HNL_moneyplot_island.png`

**⚠️ Critical:** Step 2 (combine) is REQUIRED to avoid undercounting signal at overlapping masses!

---

## Environment (Important)

- **Default:** All analysis/plotting scripts **must** run inside the `llpatcolliders` conda env.
  - Interactive shell: `conda activate llpatcolliders`
  - One-off command: `conda run -n llpatcolliders python <script>.py`
- **Exception:** `production/madgraph_production` uses its own MG5 setup; keep using `python3 scripts/run_hnl_scan.py` there (do *not* switch to the conda env).
- If Matplotlib cache complains about permissions, set `MPLCONFIGDIR` to a writable path, e.g. `export MPLCONFIGDIR=$PWD/.mplcache`.
- **HNLCalc caches:** Cache paths are anchored to `analysis_pbc/model/` (generated + gitignored) so they won't pollute the repo root. Still prefer running from inside `analysis_pbc` to keep outputs tidy and consistent.

---

## Project Structure

```
├── config_mass_grid.py           # SINGLE SOURCE OF TRUTH for mass grids
├── production/
│   ├── pythia_production/        # Meson production (m < 5 GeV)
│   │   ├── main_hnl_production.cc
│   │   └── run_full_production.sh
│   └── madgraph_production/      # EW production (m ≥ 5 GeV)
│       └── scripts/run_hnl_scan.py
├── analysis_pbc/                 # Limit calculation pipeline
│   ├── config/                   # Cross-sections (σ_K, σ_D, σ_B, σ_W, σ_Z)
│   ├── decay/                    # HNL decay simulation (REQUIRED)
│   │   ├── decay_detector.py     # Decay sampling + track separation
│   │   ├── rhn_decay_library.py  # MATHUSLA decay file loader
│   │   └── external/             # MATHUSLA_LLPfiles_RHN_U{e,mu,tau}
│   ├── geometry/                 # Ray-tracing detector mesh
│   ├── models/                   # HNLCalc wrapper
│   ├── cache/                    # Geometry + decay cache files
│   └── limits/
│       ├── combine_production_channels.py  # Combine meson + EW (RUN FIRST!)
│       ├── run.py                           # Main analysis driver
│       └── expected_signal.py               # Signal yield kernel
└── money_plot/                   # Exclusion plots
```

---

## Physics Essentials

### HNL Production

| Mass | Channel | Generator | Cross-section |
|------|---------|-----------|---------------|
| 0.2-0.5 GeV | K → ℓN | Pythia | σ(K⁺) ≈ 5×10¹⁰ pb |
| 0.5-2 GeV | D → ℓN | Pythia | σ(D⁰) ≈ 3×10¹⁰ pb |
| 2-5 GeV | B → ℓN | Pythia | σ(B⁰) ≈ 4×10⁸ pb |
| 5-80 GeV | W/Z → ℓN | MadGraph | σ(W) ≈ 2×10⁸ pb |

### Detector

- **Location:** CMS drainage gallery, z = 22m above IP
- **Volume:** ~800 m³ tube (radius 1.54m, length ~100m)
- **Sensitivity:** HNL lifetimes cτ ~ 10-100m

### Signal Formula

```
N_sig = L × Σ_parents [σ_parent × BR(parent→ℓN) × P_decay × ε_sep]
```

Where:
- **L** = 3000 fb⁻¹ (HL-LHC)
- **σ_parent** = Production cross-section (from `config/production_xsecs.py`)
- **BR(parent→ℓN)** = Branching ratio from HNLCalc (∝ |U_ℓ|²)
- **P_decay** = Probability to decay inside detector volume (depends on cτ ∝ 1/|U|²)
- **ε_sep** = Track separation efficiency (charged decay products ≥ 1mm apart)

**Exclusion threshold:** N_sig ≥ 3 (95% CL)

---

## Mass Grids

**Configuration:** `config_mass_grid.py` (single source of truth)

```bash
python config_mass_grid.py  # View all grids
```

**Production grids (base + closure):**
- Electron: 93 points (0.2-17 GeV common grid + tau fine low-mass points)
- Muon: 93 points (0.2-17 GeV common grid + tau fine low-mass points)
- Tau: 44 points with **finer grid** (0.03-0.05 GeV steps in 0.2-1.6 GeV region)
  - Finer spacing reduces jaggedness in exclusion curve
  - fromTau mode uses 1M events (10x more than direct) for better statistics

**EW grids (MadGraph):**
- Electron: 93 points (unified grid, currently up to 17 GeV)
- Muon: 93 points (unified grid, currently up to 17 GeV)
- Tau: 93 points (unified grid, currently up to 17 GeV)

---

## CSV File Format

**Pythia/MadGraph Output:**
```csv
event,weight,hnl_id,parent_pdg,pt,eta,phi,p,E,mass,prod_x_mm,prod_y_mm,prod_z_mm,beta_gamma
```

**Key columns for analysis:**
- `parent_pdg`: Parent particle (511=B⁰, 521=B±, 24=W±, etc.)
- `p`, `mass`: For boost factor β γ = p/m
- `eta`, `phi`: For ray-tracing trajectory
- `weight`: Relative MC weight (typically 1.0)

**Location:** `output/csv/simulation/`

---

## Analysis Pipeline

### Stage 0: Combine Production Channels (Required!)

**Why:** At overlapping masses (2-6 GeV), HNLs are produced from BOTH meson and EW decays:
- Meson production: B/D/K → ℓN (Pythia)
- EW production: W/Z → ℓN (MadGraph)

These are **different physics channels** that must be ADDED, not double-counted.

```bash
cd analysis_pbc
conda run -n llpatcolliders python limits/combine_production_channels.py
```

**What it does:**
- Finds masses with multiple production files (e.g., `*_beauty.csv` + `*_ew.csv`)
- Combines them into `*_combined.csv`
- Deletes original separate files (saves ~2 GB)
- Preserves all parent PDG codes for proper per-parent counting

**Example:**
```
Before:  HNL_4p00GeV_muon_beauty.csv  (86,615 HNLs from B mesons)
         HNL_4p00GeV_muon_ew.csv      (100,000 HNLs from W/Z)
After:   HNL_4p00GeV_muon_combined.csv (186,615 HNLs total)
```

**Critical:** Without combining, analysis would only see one channel → undercounts signal by ~50%!

### Stage 1: Geometry Preprocessing (Automatic, Cached)

```python
# analysis_pbc/geometry/per_parent_efficiency.py
for each HNL:
    - Ray-trace trajectory from IP to detector
    - Compute entry distance and path length
    - Calculate β γ = p/m
```

**Cache:** `analysis_pbc/cache/geom_{flavour}_{mass}.csv` (fast reload on subsequent runs)

### Stage 2: Decay + Track Separation (REQUIRED)

```python
# analysis_pbc/decay/decay_detector.py
for each HNL hitting detector:
    - Sample decay event from MATHUSLA RHN files
    - Boost charged daughters to lab frame
    - Ray-trace to detector surface
    - Check minimum track separation (default: 1mm)
```

**Cache:** `analysis_pbc/cache/decaycache_{flavour}_m{mass}_sep{sep}mm_seed{seed}.pkl`

**Key insight:** Separation cut is precomputed once (at path midpoint) and reused across all |U|² values.
This approximation is valid because track separation depends mainly on decay angles, not exact position.

**Decay file selection (flavour-aware):**
- Electron/muon: prioritize `inclDs`, `inclDD`, `inclD`, `nocharm`, `nocharmnoss`, `lightfonly`; use analytical 2/3-body files below the low-mass threshold.
- Tau: use analytical 2/3-body files below 0.42 GeV; above that, select the nearest `lightfonly` / `lightfstau` / `lightfstauK` file.

### Stage 3: Limit Calculation

```python
# analysis_pbc/limits/run.py --parallel --separation-mm 1.0
for mass, flavour, |U|² in scan:
    N_sig = Σ_parents [σ × BR(|U|²) × P_decay × ε_sep]
    if N_sig ≥ 3: excluded
```

**Output:** `output/csv/analysis/HNL_U2_limits_summary.csv`

### Stage 3: Plotting

```python
# money_plot/plot_money_island.py
- Read limits CSV
- Create 3-panel plot (electron/muon/tau)
- Shade excluded regions ("islands")
```

**Output:** `output/images/HNL_moneyplot_island.png`

---

## Key Implementation Details

### Per-Parent Counting (Critical!)

Pythia events can have **multiple HNLs from different parents**:
```csv
event,parent_pdg
44,511        # B⁰
44,521        # B⁺
44,431        # Ds
```

**Count separately** (not per-event) because each parent has different σ.

### Weight Handling

- CSV `weight` = **relative MC weight** (typically 1.0)
- **NOT absolute cross-section** (would double-count)
- Real σ applied from `config/production_xsecs.py`

### HNLCalc Integration

```python
from models.hnl_model_hnlcalc import HNLModel
model = HNLModel(mass_GeV=2.6, Ue2=0.0, Umu2=1e-6, Utau2=0.0)

# Physics outputs:
ctau0_m = model.ctau0_m              # Proper lifetime [m]
brs = model.production_brs()         # {parent_pdg: BR(parent→ℓN)}
```

**Source:** arXiv:2405.07330 (150+ production, 100+ decay channels)

---

## Tests

```bash
cd analysis_pbc

# Level 0: Math kernel (~1 sec)
python tests/closure_anubis/test_expected_signal_events_kernel.py

# Level 1: Fast single mass point (~4 min with caching)
PYTHONUNBUFFERED=1 python test_fast_muon_1gev.py

# Level 2: Physics benchmark (~5 min)
python tests/test_26gev_muon.py
```

**Validated benchmarks:**

| Mass | Flavour | Peak N_sig | Exclusion |U|² | Width |
|------|---------|------------|-----------|-------|
| 1.0 GeV | muon | 3.2×10⁸ | [4.3×10⁻⁹, 1×10⁻²] | 6.4 dec |
| 2.6 GeV | muon | 2.3×10⁵ | [7.5×10⁻⁹, 9.1×10⁻⁵] | 4.1 dec |

**Physics notes:**
- 1.0 GeV dominated by D mesons (σ_charm >> σ_bottom), longer cτ → wider island
- 2.6 GeV dominated by B mesons, shorter cτ → narrower island

---

## Common Issues

### 1. Forgot to Combine Production Channels ← MOST COMMON!
**Symptom:** Analysis results look weak, missing mass points

**Solution:**
```bash
cd analysis_pbc
conda run -n llpatcolliders python limits/combine_production_channels.py
```

**Check:** After combining, you should see `*_combined.csv` files:
```bash
ls output/csv/simulation/*_combined.csv | wc -l  # Should be ~122
```

### 2. Missing HNLCalc
```bash
cd analysis_pbc
git clone https://github.com/laroccod/HNLCalc.git
touch HNLCalc/__init__.py  # Make it a Python package
pip install sympy mpmath particle numba 'scikit-hep==0.4.0'
```

### 3. Empty CSV Files
- Check logs in `output/logs/simulation/`
- High-mass Pythia (m ≥ 5 GeV) fails → use MadGraph instead

### 4. Geometry Cache Corrupt
```bash
rm output/csv/geometry/HNL_*_geom.csv  # Force recompute
```

### 5. Wrong Parent PDGs
- Run `limits/run.py` and watch for `[WARN] ... have no HNLCalc BR` / `no cross-section` messages from `limits/expected_signal.py`
- K_L (130) now supported; K_S (310) omitted (τ_S/τ_L ≈ 1/570 suppression)

---

## File Naming Conventions

**Simulation CSVs (after combining):**
```
HNL_{mass}GeV_{flavour}_{regime}[_mode].csv

Examples (post-combine):
  HNL_2p60GeV_muon_combined.csv    ← beauty + ew merged
  HNL_10p0GeV_electron_ew.csv      ← EW only (no overlap)
  HNL_0p50GeV_tau_combined.csv     ← kaon + charm + ew merged
```

**Regimes:**
- `kaon` (m < 0.5 GeV)
- `charm` (0.5-2 GeV)
- `beauty` (2-5 GeV)
- `ew` (≥ 5 GeV)
- `combined` (multiple regimes merged)

**Tau modes:**
- `direct`: B/D/W → τN (all masses)
- `fromTau`: τ → πN cascade (m < 1.64 GeV only)

**⚠️ Tau statistics note:** The fromTau mode produces very few HNLs (~0.4% efficiency).
To get adequate statistics, `run_parallel_production.sh` uses `NEVENTS_FROMTAU=1000000` (10x more events).
Even with 1M events, only ~1-2 fromTau HNLs pass the full selection chain (geometry + decay with ≥2 charged tracks),
causing some residual jaggedness in the tau exclusion curve below 1.6 GeV.
Additional small wiggles can come from discrete decay-file mass points and the 100-step |U|² scan grid.

**Note:** After running `combine_production_channels.py`, overlapping masses will have `*_combined.csv` files instead of separate regime files.

---

## Production Tips

### Parallel Execution
```bash
# 10 cores → ~2 hours (vs 20 hours single-core)
MAX_PARALLEL=10
for mass in masses; do
    wait_for_slot
    ./main_hnl_production $mass $flavour 200000 &
done
wait
```

### Validation
```bash
# Check output
ls output/csv/simulation/*.csv | wc -l  # Expected: ~180 files
find output/csv/simulation/ -name "*.csv" -size 0  # No empty files

# Verify parent PDGs
awk -F',' 'NR>1 {print $4}' output/csv/simulation/HNL_2p60GeV_muon_combined.csv | sort | uniq -c | head
# Should include e.g. 511 (B⁰), 521 (B±), 531 (Bs) (and possibly 24/23 if EW is present).
```

---

## Reference Materials

**Code:**
- `production/pythia_production/main_hnl_production.cc` - Pythia simulation
- `production/madgraph_production/scripts/run_hnl_scan.py` - MadGraph driver
- `analysis_pbc/limits/combine_production_channels.py` - Combine meson + EW (REQUIRED!)
- `analysis_pbc/limits/run.py` - Main analysis driver
- `analysis_pbc/limits/expected_signal.py` - Signal-yield kernel (expected_signal_events)
- `config_mass_grid.py` - Mass grid definitions

**Physics:**
- HNLCalc: arXiv:2405.07330
- PBC Report: CERN-PBC-REPORT-2018-007
- MATHUSLA: arXiv:1811.00927
- ANUBIS: arXiv:1909.13022

**Validation:**
- `analysis_pbc/VALIDATION.md` - Methodology validation
- `analysis_pbc/README.md` - Analysis details

---

**Last Updated:** January 2025 (tau finer grid + fromTau 1M events)
**Status:** ✅ Production-ready (decay simulation + track separation cuts)
