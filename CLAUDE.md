# HNL LLP Detector Pipeline - Quick Reference

**Project:** Heavy Neutral Lepton (HNL) search at CMS Drainage Gallery Detector
**Physics:** Sterile neutrino displaced decays at LHC (√s = 14 TeV)
**Mass Range:** 0.2 - 80 GeV | **Luminosity:** 3000 fb⁻¹ (HL-LHC)

---

## Quick Start (3 Commands)

```bash
# 1. Generate events (Pythia meson + MadGraph EW)
cd production/pythia_production && ./run_full_production.sh
cd production/madgraph_production && python3 scripts/run_hnl_scan.py

# 2. Calculate limits
cd analysis_pbc && conda run -n llpatcolliders python limits/run_serial.py

# 3. Generate plot
cd money_plot && conda run -n llpatcolliders python plot_money_island.py
```

**Output:** `output/images/HNL_moneyplot_island.png`

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
│   ├── geometry/                 # Ray-tracing detector mesh
│   ├── models/                   # HNLCalc wrapper
│   └── limits/run_serial.py      # Main analysis driver
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
N_sig = L × Σ_parents [σ_parent × BR(parent→ℓN) × ε_geom(parent)]
```

Where:
- **L** = 3000 fb⁻¹ (HL-LHC)
- **σ_parent** = Production cross-section (from `config/production_xsecs.py`)
- **BR(parent→ℓN)** = Branching ratio from HNLCalc (∝ |U_ℓ|²)
- **ε_geom** = Geometric acceptance (ray-tracing + decay probability)

**Exclusion threshold:** N_sig ≥ 3 (95% CL)

---

## Mass Grids

**Configuration:** `config_mass_grid.py` (single source of truth)

```bash
python config_mass_grid.py  # View all grids
```

**Production grids (base + closure):**
- Electron: 50 points (0.2-5.5 GeV)
- Muon: 50 points (0.2-5.5 GeV)
- Tau: 37 points (0.5-4.4 GeV)

**EW grids (MadGraph):**
- Electron: 15 points (5-80 GeV)
- Muon: 13 points (5-80 GeV)
- Tau: 11 points (6-80 GeV)

---

## CSV File Format

**Pythia/MadGraph Output:**
```csv
event,weight,hnl_id,parent_pdg,pt,eta,phi,p,E,mass,prod_x_mm,prod_y_mm,prod_z_mm,boost_gamma
```

**Key columns for analysis:**
- `parent_pdg`: Parent particle (511=B⁰, 521=B±, 24=W±, etc.)
- `p`, `mass`: For boost factor β γ = p/m
- `eta`, `phi`: For ray-tracing trajectory
- `weight`: Relative MC weight (typically 1.0)

**Location:** `output/csv/simulation_new/`

---

## Analysis Pipeline

### Stage 1: Geometry Preprocessing (Automatic, Cached)

```python
# analysis_pbc/geometry/per_parent_efficiency.py
for each HNL:
    - Ray-trace trajectory from IP to detector
    - Compute entry distance and path length
    - Calculate β γ = p/m
```

**Cache:** `output/csv/geometry/` (30x speedup on subsequent runs)

### Stage 2: Limit Calculation

```python
# analysis_pbc/limits/run_serial.py
for mass, flavour, |U|² in scan:
    N_sig = Σ_parents [σ × BR(|U|²) × ε_geom]
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

# Level 1: Integration (~30 sec)
python tests/test_pipeline.py

# Level 2: Physics benchmark (~3 min)
python tests/test_26gev_muon.py
```

**Validation:** 2.6 GeV muon benchmark
- **Expected:** |Uμ|² ∈ [6.9×10⁻⁹, 2.4×10⁻⁵]
- **Peak signal:** 2880 events at |Uμ|² ~ 1×10⁻⁶

---

## Common Issues

### 1. Missing HNLCalc
```bash
cd analysis_pbc
git clone https://github.com/laroccod/HNLCalc.git
pip install sympy mpmath particle numba 'scikit-hep==0.4.0'
```

### 2. Empty CSV Files
- Check logs in `output/logs/simulation_new/`
- High-mass Pythia (m ≥ 5 GeV) fails → use MadGraph instead

### 3. Geometry Cache Corrupt
```bash
rm output/csv/geometry/HNL_*_geom.csv  # Force recompute
```

### 4. Wrong Parent PDGs
- Check `limits/diagnostic_pdg_coverage.py` for missing BRs
- PDG 310 (K_S⁰) not in HNLCalc → <0.1% loss

---

## File Naming Conventions

**Simulation CSVs:**
```
HNL_{mass}GeV_{flavour}_{regime}[_mode].csv

Examples:
  HNL_2p60GeV_muon_beauty.csv
  HNL_10p0GeV_electron_ew.csv
  HNL_1p00GeV_tau_charm_fromTau.csv
```

**Regimes:**
- `kaon` (m < 0.5 GeV)
- `charm` (0.5-2 GeV)
- `beauty` (2-5 GeV)
- `ew` (≥ 5 GeV)

**Tau modes:**
- `direct`: B/D/W → τN (all masses)
- `fromTau`: τ → πN cascade (m < 1.64 GeV only)

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
ls output/csv/simulation_new/*.csv | wc -l  # Expected: ~180 files
find output/csv/simulation_new/ -name "*.csv" -size 0  # No empty files

# Verify parent PDGs
awk -F',' 'NR>1 {print $4}' HNL_2p60GeV_muon_beauty.csv | sort | uniq -c
# Should see: 511 (B⁰), 521 (B±), 531 (Bs)
```

---

## Reference Materials

**Code:**
- `production/pythia_production/main_hnl_production.cc` - Pythia simulation
- `production/madgraph_production/scripts/run_hnl_scan.py` - MadGraph driver
- `analysis_pbc/limits/u2_limit_calculator.py` - Core analysis logic
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

**Last Updated:** December 2024
**Status:** ✅ Production-ready (meson + EW regimes)
