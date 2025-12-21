# LLP at Colliders - Complete Documentation

## Purpose

This repository calculates sensitivity limits for **Long-Lived Particles (LLPs)** at the LHC, targeting far detectors like the CMS drainage gallery. It currently implements **Heavy Neutral Leptons (HNLs)** and is designed to be extended to **Axion-Like Particles (ALPs)**.

The final output is a "**money plot**" showing exclusion regions in the (mass, coupling) parameter space.

---

## Repository Structure

```
llpatcolliders/
├── config_mass_grid.py                 # Central mass grid configuration
├── LLPATCOLLIDERS_DOCUMENTATION.md     # This file
│
├── production/
│   ├── pythia_production/              # Stage 1A: Meson production
│   │   ├── main_hnl_production.cc      # C++ Pythia8 generator
│   │   ├── run_parallel_production.sh  # Parallel job orchestrator
│   │   ├── load_mass_grid.sh           # Mass grid loader
│   │   └── pythia8315/                 # Pythia 8.315 library
│   │
│   └── madgraph_production/            # Stage 1B: Electroweak production
│       ├── Dockerfile                  # Docker container setup
│       ├── scripts/
│       │   ├── run_hnl_scan.py         # Main MadGraph driver
│       │   └── lhe_to_csv.py           # LHE to CSV converter
│       ├── cards/                      # MadGraph process/run cards
│       └── mg5/                        # MadGraph installation
│
├── analysis_pbc/                       # Stage 2-3: Analysis pipeline
│   ├── HNLCalc/                        # External: decay widths/BRs
│   │   └── HNLCalc.py                  # 150+ decay modes
│   ├── config/
│   │   └── production_xsecs.py         # Cross-sections database
│   ├── models/
│   │   └── hnl_model_hnlcalc.py        # HNL physics model wrapper
│   ├── geometry/
│   │   └── per_parent_efficiency.py    # Ray-tracing detector geometry
│   └── limits/
│       ├── run.py                      # Main analysis driver
│       ├── expected_signal.py          # Physics kernel (N_sig calculation)
│       └── combine_production_channels.py  # Merge production files
│
├── money_plot/                         # Stage 4: Visualization
│   └── plot_money_island.py            # Final exclusion plot
│
└── output/
    └── csv/
        ├── simulation/                 # Production CSV outputs
        ├── geometry/                   # Cached geometry calculations
        └── analysis/                   # Final limits summary
```

---

## Complete Pipeline Overview

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                         STAGE 1: PARTICLE PRODUCTION                        │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│   PYTHIA (Meson Decays)                  MADGRAPH (Electroweak)             │
│   ─────────────────────                  ───────────────────────            │
│   K± → ℓ± N                              W± → ℓ± N                          │
│   D+/D0/Ds → ℓ N (+ X)                   Z  → ν N                           │
│   B+/B0/Bs → ℓ N (+ X)                                                      │
│   τ → π N (cascade)                                                         │
│                                                                             │
│   Output: HNL kinematics (pT, η, φ, production vertex, boost factor)        │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
                                     │
                                     ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                        STAGE 2: GEOMETRY RAY-TRACING                        │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│   For each HNL trajectory:                                                  │
│   1. Convert (η, φ) → 3D direction vector                                   │
│   2. Ray-trace from production vertex through detector mesh                 │
│   3. Record: hits_tube, entry_distance, path_length                         │
│                                                                             │
│   Detector: CMS drainage gallery at z = 22m, cylindrical tube               │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
                                     │
                                     ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                        STAGE 3: LIMITS CALCULATION                          │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│   For each mass point:                                                      │
│   1. Scan coupling |U|² from 10⁻¹² to 10⁻²                                  │
│   2. Calculate N_sig = Σ [L × σ_parent × BR × P_decay × ε_geom]             │
│   3. Find exclusion range where N_sig ≥ 2.996 (3σ threshold)                │
│                                                                             │
│   Uses HNLCalc for decay widths/branching ratios                            │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
                                     │
                                     ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                         STAGE 4: MONEY PLOT                                 │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│   Plot exclusion "island" in (mass, |U|²) plane:                            │
│   - Lower boundary: HNL too long-lived (decays downstream)                  │
│   - Upper boundary: HNL too prompt (decays before detector)                 │
│   - Shaded region: Excluded parameter space                                 │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## Stage 1A: Pythia Production (Meson Decays)

### Purpose
Generate HNL kinematics from meson decays at 14 TeV LHC using Pythia 8.315.

### Source File
`production/pythia_production/main_hnl_production.cc`

### Physics Processes

| Parent | Decay Mode | Kinematic Limit (m_N) |
|--------|------------|----------------------|
| K± | K± → ℓ± N | e: <0.494 GeV, μ: <0.387 GeV |
| D± | D± → ℓ± N | e: <1.870 GeV, μ: <1.764 GeV |
| D0 | D0 → K ℓ N | e: <1.865 GeV, μ: <1.759 GeV |
| Ds | Ds → ℓ N | e: <1.968 GeV, μ: <1.862 GeV |
| B± | B± → ℓ± N | e: <5.279 GeV, μ: <5.173 GeV, τ: <3.502 GeV |
| B0 | B0 → D ℓ N | e: <5.280 GeV, μ: <5.174 GeV |
| τ | τ → π N | All flavors: <1.637 GeV (cascade mode) |

### Usage
```bash
./main_hnl_production <mass_GeV> <flavor> [nEvents] [mode]

# Examples:
./main_hnl_production 1.0 muon 100000 direct      # D/B mesons
./main_hnl_production 0.5 tau 100000 fromTau      # τ cascade
```

### Parallel Orchestration
```bash
./run_parallel_production.sh [electron|muon|tau|all]

# Configuration in script:
NEVENTS=100000
MAX_PARALLEL=8
```

### Output Format
File: `output/csv/simulation/HNL_{mass}GeV_{flavor}_{mode}.csv`

```csv
event,weight,hnl_id,parent_pdg,pt,eta,phi,p,E,mass,prod_x_mm,prod_y_mm,prod_z_mm,beta_gamma
1,1.0,9900012,511,45.2,0.123,2.456,150.5,152.3,1.0,0.001,0.002,0.003,145.5
```

| Column | Description |
|--------|-------------|
| `event` | Event number (groups multi-HNL events) |
| `parent_pdg` | Parent meson PDG ID (321=K±, 411=D±, 421=D0, 511=B0, 521=B±, etc.) |
| `pt, eta, phi` | HNL kinematics in lab frame |
| `p, E, mass` | 4-momentum components |
| `prod_x/y/z_mm` | Production vertex in mm (IP = origin) |
| `beta_gamma` | Boost factor β*γ = p/m (for decay length scaling) |

### Key Design Decisions
1. **Pythia as kinematic generator only** - cross-sections and BRs applied in analysis
2. **HNL PDG ID = 9900012** - unified with MadGraph
3. **HNL stable in Pythia** - decay handled analytically in analysis
4. **Per-parent output** - enables correct statistical treatment

---

## Stage 1B: MadGraph Production (Electroweak)

### Purpose
Generate HNL kinematics from W/Z boson decays at 14 TeV LHC using MadGraph 3.6.6.

### Physics Processes
```
W± → e± N     W± → μ± N     W± → τ± N
Z  → νe N̄    Z  → νμ N̄    Z  → ντ N̄
```

### Docker Setup
```bash
cd production/madgraph_production
docker build -t mg5-hnl:latest .
docker run -it --rm -v $(pwd)/../..:/work mg5-hnl:latest bash
```

### Usage
```bash
python3 scripts/run_hnl_scan.py --flavour muon --masses 5 10 15 --nevents 100000
```

### MadGraph Cards
- `cards/proc_card_*.dat` - Process definitions
- `cards/run_card_template.dat` - Run configuration
- `cards/param_card_template.dat` - HNL mass and mixing parameters

### Mixing Convention
Generation uses unit coupling (|U|² = 1), rescaled in analysis:
```python
MIXING_CONFIGS = {
    'electron': {'ve1': 1.0, 'vmu1': 0.0, 'vtau1': 0.0},
    'muon':     {'ve1': 0.0, 'vmu1': 1.0, 'vtau1': 0.0},
    'tau':      {'ve1': 0.0, 'vmu1': 0.0, 'vtau1': 1.0},
}
```

### Output Format
File: `output/csv/simulation/HNL_{mass}GeV_{flavor}_ew.csv`

Same CSV format as Pythia, with `parent_pdg = 24` (W±) or `23` (Z).

---

## Combining Production Channels

### Purpose
Merge Pythia (meson) and MadGraph (EW) files at overlapping masses (4-8 GeV).

### Script
`analysis_pbc/limits/combine_production_channels.py`

### Usage
```bash
# Dry run
python combine_production_channels.py --dry-run

# Combine all
python combine_production_channels.py

# Single flavor
python combine_production_channels.py --flavour muon
```

### Algorithm
1. Find all CSV files for same (mass, flavor)
2. Concatenate DataFrames (preserve all rows)
3. Write combined file: `HNL_{mass}GeV_{flavor}_combined.csv`
4. Delete originals (optional `--no-cleanup` to keep)

### Why This Matters
Without combining, analysis **undercounts by ~50%** at overlapping masses where both meson and EW production contribute.

---

## Stage 2: Geometry Ray-Tracing

### Purpose
Determine geometric acceptance: which HNLs hit the detector and with what path length.

### Script
`analysis_pbc/geometry/per_parent_efficiency.py`

### Detector Geometry
```
Location: CMS drainage gallery
Position: z = 22 m from IP
Shape: Cylindrical tube
Radius: ~1.54 m
Implementation: 40-vertex polyline, triangulated mesh
```

### Algorithm
1. Convert HNL (η, φ) → 3D direction unit vector
2. Ray-trace from production vertex
3. Find intersection with detector mesh
4. Record entry distance and path length

### Output Columns (added to simulation CSV)
| Column | Description |
|--------|-------------|
| `hits_tube` | Boolean: ray intersects detector |
| `entry_distance` | Distance from IP to detector entry (meters) |
| `path_length` | Path length inside detector (meters) |

### Caching
Files cached in `output/csv/geometry/HNL_{mass}GeV_{flavor}_geom.csv`
- Built once on first analysis run
- Reused for all subsequent runs

---

## Stage 3: Limits Calculation

### Main Driver
`analysis_pbc/limits/run.py`

### Usage
```bash
# Single-threaded
python limits/run.py

# Parallel
python limits/run.py --parallel --workers 8

# Dirac interpretation (×2 yield)
python limits/run.py --dirac
```

### Algorithm
```
For each mass × flavor:
    1. Load geometry DataFrame
    2. Scan |U|² from 10⁻¹² to 10⁻² (100 log-spaced points)
    3. For each |U|²:
        N_sig = Σ_parents [ L × σ_parent × BR(parent→ℓN) × P_decay × ε_geom ]
    4. Find eps2_min, eps2_max where N_sig ≥ 2.996
```

### Physics Kernel
`analysis_pbc/limits/expected_signal.py`

```python
def expected_signal_events(geom_df, mass_GeV, eps2, benchmark, lumi_fb, dirac=False):
    """
    N_sig = Σ_parents [ L × σ_parent × BR(parent→ℓN) × P_decay × ε_geom ]

    P_decay = exp(-entry_dist/λ) × [1 - exp(-path_length/λ)]
    λ = β*γ × c*τ₀(m, |U|²)
    """
```

### Key Formula: Decay Probability
```
P_decay = exp(-d_entry/λ) × [1 - exp(-d_path/λ)]

where:
  λ = β*γ × c × τ₀(m, |U|²)    [decay length in lab frame]

  First term: Probability to survive until detector entry
  Second term: Probability to decay while inside detector
```

### Per-Parent Counting (Critical!)
```
CORRECT:
  N_sig = Σ_parents [ L × σ_parent × BR(parent→ℓN) × ε(parent) ]

WRONG (per-event):
  P_event = 1 - Π(1 - P_i)  ← Requires single σ, loses ~75% sensitivity
```

Each parent counted independently with its own cross-section. This is the standard approach in MATHUSLA, ANUBIS, CODEX-b analyses.

### Output
File: `output/csv/analysis/HNL_U2_limits_summary.csv`

```csv
mass_GeV,flavour,benchmark,eps2_min,eps2_max,peak_events
0.20,electron,100,1.12e-07,0.01000,2.17e+10
2.60,muon,010,5.46e-09,9.54e-05,2.81e+05
```

---

## HNL Physics Model

### Wrapper
`analysis_pbc/models/hnl_model_hnlcalc.py`

### Interface
```python
class HNLModel:
    def __init__(self, mass_GeV, Ue2, Umu2, Utau2):
        """Initialize with mass and mixing parameters"""

    @property
    def ctau0_m(self):
        """Proper decay length in meters"""
        # Scales as 1/|U|²

    def production_brs(self):
        """Returns {parent_pdg: BR(parent→ℓN)}"""
```

### HNLCalc
External package from arXiv:2405.07330 implementing 150+ decay modes.

Located in: `analysis_pbc/HNLCalc/HNLCalc.py`

---

## Cross-Sections Database

### File
`analysis_pbc/config/production_xsecs.py`

### Values (14 TeV LHC)
```python
# Heavy quark production
SIGMA_CCBAR_PB = 2.4e10    # pp → cc̄
SIGMA_BBBAR_PB = 5.0e8     # pp → bb̄
SIGMA_KAON_PB  = 5.0e10    # pp → K± (soft QCD)

# Electroweak
SIGMA_W_PB = 2.0e8         # pp → W±
SIGMA_Z_PB = 6.0e7         # pp → Z

# Fragmentation fractions
FRAG_C_D0    = 0.59        # c̄ → D0
FRAG_C_DPLUS = 0.24        # c → D+
FRAG_C_DS    = 0.10        # c → Ds
FRAG_B_B0    = 0.40        # b̄ → B0
FRAG_B_BPLUS = 0.40        # b → B+
```

---

## Stage 4: Money Plot

### Script
`money_plot/plot_money_island.py`

### Input
`output/csv/analysis/HNL_U2_limits_summary.csv`

### Algorithm
```python
for flavor in ['electron', 'muon', 'tau']:
    1. Filter by flavor
    2. Remove NaN rows (no sensitivity)
    3. Deduplicate by mass (keep tightest limits)
    4. Plot exclusion region
```

### Plot Structure
```
        |U|²
          │
    10⁻² ─┼────────────────────────
          │   Upper boundary (too prompt)
          │   ╔═══════════════════╗
          │   ║  EXCLUDED REGION  ║
          │   ╚═══════════════════╝
          │   Lower boundary (too long-lived)
   10⁻¹² ─┼────────────────────────
          │
          └────┬─────┬─────┬─────┬──── m_N [GeV]
              0.1   1     10   100
```

### Physics Interpretation
- **Upper boundary**: Large |U|² → short lifetime → HNL decays before reaching detector
- **Lower boundary**: Small |U|² → long lifetime → HNL passes through detector without decaying
- **Island shape**: Detector provides "window" in lifetime space

---

## Mass Grid Configuration

### File
`config_mass_grid.py`

### Structure
```python
_COMMON_GRID = [
    # K-regime: 0.05 GeV steps
    0.20, 0.25, 0.30, ..., 0.50,

    # D-regime: 0.1-0.2 GeV steps
    0.60, 0.70, ..., 2.00,

    # B-regime: 0.2 GeV steps
    2.20, 2.40, ..., 5.00,

    # EW transition: 0.2 GeV steps
    5.20, 5.40, ..., 8.00,

    # High-mass EW: 0.5 GeV steps
    8.50, 9.00, ..., 17.00,
]
```

### Usage
```python
from config_mass_grid import get_mass_grid

masses = get_mass_grid('muon', 'meson')     # Meson production only
masses = get_mass_grid('muon', 'ew')        # Electroweak only
masses = get_mass_grid('muon', 'combined')  # Both
```

---

## Key Constants

### Luminosity
```python
lumi_fb = 3000  # HL-LHC full run
```

### Exclusion Threshold
```python
N_limit = 2.996  # 3σ on background (Poisson)
```

### Particle Masses (GeV)
```
Electron:  0.000511    Muon:  0.10566    Tau:  1.777
Pion:      0.140       Kaon:  0.494
D±:        1.870       D0:    1.865      Ds:   1.968
B±:        5.279       B0:    5.280      Bs:   5.367
```

---

## Running the Full Pipeline

### 1. Pythia Production (~25 hours)
```bash
cd production/pythia_production
./run_parallel_production.sh all
```

### 2. MadGraph Production (~12 hours)
```bash
cd production/madgraph_production
docker run -it --rm -v $(pwd)/../..:/work mg5-hnl:latest bash
cd /work/production/madgraph_production
python3 scripts/run_hnl_scan.py --flavour electron
python3 scripts/run_hnl_scan.py --flavour muon
python3 scripts/run_hnl_scan.py --flavour tau
```

### 3. Combine Production Files (~5 min)
```bash
cd analysis_pbc/limits
python combine_production_channels.py
```

### 4. Run Analysis (~45 min)
```bash
cd analysis_pbc/limits
python run.py --parallel
```

### 5. Generate Money Plot (<5 sec)
```bash
cd money_plot
python plot_money_island.py
```

### Output
`output/images/HNL_moneyplot_island.png`

---

## Extending to ALPs

To add ALP support:

### 1. Create ALP Model
Create `analysis_pbc/models/alp_model.py`:
```python
class ALPModel:
    def __init__(self, mass_GeV, f_a):
        """
        mass_GeV: ALP mass
        f_a: Decay constant (GeV)
        """

    @property
    def ctau0_m(self):
        """ALP proper decay length"""
        # Γ ∝ m³/f_a² for a → γγ

    def production_brs(self):
        """Production branching ratios"""
        # Depends on ALP-photon, ALP-gluon couplings
```

### 2. Add ALP Production
- **Pythia**: Meson radiative decays (K → π a, B → K a)
- **MadGraph**: Photon fusion (γγ → a), gluon fusion (gg → a)

### 3. Update Analysis
Modify `expected_signal.py` to support ALP model alongside HNL.

### 4. Create ALP Limits Script
Copy `run.py` to `run_alp.py`, swap model class.

### 5. Update Money Plot
Add ALP plotting alongside HNL exclusion regions.

---

## File Format Reference

### Production CSV
```csv
event,weight,hnl_id,parent_pdg,pt,eta,phi,p,E,mass,prod_x_mm,prod_y_mm,prod_z_mm,beta_gamma
```

### Geometry CSV (adds columns)
```csv
...,hits_tube,entry_distance,path_length
```

### Limits Summary CSV
```csv
mass_GeV,flavour,benchmark,eps2_min,eps2_max,peak_events
```

---

## References

- **Pythia 8.315**: http://home.thep.lu.se/~torbjorn/pythia/
- **MadGraph 3.6.6**: https://launchpad.net/mg5amcatnlo
- **HNLCalc**: https://github.com/laroccod/HNLCalc (arXiv:2405.07330)
- **Methodology**:
  - MATHUSLA: arXiv:1811.00927
  - ANUBIS: arXiv:1909.13022
  - CODEX-b: arXiv:1911.00481
  - PBC: CERN-PBC-REPORT-2018-007
