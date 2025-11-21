# Long-Lived Particle Detector Simulation

## Project Overview

This project simulates Heavy Neutral Leptons (HNLs) and other Long-Lived Particles (LLPs) for a proposed detector in the CMS "drainage gallery" at the LHC.

## Detector Concept

We propose a transverse long-lived particle detector located in the CMS drainage gallery, approximately 20m above the CMS interaction point. This location provides:

- **Shielding**: Concrete and earth shielding from both surface and collisional backgrounds
- **Large Active Volume**: Capacity for up to 800 m³ of active detector volume
- **Accessibility**: Capability to construct, power and operate the detector infrastructure

## Physics Goals

The detector aims to search for LLPs using common benchmark models and demonstrate competitive sensitivity compared to existing proposals. The primary focus includes:

- Heavy Neutral Leptons (HNLs) coupling to electrons, muons, and taus
- LLPs produced in proton-proton collisions at √s = 14 TeV
- Mass range scanning from very low masses (0.2 GeV) to high masses (60 GeV)
- Lepton flavor-dependent sensitivity studies

## Simulation Framework

### Main Components

1. **Event Generation**: Uses PYTHIA 8.315 for pp collisions and B-meson/W/Z production
2. **Multi-Lepton Support**: Systematic scans for electron, muon, and tau couplings
3. **Mass Scan**: Lepton-specific mass ranges optimized for coverage
4. **Decay Analysis**: Python-based lifetime and geometric acceptance calculations
5. **Background Estimation**: Simulation of SM backgrounds and detector response

### Project Structure

```
llpatcolliders/
├── CLAUDE.md                                  # Project documentation
├── environment.yml                            # Conda environment specification
├── tmp/                                       # Temporary configuration files (auto-cleaned)
├── output/                                    # All output data organized by type
│   ├── csv/
│   │   ├── simulation/                        # Raw simulation data
│   │   │   ├── HNL_mass_*.csv                # Event data
│   │   │   └── HNL_mass_*.meta               # Cross-section metadata
│   │   └── analysis/                          # Processed analysis results
│   │       ├── HNL_mass_*_BR_vs_ctau.csv
│   │       └── HNL_U2_limits_summary.csv
│   ├── logs/
│   │   ├── simulation/                        # Simulation execution logs
│   │   └── analysis/                          # Analysis execution logs
│   └── images/                                # Plots and visualizations
├── production/                                # Production simulation code
│   ├── main_hnl_single.cc                    # Main C++ simulation
│   ├── main_hnl_single                       # Compiled executable
│   ├── hnl_LowMass_Inclusive_Template.cmnd   # B-meson production config
│   ├── hnl_HighMass_Inclusive_Template.cmnd  # W/Z production config
│   ├── make.sh                               # Universal compilation & scan script
│   ├── scan_electron.sh                      # Electron mass scan
│   ├── scan_muon.sh                          # Muon mass scan
│   └── scan_tau.sh                           # Tau mass scan
└── analysis/                                  # Analysis and plotting scripts
    ├── decayProbPerEvent.py                  # Decay probability calculation
    ├── hnl_coupling_limit.py                 # Coupling limit analysis
    ├── plot_money_island.py                  # Island exclusion plot
    └── plot_money.py                         # Money plot generator
```

### Lepton Flavor Coverage

**Dense Mass Grid (102 total points)**

**Electrons:** 38 mass points (0.2 - 80.0 GeV)
- Dense spacing at low mass where K/D/B interplay is complex
- Mass points: 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.3, 2.6, 3.0, 3.4, 3.8, 4.2, 4.6, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 15.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0 GeV

**Muons:** 38 mass points (0.2 - 80.0 GeV)
- Same grid as electrons for consistency
- Mass points: 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.3, 2.6, 3.0, 3.4, 3.8, 4.2, 4.6, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 15.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0 GeV

**Taus:** 26 mass points (0.5 - 80.0 GeV)
- Fills previously skipped 3.6-5.2 GeV gap
- Mass points: 0.5, 0.7, 1.0, 1.3, 1.6, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 15.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0 GeV

**Rationale:**
- Quasi-logarithmic spacing: dense at low mass (0.05 GeV steps), coarser at high mass
- Captures all meson resonances (K, D, D*, Ds, B)
- Resolves boost/lifetime/geometry interplay in K-dominance (0.2-0.5 GeV), D-dominance (0.5-2 GeV), and B-dominance (2-5 GeV) regions
- Extends to 80 GeV for full EW production coverage (W/Z/top)

## Usage

### Running Production Simulations

Navigate to the `production/` directory:
```bash
cd production
```

Run mass scans:
```bash
# Single lepton flavor
./make.sh electron
./make.sh muon
./make.sh tau

# All three leptons
./make.sh all

# Or use dedicated scripts
./scan_electron.sh
./scan_muon.sh
./scan_tau.sh

# Single mass point
./main_hnl_single 1.0 muon
```

### Running Analysis

Navigate to the `analysis/` directory:
```bash
cd analysis
```

Process simulation data:
```bash
# Decay probability analysis
python decayProbPerEvent.py ../output/csv/simulation/HNL_mass_1.0_electron.csv \
       --meta ../output/csv/simulation/HNL_mass_1.0_electron.meta

# Coupling limit analysis
python hnl_coupling_limit.py

# Generate plots
python plot_money_island.py
python plot_money.py
```

## Simulation Strategy

### Two-Stage Approach

The simulation uses a clean separation between production and decay:

**Stage 1: Production in Pythia (C++)**
1. Generation of b-bbar pairs (low mass) or W/Z bosons (high mass)
2. Hadronization and B-meson production
3. Forced production decays: B → D ℓ N or W/Z → ℓ N
4. **HNL is stable in Pythia** (`mayDecay = off`)
5. Records HNL kinematics (pT, η, φ, production vertex)

**Stage 2: Decay and Acceptance in Python**
1. Python layer (`decayProbPerEvent.py`) handles all decay calculations
2. Uses proper HNL lifetime scaling with momentum/mass
3. Ray-tracing through detector geometry
4. Calculates decay probabilities within detector volume
5. Generates BR vs cτ exclusion limits

**Why this approach?**
- Pythia production cross-sections are accurate and mass-dependent
- Python provides flexible lifetime scanning without re-running Pythia
- Clean separation: production physics (Pythia) vs detector effects (Python)
- Matches methodology of ANUBIS, MATHUSLA, CODEX-b proposals

### Cross-Section Methodology

**IMPORTANT:** The simulation uses **forced decays** (BR=1.0) to decouple event kinematics from coupling strength.

**What the simulation provides:**
- CSV files contain **event kinematics** (pT, η, φ, production vertex)
- CSV files contain **parent_id** to track which meson/boson produced the HNL
- These are generated with forced BR(parent→HNL) = 1.0 (artificial, for sampling efficiency)

**What the analysis must calculate:**
The physical HNL production cross-section requires:

```
σ_HNL(m, |U|²) = σ_parent × BR_real(parent → ℓN, |U|²)

where:
- σ_parent: Parent production cross-section from literature or HNLCalc
  - σ(pp → K): ~10⁷ pb (QCD)
  - σ(pp → D): ~10⁶ pb (charm production)
  - σ(pp → B): ~10⁵ pb (bottom production)
  - σ(pp → W): ~10⁴ pb (EW)
  - σ(pp → Z): ~10³ pb (EW)

- BR_real(parent → ℓN, |U|²): Real branching ratio from HNLCalc
  - BR_real ∝ |U|² (for small mixing)
  - Typical values: 10⁻⁶ to 10⁻¹⁰ for mesons, |U|² for W/Z
```

**Analysis workflow:**
1. Use CSV kinematics to calculate **ε_geom(parent, cτ)** via ray-tracing
2. Use HNLCalc to get **σ_parent** and **BR(parent→ℓN, |U|²)**
3. Calculate signal: **N_sig = L × Σ_parents [σ_parent × BR × ε_geom]**
4. Solve for |U|² limit: **N_sig(|U|²) = N_limit** (typically 3 events at 90% CL)

**Why forced decays?**
- Generates clean kinematic samples without wasting CPU on rare events
- Identical methodology to ANUBIS, MATHUSLA, CODEX-b, PBC proposals
- Separates "production physics" (HNLCalc) from "detector effects" (geometry)

## Technical Implementation

### Parallel Execution

The simulation supports parallel execution with up to 10 concurrent processes (configurable in `make.sh`):
- **Unique Temporary Files**: Each mass point generates unique config files in `tmp/` directory
- **Independent Processes**: No shared file I/O between mass points
- **RAII Cleanup**: Temporary files deleted automatically, even on crashes

### HNL Particle Selection

The HNL is configured as a **stable particle** in Pythia (`9900015:mayDecay = off`):
- Clean separation: production (Pythia) vs decay/acceptance (Python)
- Flexible lifetime scanning without re-running Pythia
- Simple particle selection: just `isFinal()` check
- Python handles all decay calculations

### Output Data

- **102 CSV files**: Complete simulation results across all mass points and lepton flavors
  - 38 electron + 38 muon + 26 tau = 102 total
- **102 log files**: Detailed execution logs (one per CSV)
- Each CSV contains: event weights, particle IDs, kinematics (pT, η, φ), production vertices, parent_id
- Total data size: ~2 GB raw simulation output (~20 MB per mass point)

**CSV Format:**
```
event,weight,id,parent_id,pt,eta,phi,momentum,energy,mass,prod_x_m,prod_y_m,prod_z_m
0,1.0,9900015,321,2.34,1.56,0.78,5.67,6.12,1.0,0.0,0.0,0.0
```
- `parent_id`: PDG code of parent meson/boson (321=K+, 411=D+, 521=B+, 24=W, 23=Z, etc.)
- Kinematics: Used for boost calculations and geometric acceptance
- Production vertex: Used for ray-tracing through detector geometry

## Analysis Approach

The analysis strategy focuses on:
- Design optimization of the detector geometry
- Lepton flavor-dependent sensitivity comparisons
- Background characterization and mitigation
- Signal efficiency calculations across all three lepton flavors
- Sensitivity reach for various LLP benchmarks
- Comparison with existing LLP detector proposals

## PBC-Grade Analysis Pipeline

### HNLCalc Integration

**Location:** `analysis_pbc_test/` directory contains the PBC-style analysis pipeline

**HNLCalc Repository:** https://github.com/laroccod/HNLCalc
- Comprehensive HNL calculator with 150+ production modes and 100+ decay modes
- Covers HNL masses from 0.2 GeV to 10 GeV
- Arbitrary mixing patterns (ve, vmu, vtau)
- Based on arXiv:2405.07330 (Feng, Hewitt, Kling, La Rocco)

**Installation:**
```bash
cd analysis_pbc_test
git clone https://github.com/laroccod/HNLCalc.git
conda run -n llpatcolliders pip install sympy mpmath particle numba
conda run -n llpatcolliders pip install 'scikit-hep==0.4.0'  # Old version with skhep.math
```

**Key HNLCalc Methods:**
- `gen_ctau(mass)` - Proper lifetime at |U|² = 1
- `parent_br(pdg, mass)` - Production branching ratio for parent → ℓN
- `gen_widths(mass)` - Total HNL width
- `Gamma_*()` - Individual decay channel widths

### PBC Analysis Structure

```
analysis_pbc_test/
├── HNLCalc/                          # HNLCalc repository (cloned)
│   └── HNLCalc.py                   # Main module (~3000 lines)
├── models/
│   └── hnl_model_hnlcalc.py         # Wrapper around HNLCalc
├── geometry/
│   └── per_parent_efficiency.py     # ε_geom(parent, cτ) calculator
├── limits/
│   └── u2_limit_calculator.py       # Direct |U|² limit solver
└── test_hnlcalc.py                  # Installation verification script
```

### PBC Methodology

**Per-Parent Efficiency Maps:**
- Separate ε_geom(K), ε_geom(D), ε_geom(B), ε_geom(W), ε_geom(Z)
- Accounts for different boost distributions per parent species
- Uses ray-tracing through detector geometry
- Scans 100 log-spaced cτ points from 1 mm to 1 km

**Direct |U|² Limit Calculation:**
```
N_sig(m, flavour, |U|²) = L × Σ_parents [σ_p × BR(p→ℓN, |U|²) × ε_geom(p, cτ(|U|²))]

Solve: N_sig(m, flavour, |U|²) = N_limit  (typically 3 events for 90% CL)
```

**Key Features:**
- No BR_limit intermediate step (direct |U|² calculation)
- Uses real HNLCalc physics (widths, BRs, production rates)
- Per-parent tracking via `parent_id` column in CSV
- Proper lifetime scaling: cτ(|U|²) = cτ₀/|U|²

### Benchmark Curves

**BC6:** Electron coupling (ve=1, vmu=0, vtau=0)
**BC7:** Muon coupling (ve=0, vmu=1, vtau=0)
**BC8:** Tau coupling (ve=0, vmu=0, vtau=1)

Comparing drainage gallery sensitivity to PBC proposals (ANUBIS, MATHUSLA, CODEX-b, AL3X, etc.)

## Current Status

**Production:** Dense mass grid simulation in progress (102 mass points)
- 38 electron + 38 muon + 26 tau = 102 total
- Using unified K+D+B meson production (< 5 GeV) and W/Z EW production (≥ 5 GeV)
- Generating ~200k events per mass point

**Analysis:** PBC pipeline framework complete
- HNLCalc successfully installed and tested
- Per-parent efficiency calculator ready
- Direct |U|² limit solver implemented
- Awaiting simulation completion for full analysis run
