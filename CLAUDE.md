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

Complete Analysis Pipeline: File-by-File Flow

  Base directory: /Users/fredi/cernbox/Physics/llpatcolliders/llpatcolliders

  ---
  Stage 1: Pythia Event Generation (C++)

  File: production/main_hnl_single.cc
  - Input: None (generates events)
  - Configuration:
    - production/hnl_Meson_Inclusive_Template.cmnd (for m < 5 GeV)
    - production/hnl_HighMass_Inclusive_Template.cmnd (for m ≥ 5 GeV)
  - Output: output/csv/simulation/HNL_mass_{mass}_{flavour}_{regime}.csv
    - Example: output/csv/simulation/HNL_mass_2.6_muon_Meson.csv
    - Columns: event,weight,id,parent_id,pt,eta,phi,momentum,energy,mass,pro
  d_x_m,prod_y_m,prod_z_m

  Compilation script: production/make.sh

  Scan scripts:
  - production/scan_electron.sh
  - production/scan_muon.sh
  - production/scan_tau.sh

  ---
  Stage 2: Geometry Preprocessing (Python)

  File: analysis_pbc_test/geometry/per_parent_efficiency.py

  Function: preprocess_hnl_csv(csv_file, mesh, origin)

  Input: output/csv/simulation/HNL_mass_{mass}_{flavour}_{regime}.csv

  Output: output/csv/geometry/HNL_mass_{mass}_{flavour}_geom.csv
  - Example: output/csv/geometry/HNL_mass_2.6_muon_geom.csv
  - Added columns:
    - hits_tube (bool)
    - entry_distance (m)
    - path_length (m)
    - beta_gamma (dimensionless)

  Helper function: build_drainage_gallery_mesh() - Creates detector mesh at
  z=22m

  ---
  Stage 3: HNL Physics Model (Python)

  File: analysis_pbc_test/models/hnl_model_hnlcalc.py

  Class: HNLModel(mass_GeV, Ue2, Umu2, Utau2)

  Dependencies:
  - analysis_pbc_test/HNLCalc/HNLCalc.py (external repository)

  Methods:
  - .ctau0_m - Returns proper lifetime in meters
  - .production_brs() - Returns dict {parent_pdg: BR(parent→ℓN)}

  No direct input/output - Called by limit calculator

  ---
  Stage 4: Cross-Section Database (Python)

  File: analysis_pbc_test/config/production_xsecs.py

  Function: get_parent_sigma_pb(parent_pdg)

  Input: PDG code (int)

  Output: Production cross-section in pb (float)

  Examples:
  - σ(pp→D⁰) = 2.8 × 10¹⁰ pb
  - σ(pp→B⁰) = 4.0 × 10⁸ pb
  - σ(pp→K⁺) = 5.0 × 10¹⁰ pb

  No file I/O - Hardcoded values from PBC literature

  ---
  Stage 5: Limit Calculation (Python)

  File: analysis_pbc_test/limits/u2_limit_calculator.py

  Main function: run_reach_scan(flavour, benchmark, lumi_fb, n_jobs=4)

  Input files (per mass point):
  - output/csv/simulation/HNL_mass_{mass}_{flavour}_{regime}.csv
  - output/csv/geometry/HNL_mass_{mass}_{flavour}_geom.csv (created if not
  cached)

  Uses:
  - models/hnl_model_hnlcalc.py (HNLModel class)
  - config/production_xsecs.py (cross-sections)
  - geometry/per_parent_efficiency.py (geometry preprocessing)

  Output: output/csv/analysis/HNL_U2_limits_summary.csv
  - Columns:
    - mass_GeV - HNL mass
    - flavour - electron/muon/tau
    - benchmark - 100/010/001 (coupling pattern)
    - eps2_min - Lower |U|² exclusion boundary
    - eps2_max - Upper |U|² exclusion boundary
    - peak_events - Maximum signal events

  Execution:
  cd analysis_pbc_test
  python limits/u2_limit_calculator.py

  ---
  Stage 6: Plotting (Python)

  File: analysis/plot_money_island.py

  Input: output/csv/analysis/HNL_U2_limits_summary.csv

  Output: output/images/HNL_moneyplot_island.png
  - 3-panel plot (electron, muon, tau)
  - Shows exclusion islands (shaded regions)
  - Upper/lower boundaries for each flavor

  Execution:
  cd analysis
  python plot_money_island.py

  ---
  Complete File Dependency Tree

  production/main_hnl_single.cc
  ├── [reads] production/hnl_Meson_Inclusive_Template.cmnd
  ├── [reads] production/hnl_HighMass_Inclusive_Template.cmnd
  └── [writes] output/csv/simulation/HNL_mass_*.csv
              ↓
  analysis_pbc_test/geometry/per_parent_efficiency.py
  ├── [reads] output/csv/simulation/HNL_mass_*.csv
  └── [writes] output/csv/geometry/HNL_mass_*_geom.csv
              ↓
  analysis_pbc_test/limits/u2_limit_calculator.py
  ├── [reads] output/csv/simulation/HNL_mass_*.csv
  ├── [reads] output/csv/geometry/HNL_mass_*_geom.csv
  ├── [uses] analysis_pbc_test/models/hnl_model_hnlcalc.py
  │          └── [uses] analysis_pbc_test/HNLCalc/HNLCalc.py
  ├── [uses] analysis_pbc_test/config/production_xsecs.py
  └── [writes] output/csv/analysis/HNL_U2_limits_summary.csv
              ↓
  analysis/plot_money_island.py
  ├── [reads] output/csv/analysis/HNL_U2_limits_summary.csv
  └── [writes] output/images/HNL_moneyplot_island.png

  ---
  Key Intermediate Files (Example: 2.6 GeV Muon)

  1. Simulation: output/csv/simulation/HNL_mass_2.6_muon_Meson.csv
    - 8310 rows (HNLs)
    - 13 columns
  2. Geometry: output/csv/geometry/HNL_mass_2.6_muon_geom.csv
    - 8310 rows (same)
    - 17 columns (+4 geometry columns)
  3. Summary: output/csv/analysis/HNL_U2_limits_summary.csv
    - 37 rows (all mass points × flavors)
    - 6 columns (mass, flavour, benchmark, limits, peak)
  4. Final plot: output/images/HNL_moneyplot_island.png
    - 2700×750 px PNG image

  ---
  Supporting Documentation Files

  - CLAUDE.md - Main project documentation
  - analysis_pbc_test/README.md - PBC pipeline guide
  - analysis_pbc_test/VALIDATION.md - Validation report
  - analysis_pbc_test/limits/MULTI_HNL_METHODOLOGY.md - Per-parent counting
  explanation
  - analysis_pbc_test/limits/ROBUSTNESS_FIXES.md - Defensive programming
  guide

  ---
  Test Files

  - analysis_pbc_test/tests/test_pipeline.py - Pipeline smoke tests
  - analysis_pbc_test/tests/test_26gev_muon.py - Benchmark validation

