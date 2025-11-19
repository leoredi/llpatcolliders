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
- Mass range scanning from very low masses (0.2 GeV) to high masses (80 GeV)
- Lepton flavor-dependent sensitivity studies

## Simulation Framework

### Main Components

1. **Event Generation**: Uses PYTHIA 8.315 for pp collisions and B-meson production
2. **Multi-Lepton Support**: Systematic scans for electron, muon, and tau couplings
3. **Mass Scan**: Lepton-specific mass ranges optimized for kinematic constraints
4. **Decay Analysis**: Tracking of HNL production and decay within detector acceptance
5. **Background Estimation**: Simulation of SM backgrounds and detector response

### Key Files

#### Simulation Code
- `main_hnl_single.cc`: Main C++ simulation with lepton flavor parameter support

#### Configuration Templates
- `hnl_LowMass_Inclusive_Template.cmnd`: Template for low-mass HNLs (B-meson production)
- `hnl_HighMass_Inclusive_Template.cmnd`: Template for high-mass HNLs (W/Z production)

#### Mass Scan Scripts
- `make.sh`: Universal compilation and scan script (supports all leptons)
- `scan_electron.sh`: Electron-specific mass scan
- `scan_muon.sh`: Muon-specific mass scan
- `scan_tau.sh`: Tau-specific mass scan

### Lepton Flavor Capabilities

The simulation supports three lepton flavors with optimized mass ranges:

**Electrons:**
- Mass points: 0.2, 0.3, 0.4, 0.5, 1.0, 1.5, 2.0, 2.5, 2.8, 3.1, 3.3, 3.6, 4.0, 4.4, 4.8, 5.2, 10.0, 15.0, 20.0, 40.0, 80.0 GeV
- Covers masses below muon threshold
- 21 mass points

**Muons:**
- Mass points: 0.5, 1.0, 1.5, 2.0, 2.5, 2.8, 3.1, 3.3, 3.6, 4.0, 4.4, 4.8, 5.2, 10.0, 15.0, 20.0, 40.0, 80.0 GeV
- Standard range starting from muon mass
- 18 mass points

**Taus:**
- Low-mass: 0.5, 1.0, 1.5, 2.0, 2.5, 2.8, 3.1, 3.4 GeV (B-meson decays)
- High-mass: 10.0, 15.0, 20.0, 40.0, 80.0 GeV (W/Z decays)
- Skips 3.6-5.2 GeV (zero yield due to tau mass threshold)
- 13 mass points

### Usage

Run simulations using the universal script:
```bash
./make.sh electron    # Electron scan
./make.sh muon        # Muon scan (default)
./make.sh tau         # Tau scan
./make.sh all         # All three leptons
```

Or use dedicated scripts:
```bash
./scan_electron.sh
./scan_muon.sh
./scan_tau.sh
```

Single mass point:
```bash
./main_hnl_single <mass_GeV> [electron|muon|tau]
```

### Simulation Strategy

The simulation performs:
1. Generation of b-bbar pairs (low mass) or W/Z bosons (high mass)
2. Hadronization and B-meson production
3. Lepton-flavor-specific HNL production channels with proper kinematic thresholds
4. Tracking of HNL decay products (N → ℓ π)
5. Geometric acceptance cuts for the detector location
6. Statistical analysis for sensitivity projections

#### Kinematic Threshold Logic

The simulation implements lepton-flavor-dependent thresholds for transitioning from 3-body (B → D ℓ N) to 2-body (B → ℓ N) decays:

- **Electrons/Muons**: Threshold at 3.3 GeV (m_B - m_D - m_μ)
- **Taus**: Threshold at 1.65 GeV (m_B - m_D - m_τ ~ 5.28 - 1.87 - 1.77 GeV)

This ensures kinematically allowed decays across all mass points and prevents zero-event yields.

## Analysis Approach

The analysis strategy focuses on:
- Design optimization of the detector geometry
- Lepton flavor-dependent sensitivity comparisons
- Background characterization and mitigation
- Signal efficiency calculations across all three lepton flavors
- Sensitivity reach for various LLP benchmarks
- Comparison with existing LLP detector proposals

## Current Status

The simulation framework is fully operational with multi-lepton support. Mass scans can be performed for electron, muon, and tau couplings to determine flavor-dependent detector sensitivity and physics reach.

## Project Structure

```
llpatcolliders/
├── CLAUDE.md                                  # Project documentation
├── main_hnl_single.cc                         # Main C++ simulation (with RAII cleanup)
├── main_hnl_single                            # Compiled executable
├── make.sh                                    # Universal compilation & scan script
├── scan_electron.sh                           # Electron mass scan script
├── scan_muon.sh                               # Muon mass scan script
├── scan_tau.sh                                # Tau mass scan script
├── run_decay_analysis.sh                      # Run decay probability analysis
├── hnl_LowMass_Inclusive_Template.cmnd        # B-meson production config template
├── hnl_HighMass_Inclusive_Template.cmnd       # W/Z production config template
├── hnl_coupling_limit.py                      # Coupling limit analysis script
├── decayProbPerEvent.py                       # Decay probability calculation
├── environment.yml                            # Conda environment specification
├── csv/                                       # Simulation output data
│   ├── HNL_mass_0.2_electron.csv
│   ├── HNL_mass_0.5_muon.csv
│   ├── HNL_mass_1.0_tau.csv
│   └── ... (52 CSV files total)
├── logs/                                      # Execution logs per mass point
│   ├── log_electron_0.2.txt
│   ├── log_muon_0.5.txt
│   ├── log_tau_1.0.txt
│   └── ... (52 log files total, one per CSV)
└── output/                                    # Analysis results
    ├── csv/                                   # Processed analysis data
    └── images/                                # Plots and visualizations
```

### Output Data Summary

- **52 CSV files**: Complete simulation results for all mass points across three lepton flavors
  - 21 electron mass points (0.2-80.0 GeV)
  - 18 muon mass points (0.5-80.0 GeV)
  - 13 tau mass points (0.5-80.0 GeV, optimized for kinematic constraints)
- **52 log files**: Detailed execution logs matching each CSV output (one-to-one correspondence)
- Each CSV contains: event weights, particle IDs, kinematics (pT, η, φ), decay vertices, and flight distances

## Technical Implementation Details

### Parallel Execution Safety

The simulation is designed for safe parallel execution with up to 10 concurrent processes (configurable in `make.sh`). Key safety features include:

- **Unique Temporary Files**: Each mass point generates unique temporary configuration files (e.g., `hnl_LowMass_muon_1.0_temp.cmnd`) to prevent race conditions
- **Independent Processes**: Each mass point runs as a completely independent process with no shared file I/O
- **Separate Output Files**: Results are written to unique CSV files per mass and lepton flavor

### Recent Updates (2025-11-19)

#### Bug Fixes

Two critical bugs were identified and fixed:

1. **Race Condition (FIXED)**: Temporary configuration files now include the mass value in their filenames, preventing file corruption when multiple processes run in parallel. This was essential for `./make.sh all` to function correctly.

2. **Tau Kinematic Threshold (FIXED)**: The threshold for switching from 3-body to 2-body B-meson decays is now lepton-flavor-dependent:
   - Taus: 1.65 GeV (accounts for heavy tau mass)
   - Electrons/Muons: 3.3 GeV (original threshold)

   This fix prevents zero-event yields for tau mass points in the 2.0-3.4 GeV range.

#### RAII-Based Temporary File Cleanup (IMPLEMENTED)

The simulation now implements automatic temporary file cleanup using RAII (Resource Acquisition Is Initialization):

- **ScopedFileRemover Class**: Added to `main_hnl_single.cc` (lines 84-103)
  - Automatically deletes temporary `.cmnd` configuration files when the program exits
  - Works even in error cases (failed init, exceptions, early returns)
  - Prevents copying to avoid ownership confusion

- **Implementation Details** (`main_hnl_single.cc`):
  ```cpp
  // Instantiated after temporary filenames are generated
  ScopedFileRemover cleanupLow(cardLowMass);
  ScopedFileRemover cleanupHigh(cardHighMass);
  ```

- **Benefits**:
  - Crash-safe: Files deleted even if program exits early
  - No manual cleanup needed before each `return` statement
  - Exception-safe: Destructors run automatically during stack unwinding
  - Parallel-safe: Works with existing unique filename strategy

- **Status**: IMPLEMENTED AND TESTED
  - Existing temporary files from previous runs have been cleaned up
  - New simulations automatically clean up temporary files on exit
  - Zero temporary files remain after execution (verified)
