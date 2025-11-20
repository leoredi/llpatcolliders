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

### Key Files

#### Simulation Code
- `main_hnl_single.cc`: Main C++ simulation with lepton flavor parameter support

#### Configuration Templates
- `hnl_LowMass_Inclusive_Template.cmnd`: Template for low-mass HNLs (B-meson production)
- `hnl_HighMass_Inclusive_Template.cmnd`: Template for high-mass HNLs (W/Z production)

#### Analysis Scripts
- `decayProbPerEvent.py`: Decay probability calculation with detector geometry
- `hnl_coupling_limit.py`: Coupling limit analysis script

#### Mass Scan Scripts
- `make.sh`: Universal compilation and scan script (supports all leptons)
- `scan_electron.sh`: Electron-specific mass scan
- `scan_muon.sh`: Muon-specific mass scan
- `scan_tau.sh`: Tau-specific mass scan

### Lepton Flavor Capabilities

The simulation supports three lepton flavors with optimized mass ranges:

**Electrons:**
- Mass points: 0.2, 0.3, 0.4, 0.5, 1.0, 1.5, 2.0, 2.5, 2.8, 3.1, 3.3, 3.6, 4.0, 4.4, 4.8, 5.2, 10.0, 15.0, 20.0, 40.0, 60.0 GeV
- Covers masses below muon threshold
- 21 mass points

**Muons:**
- Mass points: 0.5, 1.0, 1.5, 2.0, 2.5, 2.8, 3.1, 3.3, 3.6, 4.0, 4.4, 4.8, 5.2, 10.0, 15.0, 20.0, 40.0, 60.0 GeV
- Standard range starting from muon mass
- 18 mass points

**Taus:**
- Low-mass: 0.5, 1.0, 1.5, 2.0, 2.5, 2.8, 3.1, 3.4 GeV (B-meson decays)
- High-mass: 10.0, 15.0, 20.0, 40.0, 60.0 GeV (W/Z decays)
- Skips 3.6-5.2 GeV (kinematic gap region)
- 13 mass points

### Usage

Run simulations using the universal script:
```fish
# Electron scan
./make.sh electron

# Muon scan (default)
./make.sh muon

# Tau scan
./make.sh tau

# All three leptons
./make.sh all
```

Or use dedicated scripts:
```fish
./scan_electron.sh
./scan_muon.sh
./scan_tau.sh
```

Single mass point:
```fish
# Example: 1.0 GeV muon
./main_hnl_single 1.0 muon
```

### Simulation Strategy

#### Two-Stage Approach

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

The simulation uses **Pythia-generated cross-sections** specific to each mass point, lepton flavor, and production mechanism. This approach ensures accurate predictions for HL-LHC sensitivity studies.

#### Implementation

**C++ Simulation (`main_hnl_single.cc`)**:
- After event generation, extracts `pythia.info.sigmaGen()` and `sigmaErr()`
- Converts from millibarns (mb) to picobarns (pb): σ_pb = σ_mb × 10⁹
- Writes cross-section to `.meta` file alongside each CSV output

**Python Analysis (`decayProbPerEvent.py`)**:
- Reads cross-section from `.meta` file for each mass point
- Converts pb → fb for luminosity calculations
- Uses mass-specific σ for exclusion limit calculations

#### Meta File Format

Each `.meta` file contains:
```
# Cross-section information from Pythia 8
# Generated at sqrt(s) = 14 TeV
sigma_gen_pb 2.714698e+07
sigma_err_pb 4.447674e+04
```

#### Production Mechanisms

Cross-sections vary by mass and production channel:

- **Low Mass (< 5 GeV)**: B-meson production (b-bbar → B-mesons → HNL)
  - Example: m_HNL = 0.2 GeV (e-coupling) → σ ≈ 27 nb
  - Dominated by QCD b-quark production
  - Uses 3-body B → D ℓ N decays

- **High Mass (≥ 5 GeV)**: W/Z boson production (pp → W/Z → HNL)
  - Example: m_HNL = 60 GeV → σ varies by lepton flavor
  - Electroweak production mechanism
  - Uses 2-body W → ℓ N and Z → ν N decays

#### Why Pythia Cross-Sections?

Previously, the analysis used a generic hardcoded cross-section (σ = 200 nb for W production). The new approach:

1. **Mass-dependent**: Cross-section changes significantly with HNL mass
2. **Process-specific**: B-meson vs W/Z production have different σ values
3. **Lepton-dependent**: Different kinematic cuts for e/μ/τ affect acceptance
4. **Generator-accurate**: Uses the same σ that generated the events

This matches the methodology used by ANUBIS and other LLP detector proposals, where events are generated at |U|² = 1 (unit mixing), and the true branching ratio BR(parent → HNL) is applied as a reweighting factor in the analysis.

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
├── main_hnl_single.cc                         # Main C++ simulation
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
├── tmp/                                       # Temporary configuration files (auto-cleaned)
└── output/                                    # All output data organized by type
    ├── csv/
    │   ├── simulation/                        # Raw simulation data
    │   │   ├── HNL_mass_0.2_electron.csv     # Event data
    │   │   ├── HNL_mass_0.2_electron.meta    # Cross-section metadata
    │   │   ├── HNL_mass_0.5_muon.csv
    │   │   ├── HNL_mass_0.5_muon.meta
    │   │   ├── HNL_mass_1.0_tau.csv
    │   │   ├── HNL_mass_1.0_tau.meta
    │   │   └── ... (52 CSV + 52 meta files = 104 total)
    │   └── analysis/                          # Processed analysis results
    │       ├── HNL_mass_*_particle_decay_results.csv
    │       ├── HNL_mass_*_event_decay_statistics.csv
    │       └── HNL_mass_*_exclusion_data.csv
    ├── logs/
    │   ├── simulation/                        # Simulation execution logs
    │   │   ├── log_electron_0.2.txt
    │   │   ├── log_muon_0.5.txt
    │   │   ├── log_tau_1.0.txt
    │   │   └── ... (52 log files total, one per CSV)
    │   └── analysis/                          # Analysis execution logs
    │       └── decay_analysis_*GeV_*.log
    └── images/                                # Plots and visualizations
```

### Output Data Summary

- **52 CSV files**: Complete simulation results for all mass points across three lepton flavors
  - 21 electron mass points (0.2-60.0 GeV)
  - 18 muon mass points (0.5-60.0 GeV)
  - 13 tau mass points (0.5-60.0 GeV)
- **52 meta files**: Cross-section metadata from Pythia 8 (one-to-one correspondence with CSV files)
  - Contains σ_gen (generated cross-section in pb) at √s = 14 TeV
  - Mass-specific, lepton-specific, and process-specific (B-meson vs W/Z production)
- **52 log files**: Detailed execution logs matching each CSV output (one-to-one correspondence)
- Each CSV contains: event weights, particle IDs, kinematics (pT, η, φ), production vertices
- Each meta file contains: Pythia-generated cross-section and uncertainty in picobarns

## Technical Implementation Details

### Parallel Execution Safety

The simulation is designed for safe parallel execution with up to 10 concurrent processes (configurable in `make.sh`). Key safety features include:

- **Unique Temporary Files**: Each mass point generates unique temporary configuration files in `tmp/` directory (e.g., `tmp/hnl_LowMass_muon_1.0_temp.cmnd`) to prevent race conditions
- **Independent Processes**: Each mass point runs as a completely independent process with no shared file I/O
- **Separate Output Files**: Results are written to unique CSV files per mass and lepton flavor
- **Automatic Cleanup**: RAII-based cleanup ensures temporary files are deleted even if processes crash or are interrupted

### HNL Particle Selection

The HNL is configured as a **stable particle** in Pythia (`9900015:mayDecay = off`). This design choice ensures:

1. **Clean separation**: Production (Pythia) vs decay/acceptance (Python)
2. **Flexible lifetime scanning**: No need to re-run Pythia for different lifetimes
3. **Simple particle selection**: Just `isFinal()` check selects the physical HNL
4. **No decay branching complications**: Python handles all decay calculations

The C++ code selects HNLs using:
```cpp
if (std::abs(prt.id()) != 9900015) continue;  // Select HNLs
if (!prt.isFinal()) continue;                  // Select final-state copy
```

### Recent Updates (2025-11-20)

#### Stable HNL Implementation (LATEST)

**Major architectural change**: HNL is now fully stable in Pythia, with all decay calculations moved to Python layer.

**Changes to `.cmnd` templates:**
- `9900015:mayDecay = off` - HNL does not decay in Pythia
- `9900015:onMode = off` - No decay channels defined
- Removed all `9900015:addChannel` lines (previous N → ℓ π decays)
- **Production decays unchanged**: B → D ℓ N, W → ℓ N, Z → ν N still forced

**Changes to `main_hnl_single.cc`:**
- Removed lepton-flavor-dependent kinematic threshold logic (previously lines 164-191)
- Removed 3-body → 2-body B decay switching code
- Simplified HNL selection to just `isFinal()` check
- Added `#include <sstream>` for proper compilation

**Benefits:**
- Clean two-stage design: production (Pythia) + decay/geometry (Python)
- No artificial mass thresholds or decay mode switching
- Consistent with ANUBIS/MATHUSLA/CODEX-b methodology
- Python layer has full control over lifetime and acceptance calculations
- Easy to scan multiple lifetimes without re-generating events

#### Pythia Cross-Section Implementation

The simulation uses **Pythia-generated cross-sections** instead of hardcoded values:

1. **C++ Implementation**: `main_hnl_single.cc` (lines 243-264)
   - Extracts `pythia.info.sigmaGen()` and `sigmaErr()` after event generation
   - Converts from millibarns to picobarns
   - Writes cross-section to `.meta` file alongside each CSV

2. **Python Integration**: `decayProbPerEvent.py`
   - Added `read_cross_section_from_meta()` function
   - Automatically reads mass-specific σ from `.meta` files
   - Replaces hardcoded σ = 200 nb with process-specific values

3. **Benefits**:
   - Mass-dependent cross-sections (σ varies from ~27 nb at 0.2 GeV to ~pb levels at 60 GeV)
   - Process-specific (B-meson vs W/Z production)
   - Lepton-flavor dependent (accounts for different kinematic cuts)
   - Matches ANUBIS methodology (generate at |U|² = 1, reweight by BR)

#### Mass Range Update (80 → 60 GeV)

Maximum HNL mass changed from 80 GeV to 60 GeV across all lepton flavors:
- **Reason**: 80 GeV is too close to W boson mass (80.4 GeV), causing generation issues
- **Updated**: All scan scripts (`make.sh`, `scan_*.sh`) now use 60 GeV as maximum
- **Impact**: Total mass points remain at 52, but highest mass is now 60 GeV

#### RAII-Based Temporary File Cleanup

The simulation implements automatic temporary file cleanup using RAII (Resource Acquisition Is Initialization):

- **ScopedFileRemover Class**: Added to `main_hnl_single.cc` (lines 86-104)
  - Automatically deletes temporary `.cmnd` configuration files when the program exits
  - Works even in error cases (failed init, exceptions, early returns)
  - Prevents copying to avoid ownership confusion

- **Temporary File Location**: All temporary configuration files stored in `./tmp/` directory
  - Keeps the working directory clean
  - Directory automatically created if it doesn't exist (`mkdir -p tmp`)
  - Files named: `tmp/hnl_{LowMass|HighMass}_{lepton}_{mass}_temp.cmnd`

- **Benefits**:
  - Crash-safe: Files deleted even if program exits early
  - No manual cleanup needed before each `return` statement
  - Exception-safe: Destructors run automatically during stack unwinding
  - Parallel-safe: Works with existing unique filename strategy
  - Clean workspace: All temporary files isolated in `tmp/`

#### Race Condition Fix

Temporary configuration files now include the mass value in their filenames, preventing file corruption when multiple processes run in parallel. This was essential for `./make.sh all` to function correctly.
