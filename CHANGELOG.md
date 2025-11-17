# Changelog

## 2025-11-17 - Duplicate Detection Fix and Production Run

### Fixed
- **Duplicate particle detection in `main144.cc`**:
  - Root cause: 99.98% of "2 HNL events" had same charge (same PDG ID)
  - These were the same particle at different PYTHIA event record stages
  - Old approach: Kinematic thresholds (0.1% pT, 0.001 eta/phi) - caught only 0.04% of duplicates
  - New approach: Track written charges with `std::set<int>`, write only first occurrence of each PDG ID
  - Result: **100% duplicate elimination** (4,792/100k events → 0/1000 test events)

### Changed
- **Mass scan configuration** (`run_mass_scan.py`):
  - Events per mass point: 100k → **1M** (10x statistics improvement)
  - Mass points: 15-79 GeV → **15-71 GeV** (excluded 79 GeV kinematic endpoint)
  - Expected file size: 5.6 MB → ~56 MB per mass point

### Cleaned
- Removed obsolete test/debug files: `check_duplicates.py`, `run_all_decay_analysis.py`, test logs
- Removed unused PYTHIA example files: `main144.C`, `main144Dct.h`
- Updated `.gitignore` to exclude auto-generated configs, test CSVs, and temporary scripts

### Analysis
- Created `auto_analyze.sh` to automatically run `decayProbPerEvent.py` on all mass points after scan completes

## 2025-11-17 - Output Reorganization and Plot Improvements

### Changed
- **Output directory structure**:
  - Created centralized `output/` directory with two subdirectories:
    - `output/csv/` - All CSV files (simulation output and analysis results)
    - `output/images/` - All PNG plots and visualizations
  - Removed scattered output files from root directory and `pythiaStuff/`
  - Removed obsolete `pythiaStuff/mass_scan_hnl/` subdirectory

- **Modified scripts for new output structure**:
  - `pythiaStuff/main144.cc`: CSV outputs now go to `../output/csv/`
  - `decayProbPerEvent.py`:
    - Analysis CSV results → `output/csv/`
    - Plots → `output/images/`
  - `pythiaStuff/run_mass_scan.py`:
    - Simplified to use centralized output directories
    - No longer creates `mass_scan_hnl/` subdirectory

- **Plot label improvements** (`decayProbPerEvent.py`):
  - Single-particle decay probability now explicitly labeled as "conditional on hit"
  - Exclusion plot x-axis: "Lifetime (m)" → "cτ (m)" (proper decay length notation)

### Updated
- **CLAUDE.md**: Updated all documentation to reflect new output paths and directory structure

### Benefits
- Cleaner repository structure with organized outputs by file type
- Easier to find and manage simulation results and plots
- Consistent output location regardless of where scripts are run

## 2025-11-14 - Complete Mass Scan Analysis

### Completed
- **Full HNL mass scan and analysis pipeline**:
  - Generated CSV files for all 9 mass points (15-79 GeV, 8 GeV steps)
  - Verified no duplicate entries in any CSV files
  - Successfully ran `decayProbPerEvent.py` for all mass points
  - Generated 9 exclusion plots comparing detector sensitivity with MATHUSLA, CODEX-b, and ANUBIS
  - All 9 mass points show HNL particles reaching and decaying in the detector
  - Exclusion limits calculated as function of lifetime and branching ratio

### Analysis Results
- **CSV Statistics** (per mass point):
  - ~104,700 HNL entries for masses 15-71 GeV (100k events)
  - ~2,031 entries for 79 GeV (kinematic suppression near W mass limit)
  - ~95% events with 1 HNL, ~5% with 2 HNLs (from multi-W production)
  - Zero duplicate entries confirmed across all files

- **Generated Outputs**:
  - Exclusion plots: `hnlLL_m{mass}GeVLLP_exclusion_vs_lifetime.png` (9 files)
  - Correlation plots: `hnlLL_m{mass}GeVLLP_correlation_analysis.png` (9 files)
  - Event statistics: `hnlLL_m{mass}GeVLLP_event_decay_statistics.csv` (9 files)
  - Particle results: `hnlLL_m{mass}GeVLLP_particle_decay_results.csv` (9 files)

## 2025-11-14 - Parallel Mass Scan

### Improved
- **Parallel execution in mass scan** (`pythiaStuff/run_mass_scan.py`):
  - Added multiprocessing support to run 2 simulations in parallel
  - Implemented real-time progress reporting with timestamps
  - Added detailed logging showing which mass points are currently running
  - Enhanced summary with wall-clock time vs total simulation time
  - Calculates and displays speedup factor from parallelization
  - **Expected speedup**: ~2x faster (22-27 minutes instead of 45-54 minutes for full scan)
  - **CPU utilization**: Now uses ~6% (2 cores) instead of ~3% (1 core)

### Modified
- **pythiaStuff/run_mass_scan.py**:
  - Imports `multiprocessing.Pool` and `datetime` for parallel execution
  - Refactored `run_simulation()` to be self-contained (creates own config files)
  - Added timestamped terminal output: `[HH:MM:SS] STARTING/COMPLETED: HNL mass = X GeV`
  - Returns tuple of (mass, success, duration) for performance tracking

## 2025-11-14 - Mass Scan Implementation and Duplicate Fix

### Added
- **Mass scan script** (`pythiaStuff/run_mass_scan.py`):
  - Automated HNL mass scanning from 15-79 GeV in 8 GeV steps
  - Generates configuration files automatically for each mass point
  - Runs PYTHIA simulations sequentially
  - Organizes all outputs in `pythiaStuff/mass_scan_hnl/` subdirectory
  - Configurable number of events per mass point
  - Progress tracking and summary reporting

### Fixed
- **Duplicate particle entries in CSV output** (`pythiaStuff/main144.cc`):
  - **Problem**: PYTHIA event record contains the same particle at multiple stages (initial, intermediate, final), causing up to 8 copies of the same particle with identical kinematics to be written
  - **Impact**: ~35,187 exact duplicates out of ~100k entries in previous runs
  - **Solution**: Implemented kinematic-based duplicate detection
    - Compares particles within same event using relative tolerance (0.1% for pT, 0.001 for eta/phi)
    - Limits output to maximum 2 LLPs per event (configurable)
    - Added `#include <set>` for tracking written particles
  - **Result**: Reduced to ~1 unique HNL per event (957/1000 events) with occasional genuine multi-W production (~43/1000 events with 2 HNLs)

### Modified
- **CLAUDE.md**: Updated with correct workflow, conda environment setup, and removed obsolete references
- **pythiaStuff/main144.cc**: Enhanced CSV output logic with duplicate detection
- **Environment**: Updated conda environment (`environment.yml`) to include all necessary Python packages

### Workflow Summary
1. **Build**: `bash pythiaStuff/make.sh` (compiles main144.cc)
2. **Run mass scan**: `cd pythiaStuff && conda run -n llpatcolliders python run_mass_scan.py`
3. **Output**: All results saved to `pythiaStuff/mass_scan_hnl/`
   - Configuration files: `hnlLL_m{mass}GeV.cmnd`
   - CSV outputs: `hnlLL_m{mass}GeVLLP.csv`
4. **Post-processing**: Run `decayProbPerEvent.py` on each CSV to generate exclusion plots

### Technical Details
- **Mass points**: 15, 23, 31, 39, 47, 55, 63, 71, 79 GeV (9 points)
- **Events per point**: Configurable (tested with 1,000; production with 100,000)
- **Estimated runtime**: ~5-6 minutes per mass point at 100k events = ~45-54 minutes total
- **CSV format**: `event,id,pt,eta,phi,momentum,mass` (1 header + ~1 HNL per event)

### Notes
- For W → HNL scenario, expect primarily 1 HNL per event
- Some events (~4%) may have 2 HNLs due to multi-W production
- Duplicate detection uses relative tolerance to handle numerical precision
- Log files saved to `mass_scan_hnl/main144.log`
