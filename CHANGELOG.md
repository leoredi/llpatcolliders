# Changelog

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
