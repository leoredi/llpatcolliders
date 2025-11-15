# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

This is a physics research codebase for analyzing Long-Lived Particle (LLP) detection at collider experiments. It combines Monte Carlo simulation, geometric modeling, and statistical analysis to study particle decay probabilities in detector geometries.

### Repository Structure
```
.
├── pythiaStuff/          # PYTHIA 8 simulation code and configs
│   ├── main144.cc        # C++ simulation code (with duplicate detection)
│   ├── make.sh           # Build script
│   ├── run_mass_scan.py  # Automated mass scan script for HNL
│   ├── mass_scan_hnl/    # Output directory for mass scan results
│   ├── *.cmnd            # PYTHIA configuration files
│   └── *LLP.csv          # Output CSV files from simulation
├── external/             # Experimental limits data (ANUBIS, CODEX, MATHUSLA)
├── decayProbPerEvent.py  # Main post-simulation analysis script
├── neutral3D.py          # 3D geometric analysis
├── neutralv2.py          # 2D geometric analysis
├── environment.yml       # Conda environment specification
└── CHANGELOG.md          # Development log and changes
```

## Core Architecture

### Simulation Pipeline
1. **Event Generation (PYTHIA 8)**
   - Source & configs live in: `pythiaStuff/` (repo path: `llpatcolliders/tree/main/pythiaStuff`)
   - Build command (from repo root):
     ```bash
     bash pythiaStuff/make.sh
     ```
     This compiles `pythiaStuff/main144.cc` using the locally installed PYTHIA 8.
   - Run command:
     ```bash
     ./pythiaStuff/main144 -c pythiaStuff/higgsLL.cmnd
     ```
   - Takes `.cmnd` configuration files as input and produces **CSV files** directly with particle kinematics and decay information.

2. **Post-Simulation Analysis**
   - **`decayProbPerEvent.py`** is the main analysis script run after simulation
   - Calculates event-level decay probabilities with lifetime scanning
   - Creates exclusion plots comparing with experimental limits from `external/` directory

3. **Additional Analysis**
   - `neutral3D.py`: Full 3D tube geometry with ray-casting for omnidirectional particle flux
   - `neutralv2.py`: Simplified 2D detector geometry using Shapely

### Workflow Summary
```
PYTHIA 8 Simulation (main144) → CSV particle data → decayProbPerEvent.py → Exclusion plots
```

Particle CSV format: `event,id,pt,eta,phi,momentum,mass` (~1 HNL per event)

### Complete Analysis Pipeline Example
```bash
# 1. Build PYTHIA simulation
bash pythiaStuff/make.sh

# 2. Run mass scan (generates CSV files for all mass points)
cd pythiaStuff && conda run -n llpatcolliders python run_mass_scan.py

# 3. Analyze each mass point for detector sensitivity
conda activate llpatcolliders
for mass in 15 23 31 39 47 55 63 71 79; do
    python decayProbPerEvent.py pythiaStuff/mass_scan_hnl/hnlLL_m${mass}GeVLLP.csv
done

# Output: Exclusion plots saved as hnlLL_m{mass}GeVLLP_exclusion_vs_lifetime.png
```

## Dependencies and Setup

### Environment Setup
This project uses a conda environment defined in `environment.yml`. Set up the environment with:
```bash
conda env create -f environment.yml
conda activate llpatcolliders
```

The conda environment includes:
- Python 3.11
- Scientific computing: numpy, pandas, scipy, matplotlib
- Geometry tools: shapely, trimesh, rtree
- Utilities: tqdm, jupyter

### C++ Code (PYTHIA Simulation)
Requires PYTHIA 8 installed locally (no ROOT dependency).

**Building the simulation:**
```bash
# From repo root
bash pythiaStuff/make.sh
```

**Running single simulations:**
```bash
./pythiaStuff/main144 -c pythiaStuff/higgsLL.cmnd   # Higgs → LLP pairs
# or
./pythiaStuff/main144 -c pythiaStuff/hnlLL.cmnd     # W → Heavy Neutral Leptons
```

**Running automated mass scans:**
```bash
cd pythiaStuff
conda run -n llpatcolliders python run_mass_scan.py
```

This will:
- Generate configuration files for mass points from 15-79 GeV (9 points in 8 GeV steps)
- Run PYTHIA for each mass point **in parallel (2 at a time)**
- Print timestamped progress showing which masses are currently running
- Save all outputs to `mass_scan_hnl/` subdirectory
- Take ~22-27 minutes for 100k events per mass point (with 2x parallelization speedup)

Available configuration files:
- `higgsLL.cmnd`: Higgs production with decay to long-lived particle pairs
- `hnlLL.cmnd`: W boson production with decay to Heavy Neutral Leptons
- `run_mass_scan.py`: Automated mass scanning (configurable via N_EVENTS variable)

**Note on CSV output**: The code includes duplicate detection to ensure only unique particles are written. For W → HNL, expect ~1 HNL per event (some events may have 2 due to multi-W production).

### Running Analysis

**Post-Simulation Analysis** (run with conda environment activated):
```bash
conda activate llpatcolliders

# For a single mass point
python decayProbPerEvent.py pythiaStuff/mass_scan_hnl/hnlLL_m15GeVLLP.csv

# For all mass points (can be run in parallel)
for mass in 15 23 31 39 47 55 63 71 79; do
    python decayProbPerEvent.py pythiaStuff/mass_scan_hnl/hnlLL_m${mass}GeVLLP.csv
done
```

**Analysis outputs:**
- `hnlLL_m{mass}GeVLLP_exclusion_vs_lifetime.png` - Main exclusion plot comparing detector sensitivity with MATHUSLA, CODEX-b, and ANUBIS
- `hnlLL_m{mass}GeVLLP_correlation_analysis.png` - Event correlation and probability distribution plots
- `hnlLL_m{mass}GeVLLP_event_decay_statistics.csv` - Event-level decay statistics
- `hnlLL_m{mass}GeVLLP_particle_decay_results.csv` - Particle-level results with path lengths and decay probabilities

**Additional Analysis Scripts:**
```bash
python neutral3D.py            # 3D geometric analysis
python neutralv2.py            # 2D geometric analysis
```

## Key Details

**Detector Geometry**:
- Tube detector at z=22m with radius ~1.54m (1.4 × 1.1 m)
- Follows a curved path defined by correctedVert coordinates
- Geometry implemented in 3D using trimesh for ray-tube intersection calculations
- Used for calculating decay probabilities based on path length through detector

**Important Implementation Notes**:
- **Duplicate Detection**: `main144.cc` includes kinematic-based duplicate filtering to prevent writing the same particle multiple times from different stages of the PYTHIA event record
- **Event Structure**: For W → HNL scenario, expect ~1 HNL per event (CSV has ~1000 particles for 1000 events)
- **Mass Scanning**: Use `run_mass_scan.py` to automate generation of multiple mass points for exclusion limit studies
  - Runs 2 simulations in parallel for 2x speedup
  - Prints real-time progress with timestamps showing which masses are running
- **CSV Format**: One particle per line with format `event,id,pt,eta,phi,momentum,mass`

## Physics Context

The code models scenarios where:

* Long-lived particles (LLPs) are produced in various processes:
  - Higgs decays to LLP pairs (higgsLL.cmnd)
  - W boson decays to Heavy Neutral Leptons (hnlLL.cmnd)
* Particles travel distances before decaying in the detector
* Decay probability depends on particle lifetime and path length through detector
* Goal is to set exclusion limits on particle properties (decay length vs branching ratio)

Analysis compares calculated sensitivities with existing experimental limits from MATHUSLA, CODEX-b, and ANUBIS experiments.

## Troubleshooting

**Duplicate particles in CSV output**:
- Fixed in current version of `main144.cc` with kinematic-based duplicate detection
- If you see excessive duplicates, rebuild with `bash pythiaStuff/make.sh`

**Empty CSV files after simulation**:
- Check that `Main:writeRoot = on` in your `.cmnd` file
- Verify the LLP PDG ID matches in both `.cmnd` and `main144.cc` (default: 9900015 for HNL)

**Mass scan script issues**:
- Must run from `pythiaStuff/` directory: `cd pythiaStuff && conda run -n llpatcolliders python run_mass_scan.py`
- Check that conda environment is activated: `conda activate llpatcolliders`
- Configure N_EVENTS in `run_mass_scan.py` before running

**Build errors**:
- Ensure PYTHIA 8 is properly installed and accessible
- Check that `make.sh` has correct paths to PYTHIA installation
- Verify `#include <set>` is present in `main144.cc` for duplicate detection

For detailed change history, see `CHANGELOG.md`.	
