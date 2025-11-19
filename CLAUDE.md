# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

This is a physics research codebase for analyzing Long-Lived Particle (LLP) detection at collider experiments. It combines Monte Carlo simulation, geometric modeling, and statistical analysis to study particle decay probabilities in detector geometries.

### Repository Structure
```
.
├── pythiaStuff/          # PYTHIA 8 simulation code and configs
│   ├── main144.cc        # C++ simulation code (with duplicate detection)
│   ├── main144           # Compiled executable
│   ├── make.sh           # Build script
│   ├── run_mass_scan.py  # Automated mass scan script for HNL (supports mu/tau)
│   ├── higgsLL.cmnd      # Higgs → LLP pairs configuration
│   ├── hnlLL.cmnd        # W → HNL (muon-coupled) configuration (template)
│   ├── hnlTauLL.cmnd     # W → HNL (tau-coupled) configuration (template)
│   ├── hnlTauLL_test.cmnd # Tau scenario test config (1k events, 31 GeV)
│   ├── hnlLL_m*GeV.cmnd  # Mass-specific muon HNL configurations (15-71 GeV)
│   └── hnlTauLL_m*GeV.cmnd # Mass-specific tau HNL configurations (15-71 GeV)
├── output/               # All output files organized by type
│   ├── csv/              # CSV files from simulation and analysis
│   └── images/           # PNG plots and visualizations
├── external/             # Experimental limits data
│   ├── ANUBIS.csv        # ANUBIS experiment exclusion limits
│   ├── CODEX.csv         # CODEX-b experiment exclusion limits
│   └── MATHUSLA.csv      # MATHUSLA experiment exclusion limits
├── decayProbPerEvent.py  # Main post-simulation analysis script
├── hnl_coupling_limit.py # HNL coupling limit calculator (BR→|U|² conversion)
├── quick_test_coupling.py # Quick test for coupling conversion logic
├── create_coupling_plot.sh # Automated coupling analysis pipeline
├── neutral3D.py          # 3D geometric analysis
├── neutralv2.py          # 2D geometric analysis
├── auto_analyze.sh       # Automated analysis script for all mass points
├── environment.yml       # Conda environment specification
├── CHANGELOG.md          # Development log and changes
├── CLAUDE.md             # AI assistant guidance (this file)
└── README.md             # Repository documentation
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
PYTHIA 8 Simulation (main144) → CSV particle data → decayProbPerEvent.py → Exclusion plots (BR vs lifetime)
                                                                        ↓
                                                    hnl_coupling_limit.py → Coupling plot (|U|² vs mass)
```

Particle CSV format: `event,id,pt,eta,phi,momentum,mass` (~1 HNL per event)

### Complete Analysis Pipeline Example
```bash
# 1. Build PYTHIA simulation
bash pythiaStuff/make.sh

# 2. Run mass scan (generates CSV files for all mass points)
cd pythiaStuff

# For muon-coupled HNL (default)
conda run -n llpatcolliders python run_mass_scan.py

# For tau-coupled HNL
conda run -n llpatcolliders python run_mass_scan.py --scenario tau

# 3. Analyze each mass point for detector sensitivity
conda activate llpatcolliders

# For muon-coupled results
for mass in 15 23 31 39 47 55 63 71; do
    python decayProbPerEvent.py output/csv/hnlLL_m${mass}GeVLLP.csv
done

# For tau-coupled results
for mass in 15 23 31 39 47 55 63 71; do
    python decayProbPerEvent.py output/csv/hnlTauLL_m${mass}GeVLLP.csv
done

# Output: Exclusion plots saved in output/images/

# 4. Generate coupling vs mass plot (final physics result)
python hnl_coupling_limit.py --scenario mu
# or use automated script:
bash create_coupling_plot.sh mu
```

### Coupling Limit Analysis Pipeline

After running `decayProbPerEvent.py` for all mass points, you can convert the BR vs lifetime exclusions into coupling limits (|U_ℓ|² vs mass) - the standard experimental result format for HNL searches.

**Physics:** For Heavy Neutral Leptons, production and decay are coupled through mixing:
- Production: BR(W → ℓ N) ∝ |U_ℓ|² × f(m_N)  [phase space factor]
- Decay: τ_N ∝ 1/|U_ℓ|²  [lifetime inversely proportional to coupling]

**Scripts:**
- `hnl_coupling_limit.py`: Main conversion script
  - Implements phase space factor f(m_N, m_W, m_ℓ) for W → ℓ N
  - Calculates HNL lifetime using Γ_N ≈ C × |U_ℓ|² × G_F² × m_N⁵ (C ≈ 1.7×10⁻³)
  - Converts each mass point's (BR, τ) exclusion to |U_ℓ|² limit
  - Creates coupling vs mass plot with all mass points
- `quick_test_coupling.py`: Test with synthetic data
- `create_coupling_plot.sh`: Automated pipeline (checks for missing data, runs analysis, generates plot)

**Usage:**
```bash
# After running decayProbPerEvent.py for all mass points:

# Generate coupling plot for muon scenario
python hnl_coupling_limit.py --scenario mu

# Generate coupling plot for tau scenario
python hnl_coupling_limit.py --scenario tau

# Test mode (uses synthetic data)
python hnl_coupling_limit.py --scenario mu --test

# Or use automated script that runs missing analyses first
bash create_coupling_plot.sh mu
```

**Outputs:**
- `output/images/hnl_coupling_vs_mass_mu.png` - Muon-coupled HNL sensitivity
- `output/images/hnl_coupling_vs_mass_tau.png` - Tau-coupled HNL sensitivity
- Log-log plot showing |U_ℓ|² limits from 10⁻¹⁰ to 10⁻² over mass range 10-100 GeV

**Requirements:**
- Needs `{scenario}_m{mass}GeVLLP_exclusion_data.csv` for each mass point
- These are created automatically by `decayProbPerEvent.py` (added in 2025-11-18)
- If missing, script will print which analyses need to be run

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

# Muon-coupled HNL (default)
conda run -n llpatcolliders python run_mass_scan.py

# Tau-coupled HNL
conda run -n llpatcolliders python run_mass_scan.py --scenario tau
```

This will:
- Generate configuration files for mass points from 15-71 GeV (8 points in 8 GeV steps)
- Run PYTHIA for each mass point **in parallel (2 at a time)**
- Print timestamped progress showing which masses are currently running
- Save all CSV outputs to `output/csv/` directory
- Default configuration: **1M events per mass point** (N_EVENTS = 1_000_000 in run_mass_scan.py)
- Runtime: ~3.5-4.5 hours for full scan at 1M events (with 2x parallelization speedup)

Available configuration files:
- `higgsLL.cmnd`: Higgs production with decay to long-lived particle pairs
- `hnlLL.cmnd`: W boson production with decay to Heavy Neutral Leptons (muon-coupled)
- `hnlTauLL.cmnd`: W boson production with decay to Heavy Neutral Leptons (tau-coupled)
- `hnlTauLL_test.cmnd`: Quick tau scenario test (1k events, 31 GeV)
- `run_mass_scan.py`: Automated mass scanning (configurable via N_EVENTS variable and --scenario flag)

**Note on CSV output**: The code includes duplicate detection to ensure only unique particles are written. For W → HNL, expect ~1 HNL per event (some events may have 2 due to multi-W production).

### Running Analysis

**Post-Simulation Analysis** (run with conda environment activated):
```bash
conda activate llpatcolliders

# For a single mass point (muon scenario)
python decayProbPerEvent.py output/csv/hnlLL_m15GeVLLP.csv

# For all muon-coupled mass points (can be run in parallel)
for mass in 15 23 31 39 47 55 63 71; do
    python decayProbPerEvent.py output/csv/hnlLL_m${mass}GeVLLP.csv
done

# For all tau-coupled mass points (can be run in parallel)
for mass in 15 23 31 39 47 55 63 71; do
    python decayProbPerEvent.py output/csv/hnlTauLL_m${mass}GeVLLP.csv
done
```

**Analysis outputs:**

For muon-coupled HNL:
- `output/images/hnlLL_m{mass}GeVLLP_exclusion_vs_lifetime.png` - Main exclusion plot comparing detector sensitivity with MATHUSLA, CODEX-b, and ANUBIS
- `output/images/hnlLL_m{mass}GeVLLP_correlation_analysis.png` - Event correlation and probability distribution plots
- `output/csv/hnlLL_m{mass}GeVLLP_event_decay_statistics.csv` - Event-level decay statistics
- `output/csv/hnlLL_m{mass}GeVLLP_particle_decay_results.csv` - Particle-level results with path lengths and decay probabilities

For tau-coupled HNL:
- `output/images/hnlTauLL_m{mass}GeVLLP_exclusion_vs_lifetime.png` - Main exclusion plot
- `output/images/hnlTauLL_m{mass}GeVLLP_correlation_analysis.png` - Event correlation plots
- `output/csv/hnlTauLL_m{mass}GeVLLP_event_decay_statistics.csv` - Event-level statistics
- `output/csv/hnlTauLL_m{mass}GeVLLP_particle_decay_results.csv` - Particle-level results

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
- **Event Structure**: For W → HNL scenario, expect ~1 HNL per event (CSV has ~1M particles for 1M events)
- **Mass Scanning**: Use `run_mass_scan.py` to automate generation of multiple mass points for exclusion limit studies
  - Currently configured for 1M events per mass point (modify N_EVENTS if different statistics needed)
  - Runs 2 simulations in parallel for 2x speedup
  - Prints real-time progress with timestamps showing which masses are running
- **CSV Format**: One particle per line with format `event,id,pt,eta,phi,momentum,mass`

## Physics Context

The code models scenarios where:

* Long-lived particles (LLPs) are produced in various processes:
  - Higgs decays to LLP pairs (higgsLL.cmnd)
  - W boson decays to Heavy Neutral Leptons:
    - Muon-coupled: W → μ N, N → μ± + jets (hnlLL.cmnd)
    - Tau-coupled: W → τ N, N → τ± + jets (hnlTauLL.cmnd)
* Particles travel distances before decaying in the detector
* Decay probability depends on particle lifetime and path length through detector
* Goal is to set exclusion limits on particle properties (decay length vs branching ratio)

**Scenario comparison:**
- Both muon and tau scenarios use the same HNL PDG ID (9900015) and physics model
- Production differs only in the lepton coupling (PDG 13/-13 for muon, 15/-15 for tau)
- HNL decay channels mirror production: N → μ±/τ± + jets (same quark content: u d̄, ū d)
- Analysis framework is identical for both scenarios

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
- Current default: N_EVENTS = 1_000_000 (modify in `run_mass_scan.py` if different statistics needed for testing)

**Build errors**:
- Ensure PYTHIA 8 is properly installed and accessible
- Check that `make.sh` has correct paths to PYTHIA installation
- Verify `#include <set>` is present in `main144.cc` for duplicate detection

For detailed change history, see `CHANGELOG.md`.	
