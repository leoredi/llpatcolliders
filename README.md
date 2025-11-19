# llpatcolliders

Monte Carlo simulation and analysis framework for Long-Lived Particle (LLP) detection at collider experiments. This codebase studies Heavy Neutral Lepton (HNL) production and decay at the HL-LHC, calculating detector sensitivity and exclusion limits.

## Overview

This repository implements a complete physics analysis pipeline:
1. **PYTHIA 8 simulation** - Generate HNL particles from W boson decays
2. **Geometric analysis** - Calculate decay probabilities in detector geometry
3. **Statistical analysis** - Compute exclusion limits as function of lifetime/branching ratio
4. **Coupling limits** - Convert to standard |U_ℓ|² vs mass format for experimental comparison

## Quick Start

```bash
# Setup conda environment
conda env create -f environment.yml
conda activate llpatcolliders

# Build PYTHIA simulation
bash pythiaStuff/make.sh

# Run full analysis for muon-coupled HNL
cd pythiaStuff
conda run -n llpatcolliders python run_mass_scan.py --scenario mu
cd ..
bash create_coupling_plot.sh mu
```

This will:
- Generate 1M events for each of 8 mass points (15-71 GeV)
- Analyze detector sensitivity vs lifetime for each mass
- Create final coupling limit plot: |U_μ|² vs mass

## Repository Structure

```
.
├── pythiaStuff/              # PYTHIA 8 simulation
│   ├── main144.cc            # C++ event generator
│   ├── run_mass_scan.py      # Automated mass scanning
│   ├── hnlLL.cmnd            # Muon-coupled HNL config
│   └── hnlTauLL.cmnd         # Tau-coupled HNL config
├── output/
│   ├── csv/                  # Simulation and analysis data
│   └── images/               # Plots and visualizations
├── external/                 # Experimental limit data (MATHUSLA, CODEX-b, ANUBIS)
├── decayProbPerEvent.py      # Main analysis script (BR vs lifetime)
├── hnl_coupling_limit.py     # Coupling limit calculator (|U|² vs mass)
├── create_coupling_plot.sh   # Automated pipeline script
└── environment.yml           # Conda dependencies
```

## Physics Scenarios

### Heavy Neutral Leptons (HNL)

**Muon-coupled:** W → μ N, with N → μ± + jets
**Tau-coupled:** W → τ N, with N → τ± + jets

Both scenarios use the same HNL PDG ID (9900015) and mass grid, differing only in lepton coupling.

## Detailed Workflow

### 1. Simulation (PYTHIA 8)

Generate HNL particles at different mass points:

```bash
cd pythiaStuff

# Muon scenario (default)
conda run -n llpatcolliders python run_mass_scan.py

# Tau scenario
conda run -n llpatcolliders python run_mass_scan.py --scenario tau
```

**Configuration:**
- Mass points: 15, 23, 31, 39, 47, 55, 63, 71 GeV (8 points)
- Events per mass: 1M (configurable via `N_EVENTS` in script)
- Parallel execution: 2 simulations at once
- Runtime: ~3.5-4.5 hours for full scan

**Output:** CSV files with format `event,id,pt,eta,phi,momentum,mass`
- Muon: `output/csv/hnlLL_m{mass}GeVLLP.csv`
- Tau: `output/csv/hnlTauLL_m{mass}GeVLLP.csv`

### 2. Detector Sensitivity Analysis

Calculate decay probabilities in detector geometry:

```bash
# For each mass point
python decayProbPerEvent.py output/csv/hnlLL_m15GeVLLP.csv

# Or batch process
for mass in 15 23 31 39 47 55 63 71; do
    python decayProbPerEvent.py output/csv/hnlLL_m${mass}GeVLLP.csv
done
```

**Analysis:**
- Detector: Tube at z=22m, radius ~1.54m
- Scans over 20 lifetime values (log-spaced)
- Calculates event-level and particle-level decay probabilities
- Compares with MATHUSLA, CODEX-b, ANUBIS experiments

**Outputs per mass point:**
- `{name}_exclusion_vs_lifetime.png` - Main exclusion plot
- `{name}_correlation_analysis.png` - Event correlation plots
- `{name}_exclusion_data.csv` - BR vs lifetime limits (for coupling conversion)
- `{name}_event_decay_statistics.csv` - Event-level statistics
- `{name}_particle_decay_results.csv` - Particle-level results

### 3. Coupling Limit Calculation

Convert BR vs lifetime to |U_ℓ|² vs mass (standard HNL format):

```bash
# Generate coupling plot
python hnl_coupling_limit.py --scenario mu

# Or use automated script (runs missing analyses first)
bash create_coupling_plot.sh mu
```

**Physics conversion:**
- BR(W → ℓ N) = |U_ℓ|² × f(m_N) / 9  (f = phase space factor)
- τ_N = ℏ / (C × |U_ℓ|² × G_F² × m_N⁵)  (C ≈ 1.7×10⁻³)
- Exclusion: |U_ℓ|² excluded if BR(|U|²) > BR_limit(τ(|U|²))

**Output:**
- `output/images/hnl_coupling_vs_mass_mu.png` - Muon-coupled sensitivity
- `output/images/hnl_coupling_vs_mass_tau.png` - Tau-coupled sensitivity

## Key Features

### Duplicate Detection
The PYTHIA simulation includes sophisticated duplicate filtering to handle the same particle appearing at multiple event record stages. Uses charge-based tracking (`std::set<int>`) to ensure each particle is written only once.

### Mass Scanning
Automated mass grid generation with:
- Parallel execution (2× speedup)
- Real-time progress reporting
- Automatic configuration file creation
- Support for both muon and tau scenarios

### Geometric Modeling
3D ray-tube intersection calculations for accurate path length determination:
- Curved detector geometry
- Event-by-event particle tracking
- Lifetime-dependent decay probability

### Statistical Framework
Proper event-level analysis accounting for:
- Multiple particles per event
- Correlated decays
- Background-free assumption
- Poisson statistics for 3 signal events

## Analysis Outputs

### Exclusion Plots (BR vs Lifetime)
Shows detector sensitivity compared to:
- MATHUSLA (surface detector)
- CODEX-b (near LHCb interaction point)
- ANUBIS (ATLAS side shaft)

### Coupling Plots (|U|² vs Mass)
Standard experimental result format showing:
- Mixing parameter limits across mass range
- Comparison with existing experiments (future)
- Both muon and tau coupling scenarios

## Requirements

- **PYTHIA 8** - Event generation (no ROOT dependency)
- **Python 3.11** - Analysis scripts
- **Conda** - Environment management

Python packages (in `environment.yml`):
- Scientific: numpy, pandas, scipy, matplotlib
- Geometry: shapely, trimesh, rtree
- Utilities: tqdm, jupyter

## Troubleshooting

**Empty CSV after simulation:**
- Check `Main:writeRoot = on` in `.cmnd` file
- Verify HNL PDG ID (9900015) matches in config and `main144.cc`

**Missing exclusion data for coupling plot:**
- Run `decayProbPerEvent.py` for all mass points first
- Or use `bash create_coupling_plot.sh {scenario}` to automate

**Build errors:**
- Ensure PYTHIA 8 is installed and paths in `make.sh` are correct
- Check `#include <set>` present in `main144.cc`

**Mass scan issues:**
- Must run from `pythiaStuff/`: `cd pythiaStuff && conda run -n llpatcolliders python run_mass_scan.py`
- Default is 1M events per mass (modify `N_EVENTS` for testing)

## Documentation

- **CLAUDE.md** - Detailed guide for AI assistant (full technical documentation)
- **CHANGELOG.md** - Development history and feature additions
- **README.md** - This file (user-facing documentation)

## Citation

This analysis framework was developed for studying HNL detection at the HL-LHC using the milliQan detector geometry.

## License

Research code - see repository for details.
