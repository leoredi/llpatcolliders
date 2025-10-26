# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

This is a physics research codebase for analyzing Long-Lived Particle (LLP) detection at collider experiments. It combines Monte Carlo simulation, geometric modeling, and statistical analysis to study particle decay probabilities in detector geometries.

**Two Physics Scenarios Supported:**
1. **Higgs → Scalar LLP Pairs** (Original) - Uses `higgsLL.cmnd`
2. **W → Heavy Neutral Leptons (HNLs)** (New) - Uses `hnlLL.cmnd` ← **See HNL_WORKFLOW.md for complete guide**

## Core Architecture

### Simulation Pipeline
1. **Event Generation (PYTHIA 8)**
   - Source & configs live in: `pythiaStuff/` (repo path: `llpatcolliders/tree/main/pythiaStuff`)
   - Run commands (from within `pythiaStuff/` or with correct working dir):
     ```bash
     # Higgs analysis (original)
     ./main144 -c higgsLL.cmnd

     # HNL analysis (new)
     ./main144 -c hnlLL.cmnd
     ```
   - Produces ROOT files with particle kinematics and decay information.

2. **Data Processing**
   - Jupyter notebooks (`Investigation.ipynb`, `Pythia_simulation.ipynb`) analyze ROOT simulation output
   - **HNL-specific**: `pythiaStuff/convertRootToCsv_HNL.py` extracts HNL particles (PDG ±9900014) from ROOT

3. **Geometric Analysis**
   Python scripts model detector geometries and calculate decay probabilities.

### Key Components

**Particle Data Flow**:
- PYTHIA 8 simulation → ROOT files → CSV particle data → Analysis scripts
- Particle data includes: event ID, particle ID, kinematics (pt, eta, phi, momentum, mass)

**Geometric Modeling**:
- `neutral3D.py`: Full 3D tube geometry with ray-casting for omnidirectional particle flux
- `neutralv2.py`: Simplified 2D detector geometry using Shapely
- Both calculate decay probabilities using exponential decay laws

**Analysis Scripts**:
- **`decayProbPerEvent.py`** *(main output/plots generator)*: Event-level decay probability calculations with lifetime scanning
  - Repo path: `llpatcolliders/blob/main/decayProbPerEvent.py` (file is at repo root as `decayProbPerEvent.py`)
  - Creates exclusion plots comparing with experimental limits from `external/` directory
  - **Updated for HNL analysis** with HNL-specific labels, titles, and output files

## Dependencies and Setup

### C++ Code (PYTHIA Simulation)
```bash
# Requires PYTHIA 8 and ROOT framework installed
# Simulation code & configs are under: pythiaStuff/
# Typical run (adjust working dir/environment as needed):
cd pythiaStuff
./main144 -c higgsLL.cmnd
````

### Python Analysis

```bash
# Core dependencies
pip install numpy pandas matplotlib scipy trimesh shapely tqdm

# Also requires ROOT with Python bindings (PyROOT)
# Installation varies by system - typically via conda-forge or system package manager
```

### Running Analysis

```bash
# HNL Analysis (NEW - Complete Workflow)
cd pythiaStuff
./main144 -c hnlLL.cmnd                        # Generate HNL events
python convertRootToCsv_HNL.py main144.root ../LLP.csv  # Convert to CSV
cd ..
python decayProbPerEvent.py                    # Generate HNL exclusion plots

# Original Higgs Analysis
cd pythiaStuff
./main144 -c higgsLL.cmnd                      # Generate Higgs events
# (use appropriate ROOT to CSV converter)
cd ..
python decayProbPerEvent.py                    # Generate exclusion plots

# Alternative geometric analyses
python neutral3D.py            # 3D geometric analysis
python neutralv2.py            # 2D geometric analysis

# Jupyter notebooks for ROOT data analysis
jupyter notebook Investigation.ipynb
jupyter notebook Pythia_simulation.ipynb
jupyter notebook Pythia_time_intersection_calculation.ipynb
```

## Key Data Structures

**Particle CSV Format**: `event,id,pt,eta,phi,momentum,mass`

* Two particles per event (particle/antiparticle pairs)
* Used by decay probability calculations

**Tube Geometry**:

* Defined by corrected vertex paths in 3D space at z=22m
* Radius \~1.54m, used for ray-detector intersections

**ROOT Tree Structure**:

* Branches: energy, x, y, z, t, pid, phi, theta, px, py, pz, MC\_event
* Includes mother/sister/daughter particle relationships

## Physics Context

The code models scenarios where:

**Higgs Scenario (Original)**:
* Long-lived scalar particles are produced in Higgs decays (H → LL + LL̄)
* Scalars decay to b-quarks, mass ~15 GeV

**HNL Scenario (New)**:
* Heavy Neutral Leptons produced from W boson decays (W → HNL + μ)
* HNLs decay to leptons + quarks, mass ~5 GeV (configurable 1-10 GeV)
* Represents mixing between HNLs and Standard Model neutrinos

**Common Features**:
* Long-lived particles travel macroscopic distances before decaying
* Decay probability depends on particle lifetime and path length through detector
* Goal is to set exclusion limits on particle properties (mass vs decay length/branching ratio)
* Analysis compares calculated sensitivities with existing experimental limits from MATHUSLA, CODEX-b, and ANUBIS experiments	
