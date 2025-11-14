# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

This is a physics research codebase for analyzing Long-Lived Particle (LLP) detection at collider experiments. It combines Monte Carlo simulation, geometric modeling, and statistical analysis to study particle decay probabilities in detector geometries.

### Repository Structure
```
.
├── pythiaStuff/          # PYTHIA 8 simulation code and configs
│   ├── main144.cc        # C++ simulation code
│   ├── make.sh           # Build script
│   ├── *.cmnd            # PYTHIA configuration files
│   └── *LLP.csv          # Output CSV files from simulation
├── external/             # Experimental limits data (ANUBIS, CODEX, MATHUSLA)
├── decayProbPerEvent.py  # Main post-simulation analysis script
├── neutral3D.py          # 3D geometric analysis
├── neutralv2.py          # 2D geometric analysis
└── environment.yml       # Conda environment specification
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

Particle CSV format: `event,id,pt,eta,phi,momentum,mass` (two particles per event)

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

**Running the simulation:**
```bash
./pythiaStuff/main144 -c pythiaStuff/higgsLL.cmnd   # Higgs → LLP pairs
# or
./pythiaStuff/main144 -c pythiaStuff/hnlLL.cmnd     # W → Heavy Neutral Leptons
```

Available configuration files:
- `higgsLL.cmnd`: Higgs production with decay to long-lived particle pairs
- `hnlLL.cmnd`: W boson production with decay to Heavy Neutral Leptons

This runs the compiled simulation with the specified `.cmnd` configuration file, outputting CSV files directly to `pythiaStuff/`.

### Running Analysis

**Post-Simulation Step** (run with conda environment activated):
```bash
conda activate llpatcolliders
python decayProbPerEvent.py pythiaStuff/higgsLLLLP.csv    # Main post-simulation analysis
# Calculates decay probabilities with lifetime scanning and generates exclusion plots
```

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

## Physics Context

The code models scenarios where:

* Long-lived particles (LLPs) are produced in various processes:
  - Higgs decays to LLP pairs (higgsLL.cmnd)
  - W boson decays to Heavy Neutral Leptons (hnlLL.cmnd)
* Particles travel distances before decaying in the detector
* Decay probability depends on particle lifetime and path length through detector
* Goal is to set exclusion limits on particle properties (decay length vs branching ratio)

Analysis compares calculated sensitivities with existing experimental limits from MATHUSLA, CODEX-b, and ANUBIS experiments.	
