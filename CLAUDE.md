# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

This is a physics research codebase for analyzing Long-Lived Particle (LLP) detection at collider experiments. It combines Monte Carlo simulation, geometric modeling, and statistical analysis to study particle decay probabilities in detector geometries.

## Core Architecture

### Simulation Pipeline
1. **Event Generation (PYTHIA 8)**  
   - Source & configs live in: `pythiaStuff/` (repo path: `llpatcolliders/tree/main/pythiaStuff`)  
   - Run command (from within `pythiaStuff/` or with correct working dir):  
     ```bash
     ./main144 -c higgsLL.cmnd
     ```  
   - Produces ROOT files with particle kinematics and decay information.

2. **Data Processing**  
   Jupyter notebooks (`Investigation.ipynb`, `Pythia_simulation.ipynb`) analyze ROOT simulation output.

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
# Direct execution of analysis scripts
python decayProbPerEvent.py    # Event-level decay analysis (generates the main "nice" output)
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

* Long-lived particles are produced in Higgs decays
* Particles travel distances before decaying
* Decay probability depends on particle lifetime and path length through detector
* Goal is to set exclusion limits on particle properties (mass vs decay length)

Analysis compares calculated sensitivities with existing experimental limits from MATHUSLA, CODEX-b, and ANUBIS experiments.	
