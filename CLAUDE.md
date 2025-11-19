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

- Heavy Neutral Leptons (HNLs) from B-meson decays
- LLPs produced in proton-proton collisions at √s = 13.6 TeV
- Mass range scanning from low masses (1 GeV) to higher masses

## Simulation Framework

### Main Components

1. **Event Generation**: Uses PYTHIA 8.315 for pp collisions and B-meson production
2. **Mass Scan**: Systematic scan across HNL mass points
3. **Decay Analysis**: Tracking of HNL production and decay within detector acceptance
4. **Background Estimation**: Simulation of SM backgrounds and detector response

### Key Files

- `pythiaStuff/main_hnl_scan.cc`: Main C++ simulation code for HNL mass scanning
- `pythiaStuff/make.sh`: Compilation script for the simulation
- `pythiaStuff/hnl_LowMass_Inclusive_Complete.cmnd`: PYTHIA configuration for low-mass HNLs
- `pythiaStuff/hnl_HighMass_Inclusive_Complete.cmnd`: PYTHIA configuration for high-mass HNLs

### Simulation Strategy

The simulation performs:
1. Generation of b-bbar pairs from hard QCD processes
2. Hadronization and B-meson production
3. B-meson decays including HNL production channels
4. Tracking of HNL decay products
5. Geometric acceptance cuts for the detector location
6. Statistical analysis for sensitivity projections

## Analysis Approach

The analysis strategy focuses on:
- Design optimization of the detector geometry
- Background characterization and mitigation
- Signal efficiency calculations
- Sensitivity reach for various LLP benchmarks
- Comparison with existing LLP detector proposals

## Current Status

The simulation framework is operational and performs mass scans over the HNL parameter space to determine detector sensitivity and physics reach.
