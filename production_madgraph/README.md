# MadGraph5 Production Pipeline for HNLs

This directory contains the pipeline for generating Heavy Neutral Lepton (HNL) events using MadGraph5_aMC@NLO.
This pipeline focuses on the **Electroweak (EW) regime** ($m_{HNL} \gtrsim 5$ GeV), where HNLs are produced via $W$ and $Z$ boson decays ($pp \to W/Z \to \ell N$).

## Motivation

The existing Pythia-only pipeline works well for meson-mediated production ($m_{HNL} < 5$ GeV). However, for heavier HNLs, Pythia struggles with off-shell W/Z bosons when forcing specific decay channels (kinematically forbidden if the W is off-shell with mass $< m_{HNL} + m_\ell$).

MadGraph5_aMC@NLO correctly handles:
- Off-shell W/Z contributions.
- Full matrix element calculations.
- Interference effects.

## Prerequisites

- **Python 3**
- **MadGraph5_aMC@NLO** installed.
  - Set the `MG5_PATH` environment variable to your MadGraph installation directory, or ensure the script defaults point to a valid location.
- **HNL UFO Model**: `SM_HeavyN_CKM_AllMasses_LO` (must be available to MadGraph).

## Directory Structure

- `cards/`: Templates for MadGraph process, run, and parameter cards.
- `scripts/`: Python scripts to drive the generation and analysis.
- `lhe/`: Generated LHE (Les Houches Event) files, organized by flavor and mass.
- `csv/`: Output CSV files converted from LHE, compatible with the project's analysis pipeline.

## Usage

### Running the Scan

The main driver script is `scripts/run_hnl_scan.py`. It automates the generation of cards, running MadGraph, and converting outputs.

```bash
# Run the full grid
python production_madgraph/scripts/run_hnl_scan.py

# Run a quick test (1 mass point, fewer events)
python production_madgraph/scripts/run_hnl_scan.py --test
```

### Configuration

You can configure the scan in `scripts/run_hnl_scan.py` (mass grid, flavors, events per point) or via environment variables.

## Output

- **LHE Files**: Stored in `production_madgraph/lhe/<flavour>/m_<mass>GeV/`.
- **CSV Files**: Per-event data stored in `production_madgraph/csv/<flavour>/`.
- **Summary**: A summary of cross-sections is appended to `production_madgraph/csv/summary_HNL_EW_production.csv`.
