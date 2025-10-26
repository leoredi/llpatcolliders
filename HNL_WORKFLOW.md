# Heavy Neutral Lepton (HNL) Analysis Workflow

This guide describes how to run the complete HNL analysis pipeline, from simulation to exclusion plots.

## Overview

The HNL analysis pipeline consists of three main steps:

1. **Generate HNL events with PYTHIA 8** - Simulate pp collisions producing W bosons that decay to HNLs
2. **Convert ROOT to CSV** - Extract HNL particle kinematics from ROOT output
3. **Analyze decay probabilities** - Calculate exclusion limits based on detector geometry

## Prerequisites

### Software Requirements
- PYTHIA 8 (with ROOT output support)
- ROOT framework with Python bindings (PyROOT)
- Python 3.x with packages: numpy, pandas, matplotlib, scipy, trimesh, shapely, tqdm

### File Structure
```
llpatcolliders/
├── pythiaStuff/
│   ├── hnlLL.cmnd              # HNL PYTHIA configuration (NEW)
│   ├── main144                 # PYTHIA executable
│   ├── convertRootToCsv_HNL.py # ROOT to CSV converter (NEW)
│   └── make.sh                 # Build script
├── decayProbPerEvent.py        # Main analysis script (UPDATED for HNL)
├── external/                   # Experimental limits
│   ├── MATHUSLA.csv
│   ├── CODEX.csv
│   └── ANUBIS.csv
└── HNL_WORKFLOW.md            # This file
```

## Step-by-Step Instructions

### Step 1: Generate HNL Events

Navigate to the pythiaStuff directory and run PYTHIA with the HNL configuration:

```bash
cd pythiaStuff
./main144 -c hnlLL.cmnd
```

This will:
- Generate pp collision events at √s = 13.6 TeV
- Produce W bosons that decay to HNL + lepton
- Output results to `main144.root`
- Create a log file with event information

**Configuration options** (edit `hnlLL.cmnd`):
- `9900014:m0` - HNL mass (default: 5 GeV, range: 1-10 GeV)
- `9900014:tau0` - HNL proper lifetime (will be scanned in analysis)
- Decay channels and branching ratios for HNL → lepton + quarks

### Step 2: Convert ROOT to CSV

Extract HNL particle data from the ROOT file:

```bash
python convertRootToCsv_HNL.py main144.root ../LLP.csv
```

This creates a CSV file with columns:
- `event` - Event number
- `id` - PDG particle ID (9900014 for HNL, -9900014 for anti-HNL)
- `pt` - Transverse momentum [GeV/c]
- `eta` - Pseudorapidity
- `phi` - Azimuthal angle [radians]
- `momentum` - Total momentum [GeV/c]
- `mass` - Particle mass [GeV/c²]

**Note**: The script extracts only HNL particles (PDG ID ±9900014). Each event typically contains 2 HNLs (particle + antiparticle).

### Step 3: Run Decay Probability Analysis

Return to the main directory and run the analysis:

```bash
cd ..
python decayProbPerEvent.py
```

This will:
1. Load HNL particle data from `LLP.csv`
2. Calculate decay probabilities for a single lifetime (10 ns default)
3. Scan over lifetime parameter space (0.1 ns to 1 μs)
4. Generate exclusion plots comparing with experimental limits
5. Save output plots and statistics

**Outputs:**
- `HNL_exclusion_vs_lifetime.png` - Main exclusion plot with experimental comparisons
- `HNL_particle_decay_results.csv` - Particle-level decay probabilities
- `HNL_event_decay_statistics.csv` - Event-level statistics

## Physics Details

### HNL Production Mechanism

In this analysis, HNLs are produced via W boson decays:
```
pp → W± + X
W+ → HNL + μ+
W- → anti-HNL + μ-
```

This represents mixing between the HNL and Standard Model neutrinos.

### HNL Decay Channels

The HNL decays to leptons and quarks:
- HNL → μ⁻ + d + d̄  (33%)
- HNL → μ⁻ + s + s̄  (33%)
- HNL → μ⁻ + u + d̄  (34%)

### Decay Probability Calculation

For each HNL, the code calculates:

1. **Decay length**: `λ = γβcτ`
   - γ = Lorentz boost factor
   - β = velocity as fraction of c
   - τ = proper lifetime

2. **Survival probability to detector**: `P(survive to d) = exp(-d/λ)`

3. **Decay probability in detector**: `P(decay in L) = exp(-d_entry/λ) × [1 - exp(-L/λ)]`
   - d_entry = distance to detector entry
   - L = path length through detector

4. **Event-level probability**: For events with 2 HNLs:
   - P(at least one decays) = 1 - (1-p₁)(1-p₂)
   - P(both decay) = p₁ × p₂

### Detector Geometry

The analysis uses a cylindrical tube detector at z = 22 m with:
- Radius: ~1.54 m
- Complex 2D path following a curved trajectory
- Ray-tracing to determine HNL intersections with detector volume

## Customization

### Adjusting HNL Mass

Edit `pythiaStuff/hnlLL.cmnd`:
```
9900014:m0 = 3.0  # Change to desired mass in GeV
```
Then re-run steps 1-3.

### Changing Lifetime Scan Range

Edit `decayProbPerEvent.py` line ~482:
```python
lifetimes = np.logspace(-9.5, -4.5, 20)  # Adjust exponent range
```

### Modifying Production Mechanism

Alternative HNL production modes can be implemented in `hnlLL.cmnd`:
- Meson decays (K, D, B → HNL + ...)
- Direct Drell-Yan production
- Higgs decays (for different physics scenario)

## Troubleshooting

**Problem**: `main144.root` not created
- **Solution**: Check PYTHIA installation and that ROOT output is enabled in build

**Problem**: No HNLs found in ROOT file
- **Solution**: Verify PDG codes match in both `hnlLL.cmnd` and `convertRootToCsv_HNL.py`

**Problem**: Analysis crashes with "file not found"
- **Solution**: Ensure `LLP.csv` exists and paths in scripts are correct

**Problem**: Empty plots or zero probabilities
- **Solution**: Check that HNLs actually hit the detector tube (eta/phi coverage)

## Comparing with Higgs Analysis

The original Higgs analysis can still be run:
```bash
cd pythiaStuff
./main144 -c higgsLL.cmnd  # Use Higgs config instead
```

Key differences:
- **Production**: Higgs → scalar LLP pairs vs W → HNL + lepton
- **Mass**: ~15 GeV (scalar) vs ~5 GeV (HNL)
- **Decay modes**: b-quarks vs leptons + quarks
- **Experimental limits**: Different interpretations

## References

- PYTHIA 8 manual: https://pythia.org/
- HNL phenomenology: arXiv:1805.08567, arXiv:1909.12524
- MATHUSLA: arXiv:1811.00927
- CODEX-b: arXiv:1911.00459
- ANUBIS: arXiv:1909.13022

## Next Steps

1. **Systematic studies**: Vary HNL mass, scan production mechanisms
2. **Detector optimization**: Test different geometries, positions
3. **Background estimation**: Add SM backgrounds to simulation
4. **Statistical analysis**: Implement proper limit-setting procedures
