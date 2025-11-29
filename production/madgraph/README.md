# MadGraph HNL Production Pipeline

**Electroweak Heavy Neutral Lepton (HNL) production at LHC 14 TeV**

This pipeline generates HNL events via W/Z boson decays using MadGraph5_aMC@NLO with the HeavyN UFO model.

---

## Table of Contents

1. [Why MadGraph?](#why-madgraph)
2. [Quick Start](#quick-start)
3. [Installation](#installation)
4. [Pipeline Overview](#pipeline-overview)
5. [Usage](#usage)
6. [Output Files](#output-files)
7. [Physics Details](#physics-details)
8. [Troubleshooting](#troubleshooting)

---

## Why MadGraph?

### The Pythia Problem

The existing Pythia-based EW HNL production (`production/main_hnl_production.cc`) attempts to generate W/Z → ℓ N by:
- Forcing W/Z production via Drell-Yan
- Hacking W/Z decay tables to enable W → ℓ N and Z → ν N
- Setting BR(W/Z → ℓN) = 1.0 for efficient sampling

**This fails for HNL masses ≳ 13 GeV because:**

1. **Off-shell W production**: At LHC, many W bosons are produced off-shell (W* with mass < m_W)
2. **Kinematic blocking**: When m(W*) < m_HNL + m_ℓ, the decay W* → ℓ N is kinematically forbidden
3. **No fallback channels**: Since we set BR(W → all other modes) = 0, Pythia has *no open decay channels*
4. **Pythia aborts**: "Error: No allowed decay modes for particle W+"

### The MadGraph Solution

State-of-the-art LLP studies (MATHUSLA, CODEX-b, ANUBIS) use **MadGraph** because:

✅ **Matrix element generation**: MadGraph computes pp → ℓ N via W/Z exchange at the matrix element level
✅ **Off-shell handling**: Automatically handles W*/Z* propagators with proper kinematics
✅ **Phase space**: Correctly samples full phase space including threshold regions
✅ **PDFs**: Uses proper parton distribution functions (PDFs) for LHC
✅ **No decay table hacks**: Process is pp → ℓ N, not pp → W → ℓ N (W is internal propagator)

**Result**: Works seamlessly from threshold (m_HNL ~ 5 GeV) to high mass (80+ GeV)

---

## Quick Start

### Test Run (Single Point)

```bash
# From project root
cd production/madgraph

# Test with 15 GeV muon (1000 events, ~5-10 minutes)
conda run -n mg5env python scripts/run_hnl_scan.py --test

# Check output
ls csv/muon/
cat csv/summary_HNL_EW_production.csv
```

### Full Production Scan

```bash
# All flavours, all masses (5-80 GeV, 50k events each)
# WARNING: This takes ~10-20 hours
conda run -n mg5env python scripts/run_hnl_scan.py

# Single flavour (e.g. muon only)
conda run -n mg5env python scripts/run_hnl_scan.py --flavour muon

# Custom mass points
conda run -n mg5env python scripts/run_hnl_scan.py --masses 10 15 20 25 --nevents 10000
```

---

## Installation

### 1. MadGraph Installation

MadGraph is already installed at `production/madgraph/mg5/`.

If you need to reinstall or update:

```bash
cd production/madgraph

# Download MadGraph 3.6.6 (or latest)
wget https://launchpad.net/mg5amcnlo/3.0/3.6.x/+download/MG5_aMC_v3.6.6.tar.gz
tar -xzf MG5_aMC_v3.6.6.tar.gz
mv MG5_aMC_v3_6_6 mg5
```

### 2. Install HeavyN UFO Model

The HeavyN model (`SM_HeavyN_CKM_AllMasses_LO`) is required for HNL processes.

**Option A: Install via MadGraph CLI**

```bash
conda run -n mg5env mg5/bin/mg5_aMC

# In MadGraph prompt:
install SM_HeavyN_CKM_AllMasses_LO
exit
```

**Option B: Manual Installation**

Download from FeynRules database or HEPForge and place in `mg5/models/`.

**Verify Installation:**

```bash
ls mg5/models/ | grep -i heavy
# Should show: SM_HeavyN_CKM_AllMasses_LO
```

### 3. Conda Environment

Ensure `mg5env` conda environment exists and is activated:

```bash
# Check available environments
conda env list

# If mg5env doesn't exist, create it
conda create -n mg5env python=3.11

# Activate
conda activate mg5env
```

### 4. Verify Setup

```bash
# Test MadGraph executable
conda run -n mg5env mg5/bin/mg5_aMC --help

# Test Python scripts
conda run -n mg5env python scripts/lhe_to_csv.py --help
```

---

## Pipeline Overview

### Architecture

```
User Script
    │
    ├─> prepare_cards()
    │     ├─ Read templates from cards/
    │     ├─ Fill placeholders (mass, mixing, n_events)
    │     └─ Write to work/hnl_{flavour}_{mass}/
    │
    ├─> run_madgraph()
    │     ├─ Execute MadGraph via conda
    │     ├─ Generate events (LHE format)
    │     └─ Store in work/hnl_{flavour}_{mass}/Events/
    │
    ├─> extract_cross_section()
    │     ├─ Parse MadGraph banner/log
    │     └─ Extract σ_LO ± δσ [pb]
    │
    ├─> convert_lhe_to_csv()
    │     ├─ Parse LHE (XML format)
    │     ├─ Extract HNL and parent W/Z 4-momenta
    │     └─ Write CSV to csv/{flavour}/
    │
    └─> append_to_summary()
          └─ Update csv/summary_HNL_EW_production.csv
```

### Key Components

| Component | File | Purpose |
|-----------|------|---------|
| **Driver** | `scripts/run_hnl_scan.py` | Main orchestration script |
| **LHE Parser** | `scripts/lhe_to_csv.py` | Converts LHE → CSV |
| **Process Cards** | `cards/proc_card_{flavour}.dat` | Define pp → ℓ N processes |
| **Run Card** | `cards/run_card_template.dat` | Collider settings (14 TeV, PDF) |
| **Param Card** | `cards/param_card_template.dat` | HNL mass and mixing |

---

## Usage

### Command-Line Interface

```bash
python scripts/run_hnl_scan.py [OPTIONS]

Options:
  --test              Test mode (15 GeV muon, 1000 events)
  --flavour {electron,muon,tau}
                      Run single flavour only
  --masses M1 M2 ...  Custom mass points in GeV
  --nevents N         Events per mass point (default: 50000)
  -h, --help          Show help message
```

### Examples

**1. Quick Test (Recommended First Run)**

```bash
# 15 GeV muon, 1000 events (~5 min)
python scripts/run_hnl_scan.py --test
```

**2. Single Mass Point**

```bash
# 20 GeV electron, 10k events
python scripts/run_hnl_scan.py --flavour electron --masses 20 --nevents 10000
```

**3. Low-Mass Scan (5-20 GeV)**

```bash
# All flavours, 5-20 GeV in 5 GeV steps
python scripts/run_hnl_scan.py --masses 5 10 15 20
```

**4. Full Production Scan**

```bash
# All masses (5-80 GeV), all flavours, 50k events each
# WARNING: ~10-20 hours total
python scripts/run_hnl_scan.py
```

**5. Resume After Failure**

The pipeline appends to summary CSV and skips existing files. To resume:

```bash
# Just re-run the same command
python scripts/run_hnl_scan.py --flavour muon
```

### Environment Variables

- `MG5_PATH`: Path to MadGraph executable (default: `../mg5/bin/mg5_aMC`)
- `MG5_CONDA_ENV`: Conda environment name (default: `mg5env`)

Example:

```bash
export MG5_PATH=/custom/path/to/mg5_aMC
export MG5_CONDA_ENV=my_mg5_env
python scripts/run_hnl_scan.py
```

---

## Output Files

### Directory Structure

```
production/madgraph/
├── csv/
│   ├── electron/
│   │   └── HNL_mass_15GeV_electron_EW.csv
│   ├── muon/
│   │   └── HNL_mass_15GeV_muon_EW.csv
│   ├── tau/
│   │   └── HNL_mass_15GeV_tau_EW.csv
│   └── summary_HNL_EW_production.csv  ← Cross-sections & metadata
│
├── lhe/
│   └── {flavour}/
│       └── m_{mass}GeV/
│           └── unweighted_events.lhe.gz
│
└── work/
    └── hnl_{flavour}_{mass}GeV/
        ├── proc_card.dat
        ├── run_card.dat
        ├── param_card.dat
        └── madgraph.log
```

### CSV Format (Per-Event Files)

**File**: `csv/{flavour}/HNL_mass_{mass}GeV_{flavour}_EW.csv`

**Header** (EXACT format required by analysis):
```csv
event_id,parent_pdgid,hnl_pdgid,mass_hnl_GeV,weight,parent_E_GeV,parent_px_GeV,parent_py_GeV,parent_pz_GeV,hnl_E_GeV,hnl_px_GeV,hnl_py_GeV,hnl_pz_GeV
```

**Example Row:**
```csv
1,24,9900012,15.0,1.234e-03,87.3,12.4,-8.7,86.5,45.2,8.1,-3.4,44.6
```

**Column Definitions:**

| Column | Type | Description |
|--------|------|-------------|
| `event_id` | int | Event number (1-indexed) |
| `parent_pdgid` | int | Parent boson PDG (24=W+, -24=W-, 23=Z) |
| `hnl_pdgid` | int | HNL PDG (9900012 = N1) |
| `mass_hnl_GeV` | float | HNL mass [GeV] |
| `weight` | float | Event weight from MadGraph |
| `parent_E_GeV` | float | Parent W/Z energy [GeV] |
| `parent_px/y/z_GeV` | float | Parent W/Z momentum [GeV] |
| `hnl_E_GeV` | float | HNL energy [GeV] |
| `hnl_px/y/z_GeV` | float | HNL momentum [GeV] |

### Summary CSV

**File**: `csv/summary_HNL_EW_production.csv`

**Header:**
```csv
mass_hnl_GeV,flavour,xsec_pb,xsec_error_pb,k_factor,n_events_generated,csv_path,timestamp
```

**Example:**
```csv
15.0,muon,4.532e+01,2.1e-01,1.30,50000,csv/muon/HNL_mass_15GeV_muon_EW.csv,2025-11-29 15:30:45
```

**Usage in Analysis:**

```python
import pandas as pd

# Load summary
summary = pd.read_csv('csv/summary_HNL_EW_production.csv')

# Get NLO cross-section for 15 GeV muon
row = summary[(summary['mass_hnl_GeV'] == 15) & (summary['flavour'] == 'muon')]
xsec_nlo = row['xsec_pb'].values[0] * row['k_factor'].values[0]  # pb

# Cross-section scales with |U|²
U2 = 1e-6
xsec_actual = xsec_nlo * U2  # pb
```

---

## Physics Details

### Production Processes

For each flavour ℓ ∈ {e, μ, τ}, we generate:

**Charged Current (W):**
```
pp → W+ → ℓ+ N1
pp → W- → ℓ- N1
```

**Neutral Current (Z):**
```
pp → Z → νℓ N1
pp → Z → ν̄ℓ N1
```

These are combined into a single MadGraph process per flavour.

### Mass Range

**EW Regime:** 5-80 GeV

- **Lower bound (5 GeV)**: Below this, meson production (K/D/B → ℓ N) dominates
- **Upper bound (80 GeV)**: Near W mass; above this, cross-sections become very small

**Default Mass Grid:**
```
5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 25, 28, 30,
32, 35, 38, 40, 45, 50, 55, 60, 65, 70, 75, 80  (GeV)
```

### Mixing Parameters

At generation time, we set **|U_ℓ|² = 1** to get maximum cross-section σ_max.

**Benchmark Configurations:**

| Benchmark | Ve1 | Vmu1 | Vtau1 | Interpretation |
|-----------|-----|------|-------|----------------|
| 100 (electron) | 1 | 0 | 0 | Pure e coupling |
| 010 (muon) | 0 | 1 | 0 | Pure μ coupling |
| 001 (tau) | 0 | 0 | 1 | Pure τ coupling |

**Cross-Section Scaling:**

Since σ ∝ |U_ℓ|², the actual cross-section at mixing |U_ℓ|² is:

```
σ(m_HNL, |U_ℓ|²) = σ_max(m_HNL) × |U_ℓ|²
```

where σ_max is the value in the summary CSV.

**Analysis rescaling:**

```python
# From summary CSV
sigma_max_pb = summary['xsec_pb'] * summary['k_factor']  # NLO

# For |U_mu|² = 1e-6
U_mu_sq = 1e-6
sigma_actual_pb = sigma_max_pb * U_mu_sq
```

### K-Factor

We apply a **K-factor = 1.3** to approximate NLO corrections:

- LO cross-sections from MadGraph
- NLO corrections typically ~20-40% for W/Z production
- Conservative K = 1.3 based on NNLO W/Z studies

**Reference**: CMS/ATLAS W/Z production measurements show K_NLO ~ 1.2-1.4

### Cross-Section Expectations

Typical EW HNL production cross-sections at |U_ℓ|² = 1:

| Mass [GeV] | σ_W [pb] | σ_Z [pb] | σ_total [pb] |
|------------|----------|----------|--------------|
| 5 | ~80 | ~20 | ~100 |
| 15 | ~45 | ~12 | ~57 |
| 30 | ~15 | ~5 | ~20 |
| 60 | ~2 | ~0.8 | ~2.8 |
| 80 | ~0.3 | ~0.1 | ~0.4 |

**Note**: W production dominates over Z by factor ~3-4

---

## Troubleshooting

### 1. "MadGraph executable not found"

**Symptom:**
```
ERROR: MadGraph not found at production/madgraph/mg5/bin/mg5_aMC
```

**Solution:**
```bash
# Check MG5_PATH
echo $MG5_PATH

# Set explicitly
export MG5_PATH=/path/to/mg5/bin/mg5_aMC

# Or install MadGraph (see Installation section)
```

### 2. "Model SM_HeavyN_CKM_AllMasses_LO not found"

**Symptom:**
```
ImportError: Model SM_HeavyN_CKM_AllMasses_LO not found
```

**Solution:**
```bash
# Install model via MadGraph
conda run -n mg5env mg5/bin/mg5_aMC
# In MadGraph prompt:
install SM_HeavyN_CKM_AllMasses_LO
exit
```

### 3. "Conda environment mg5env not found"

**Symptom:**
```
CondaEnvironmentNotFoundError: Could not find conda environment: mg5env
```

**Solution:**
```bash
# Create environment
conda create -n mg5env python=3.11
conda activate mg5env

# Or use different environment
export MG5_CONDA_ENV=base
```

### 4. MadGraph hangs or times out

**Symptom:**
- MadGraph runs for >1 hour
- Process killed with timeout

**Possible causes:**
- Very high multiplicity processes
- Too many events requested
- PDF download issues (first run)

**Solutions:**
```bash
# Reduce events
python scripts/run_hnl_scan.py --test --nevents 100

# Check MadGraph log
cat work/hnl_muon_15GeV/madgraph.log

# Increase timeout in run_hnl_scan.py (line ~400)
```

### 5. "No LHE file found"

**Symptom:**
```
✗ No LHE file found in work/hnl_muon_15GeV
```

**Causes:**
- MadGraph failed silently
- Process generation failed
- Output directory wrong

**Debug:**
```bash
# Check MadGraph log
cat work/hnl_muon_15GeV/madgraph.log

# Look for errors
grep -i error work/hnl_muon_15GeV/madgraph.log

# Check directory structure
ls -R work/hnl_muon_15GeV/
```

### 6. Cross-section extraction fails

**Symptom:**
```
Warning: Could not extract cross-section from work/...
```

**Solution:**
- Check banner file exists: `ls work/hnl_*/Events/*/run_*_banner.txt`
- Manually inspect banner for cross-section
- Update regex pattern in `extract_cross_section()` if MG version differs

### 7. Low event counts in CSV

**Symptom:**
- Requested 10k events
- CSV contains only 100 events

**Causes:**
- MadGraph generated weighted events, not unweighted
- Unweighting efficiency very low
- Process phase space very restrictive

**Solution:**
```bash
# Check run_card.dat: should request unweighted events
grep nevents work/hnl_muon_15GeV/run_card.dat

# Increase event request
python scripts/run_hnl_scan.py --nevents 100000
```

### 8. "Permission denied" when running scripts

**Solution:**
```bash
chmod +x scripts/*.py
```

---

## Comparison with Pythia Pipeline

| Aspect | Pythia | MadGraph (this pipeline) |
|--------|--------|--------------------------|
| **Mass range** | 0.2-13 GeV | 5-80 GeV |
| **Production** | Meson decays (K/D/B) | EW bosons (W/Z) |
| **Method** | Forced decays in Pythia | Matrix element generation |
| **Off-shell** | Fails for m > 13 GeV | Handles automatically |
| **Cross-sections** | Hardcoded from literature | Computed by MadGraph |
| **Output** | production/csv/simulation/ | production/madgraph/csv/ |
| **Use case** | Low-mass regime | High-mass (EW) regime |

**Both pipelines are needed** to cover full 0.2-80 GeV mass range!

---

## References

**MadGraph:**
- MadGraph5_aMC@NLO: [arXiv:1405.0301](https://arxiv.org/abs/1405.0301)
- Download: https://launchpad.net/mg5amcnlo

**HeavyN UFO Model:**
- Based on: [arXiv:1802.02537](https://arxiv.org/abs/1802.02537) (HNL phenomenology)
- Implementation: See `mg5/models/SM_HeavyN_CKM_AllMasses_LO/`

**LLP Detector Studies:**
- MATHUSLA: [arXiv:1811.00927](https://arxiv.org/abs/1811.00927)
- CODEX-b: [arXiv:1911.00481](https://arxiv.org/abs/1911.00481)
- ANUBIS: [arXiv:1909.13022](https://arxiv.org/abs/1909.13022)

**Analysis Pipeline:**
- See `../../analysis_pbc_test/` for downstream analysis
- HNLCalc physics model: [arXiv:2405.07330](https://arxiv.org/abs/2405.07330)

---

## Contact & Contribution

For questions or issues:
- Check MadGraph log files: `work/hnl_*/madgraph.log`
- Verify HeavyN model installation: `ls mg5/models/`
- Test with `--test` mode first
- See main project documentation: `../../CLAUDE.md`

**Pipeline implemented**: 2025-11-29
**MadGraph version**: 3.6.6
**Python**: 3.11+
