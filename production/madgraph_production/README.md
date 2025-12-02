# MadGraph HNL Production Pipeline (Docker)

**Electroweak Heavy Neutral Lepton (HNL) production at LHC 14 TeV**

This pipeline generates HNL events via W/Z boson production using MadGraph5_aMC@NLO with the HeavyN UFO model, running in a Docker container for reproducibility and portability.

## ✅ Production Status

**Complete:** 96 CSV files generated covering 5-80 GeV (32 mass points × 3 flavors)
- Output location: `../../output/csv/simulation_new/HNL_*_EW.csv`
- Total events: ~4.8 million (50,000 per mass point)
- See: `EW_PRODUCTION_SUCCESS.md` for full report
- See: `QUICKSTART.md` for user guide

---

## Quick Start

### 1. Build the Docker Image (Once)

```bash
cd /path/to/llpatcolliders/production/madgraph
docker build -t mg5-hnl .
```

This creates a container with:
- MadGraph5_aMC@NLO 3.6.6 at `/opt/MG5_aMC_v3_6_6`
- Python 3 + required dependencies (numpy, pandas, pylhe)
- Pythia8 and LHAPDF6 (auto-installed by MadGraph)

### 2. Run the Container

```bash
# From llpatcolliders repository root:
cd /path/to/llpatcolliders

# Launch container with repo mounted at /work
docker run --rm -it -v "$(pwd)":/work mg5-hnl bash
```

Inside the container, navigate to the MadGraph pipeline:

```bash
cd /work/production/madgraph
```

### 3. Test with Single Point

```bash
# Inside container:
python3 scripts/run_hnl_scan.py --test
```

This generates 1000 events for a 15 GeV muon-coupled HNL (~2-5 minutes).

**Expected output:**
```
✓ Process directory created: /work/production/madgraph/work/hnl_muon_15.0GeV
✓ Events generated: .../Events/run_01/unweighted_events.lhe.gz
✓ Cross-section: ~45 pb
✓ Converted 1000 events to CSV
```

---

## Why MadGraph + Docker?

### The Pythia Problem

The existing Pythia pipeline (`production/main_hnl_production.cc`) works well for **meson production** (m < 5 GeV) but **fails** for electroweak production (m ≥ 5 GeV) because:

1. **Off-shell W/Z**: Many bosons produced below W/Z mass threshold
2. **Kinematic blocking**: When m(W*) < m_HNL + m_ℓ, decay is forbidden
3. **No fallback**: Setting BR(W → all other modes) = 0 leaves no open channels
4. **Pythia aborts**: "Error: No allowed decay modes"

### The MadGraph Solution

MadGraph computes **pp → ℓ N** via W/Z exchange at the matrix element level:

✅ **Handles off-shell**: W*/Z* propagators with proper kinematics
✅ **Phase space**: Correct sampling including threshold regions
✅ **No decay hacks**: Process is pp → ℓ N (W is internal)
✅ **Works 5-80 GeV**: Seamless across full EW regime

### Why Docker?

✅ **Reproducibility**: Same environment on macOS/Linux/clusters
✅ **No conda conflicts**: Isolated from host Python/Fortran
✅ **Portable**: Share exact setup via Dockerfile
✅ **Clean**: No permanent installation on host system

---

## Usage

### Command-Line Interface

```bash
# Inside Docker container at /work/production/madgraph:

python3 scripts/run_hnl_scan.py [OPTIONS]

Options:
  --test                Test mode (15 GeV muon, 1000 events)
  --flavour {electron,muon,tau}
                        Run single flavour only
  --masses M1 M2 ...    Custom mass points in GeV
  --nevents N           Events per mass point (default: 50000)
  -h, --help            Show help message
```

### Examples

**1. Quick Test (5 minutes)**

```bash
python3 scripts/run_hnl_scan.py --test
```

**2. Single Mass Point**

```bash
# 20 GeV electron, 10k events (~10 min)
python3 scripts/run_hnl_scan.py --flavour electron --masses 20 --nevents 10000
```

**3. Multiple Mass Points**

```bash
# All flavours, 5-20 GeV in 5 GeV steps (~2 hours)
python3 scripts/run_hnl_scan.py --masses 5 10 15 20
```

**4. Single Flavour Scan**

```bash
# Muon only, default mass grid, 50k events each (~8 hours)
python3 scripts/run_hnl_scan.py --flavour muon
```

**5. Full Production Run**

```bash
# All masses (5-80 GeV), all flavours, 50k events each
# WARNING: ~20-40 hours total
python3 scripts/run_hnl_scan.py
```

---

## Output Files

### Directory Structure

```
production/madgraph/
├── csv/
│   ├── electron/
│   │   └── HNL_mass_15.0GeV_electron_EW.csv
│   ├── muon/
│   │   └── HNL_mass_15.0GeV_muon_EW.csv
│   ├── tau/
│   │   └── HNL_mass_15.0GeV_tau_EW.csv
│   └── summary_HNL_EW_production.csv  ← Cross-sections & metadata
│
├── work/
│   └── hnl_{flavour}_{mass}GeV/
│       ├── Cards/
│       │   ├── param_card.dat
│       │   └── run_card.dat
│       ├── Events/
│       │   └── run_01/
│       │       └── unweighted_events.lhe.gz
│       └── SubProcesses/
```

### CSV Format (Per-Event Files)

**File**: `csv/{flavour}/HNL_mass_{mass}GeV_{flavour}_EW.csv`

**Header**:
```csv
event_id,parent_pdgid,hnl_pdgid,mass_hnl_GeV,weight,parent_E_GeV,parent_px_GeV,parent_py_GeV,parent_pz_GeV,hnl_E_GeV,hnl_px_GeV,hnl_py_GeV,hnl_pz_GeV
```

**Example Row**:
```csv
1,24,9900012,15.0,1.234e-03,87.3,12.4,-8.7,86.5,45.2,8.1,-3.4,44.6
```

**Columns**:
- `event_id`: Event number (1-indexed)
- `parent_pdgid`: Parent boson (24=W⁺, -24=W⁻, 23=Z⁰)
- `hnl_pdgid`: HNL PDG code (9900012 = n1)
- `mass_hnl_GeV`: HNL mass [GeV]
- `weight`: Event weight from MadGraph
- `parent_E/px/py/pz_GeV`: Parent W/Z 4-momentum [GeV]
- `hnl_E/px/py/pz_GeV`: HNL 4-momentum [GeV]

**Note**: Parent W/Z is **reconstructed** from ℓ + N 4-vectors (W/Z is internal propagator in MadGraph).

### Summary CSV

**File**: `csv/summary_HNL_EW_production.csv`

**Header**:
```csv
mass_hnl_GeV,flavour,xsec_pb,xsec_error_pb,k_factor,n_events_generated,csv_path,timestamp
```

**Example**:
```csv
15.0,muon,4.532e+01,2.1e-01,1.30,50000,csv/muon/HNL_mass_15.0GeV_muon_EW.csv,2025-11-29 15:30:45
```

**Usage in Analysis**:

```python
import pandas as pd

# Load summary
summary = pd.read_csv('csv/summary_HNL_EW_production.csv')

# Get NLO cross-section for 15 GeV muon at |Uμ|² = 1
row = summary[(summary['mass_hnl_GeV'] == 15) & (summary['flavour'] == 'muon')]
xsec_max = row['xsec_pb'].values[0] * row['k_factor'].values[0]  # pb

# Scale to actual mixing |Uμ|²
U_mu_sq = 1e-6
xsec_actual = xsec_max * U_mu_sq  # pb
```

---

## Physics Details

### Production Processes

For each flavour ℓ ∈ {e, μ, τ}, we generate:

**Charged Current (W)**:
```
pp → W⁺ → ℓ⁺ N₁
pp → W⁻ → ℓ⁻ N₁
```

**Neutral Current (Z)**:
```
pp → Z⁰ → νℓ N₁
pp → Z⁰ → ν̄ℓ N₁
```

**MadGraph implementation**: These are treated as **2→2 hard processes** at the matrix element level (W/Z is internal propagator, not produced on-shell).

### Mass Range & Grid

**EW Regime**: 5-80 GeV

- **Lower bound (5 GeV)**: Below this, meson production (K/D/B → ℓ N) dominates (use Pythia pipeline)
- **Upper bound (80 GeV)**: Near W mass; above this, cross-sections become tiny

**Default Mass Grid** (32 points):
```
5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
22, 25, 28, 30, 32, 35, 38, 40, 45, 50, 55, 60, 65, 70, 75, 80
```

### Mixing Parameters

At generation time, we set **|U_ℓ|² = 1** to get maximum cross-section σ_max.

**Benchmark Configurations**:

| Benchmark | ve1 | vmu1 | vtau1 | Interpretation |
|-----------|-----|------|-------|----------------|
| 100 (electron) | 1 | 0 | 0 | Pure e coupling |
| 010 (muon) | 0 | 1 | 0 | Pure μ coupling |
| 001 (tau) | 0 | 0 | 1 | Pure τ coupling |

**Cross-Section Scaling**:

Since σ ∝ |U_ℓ|², the actual cross-section at mixing |U_ℓ|² is:

```
σ(m_HNL, |U_ℓ|²) = σ_max(m_HNL) × |U_ℓ|²
```

where σ_max is stored in `summary_HNL_EW_production.csv`.

### K-Factor

We apply **K-factor = 1.3** to approximate NLO corrections:

- MadGraph provides LO cross-sections
- NLO corrections for W/Z production typically ~20-40%
- Conservative K = 1.3 based on CMS/ATLAS measurements

**Reference**: NNLO W/Z studies show K_NLO ~ 1.2-1.4 at LHC 14 TeV.

### Cross-Section Expectations

Typical EW HNL production cross-sections at |U_ℓ|² = 1:

| Mass [GeV] | σ_total [pb] | Dominant Channel |
|------------|--------------|------------------|
| 5          | ~100         | W⁺/W⁻            |
| 15         | ~57          | W⁺/W⁻            |
| 30         | ~20          | W⁺/W⁻            |
| 60         | ~2.8         | W⁺/W⁻            |
| 80         | ~0.4         | W⁺/W⁻            |

**Note**: W production dominates over Z by factor ~3-4.

---

## Model & PDG Codes

**UFO Model**: `SM_HeavyN_CKM_AllMasses_LO`

The Docker image includes this model with 3 heavy neutrinos:

| Particle | PDG Code | Mass | Usage |
|----------|----------|------|-------|
| n1 | 9900012 | Scan parameter | **Active HNL** (used in processes) |
| n2 | 9900014 | 1000 GeV | Decoupled (heavy) |
| n3 | 9900016 | 1000 GeV | Decoupled (heavy) |

**Why 3 generations?** The UFO model supports arbitrary mixing patterns. We set n2/n3 very heavy to decouple them, focusing on n1.

**Mixing parameters** (in `param_card.dat`):
- `Block numixing`: 3×3 matrix for (n1, n2, n3) × (e, μ, τ)
- We set only one non-zero entry per flavour benchmark

---

## Pipeline Architecture

### Three-Step Workflow

```
1. Generate Process
   ├─ Call MadGraph: import model + generate p p > ℓ n1
   ├─ Create process directory: work/hnl_{flavour}_{mass}GeV/
   └─ Output: SubProcesses/, bin/, Cards/ structure

2. Write Cards
   ├─ Copy param_card_template.dat → Cards/param_card.dat
   ├─ Replace placeholders: MASS_N1_PLACEHOLDER → {mass}
   ├─ Replace mixing: VE1_PLACEHOLDER → 0 or 1
   ├─ Copy run_card_template.dat → Cards/run_card.dat
   └─ Replace N_EVENTS_PLACEHOLDER → {nevents}

3. Run Event Generation
   ├─ Execute: bin/generate_events -f --laststep=parton
   ├─ MadGraph integrates matrix element, generates unweighted events
   └─ Output: Events/run_01/unweighted_events.lhe.gz

4. Convert LHE → CSV
   ├─ Parse LHE (XML format)
   ├─ Extract HNL (PDG 9900012) and charged lepton
   ├─ Reconstruct parent W/Z from ℓ + N 4-vectors
   └─ Write CSV to csv/{flavour}/

5. Extract Cross-Section & Update Summary
   ├─ Parse Events/run_01/run_01_tag_1_banner.txt
   ├─ Extract σ_LO ± δσ [pb]
   └─ Append to csv/summary_HNL_EW_production.csv
```

---

## Troubleshooting

### 1. Docker Image Not Found

**Symptom**:
```
docker: image mg5-hnl:latest not found
```

**Solution**:
```bash
cd /path/to/llpatcolliders/production/madgraph
docker build -t mg5-hnl .
```

### 2. MadGraph Hangs During Event Generation

**Symptom**: Process runs for >1 hour without output

**Possible Causes**:
- Very high mass (m > 70 GeV) → low cross-section → slow integration
- First run → PDF download (LHAPDF)

**Solutions**:
```bash
# Reduce events for testing
python3 scripts/run_hnl_scan.py --masses 15 --nevents 100

# Check MadGraph log
cat work/hnl_muon_15.0GeV/Events/run_01/*.log

# If stuck downloading PDFs, wait ~10 min (first run only)
```

### 3. "No LHE File Found"

**Symptom**:
```
✗ No LHE file found in work/hnl_muon_15.0GeV
```

**Debug**:
```bash
# Check MadGraph output
ls -R work/hnl_muon_15.0GeV/Events/

# Look for errors in logs
grep -i error work/hnl_muon_15.0GeV/Events/run_01/*.log

# Check if matrix elements compiled
ls work/hnl_muon_15.0GeV/SubProcesses/P*/
```

### 4. Low Event Counts in CSV

**Symptom**: Requested 10k events, CSV contains only 100

**Cause**: MadGraph's `run_card.dat` has incorrect `nevents` setting

**Solution**:
```bash
# Check run card
cat work/hnl_muon_15.0GeV/Cards/run_card.dat | grep nevents

# Should show:
#   10000 = nevents
```

If template substitution failed, edit `cards/run_card_template.dat`.

### 5. Cross-Section Extraction Fails

**Symptom**:
```
Warning: Could not extract cross-section from work/...
```

**Solution**:
```bash
# Manually check banner
cat work/hnl_muon_15.0GeV/Events/run_01/run_01_tag_1_banner.txt | grep -A5 "Integrated weight"

# Cross-section should be printed there
```

---

## Integration with Analysis Pipeline

This MadGraph pipeline is **complementary** to the Pythia meson production pipeline:

| Aspect | Pythia (`production/`) | MadGraph (this pipeline) |
|--------|------------------------|--------------------------|
| **Mass range** | 0.2-5 GeV | 5-80 GeV |
| **Production** | Meson decays (K/D/B) | EW bosons (W/Z) |
| **Method** | Forced decays in Pythia | Matrix element generation |
| **Output** | `production/csv/simulation/` | `production/madgraph/csv/` |
| **Analysis** | `analysis_pbc_test/` | Same (compatible CSV format) |

**Both pipelines are needed** to cover the full 0.2-80 GeV mass range!

The analysis code (`analysis_pbc_test/limits/u2_limit_calculator.py`) can consume CSV files from **both** sources (they use the same format).

---

## Development & Customization

### Modifying the UFO Model

If you need a different HNL model:

1. **Install model** inside container:
   ```bash
   # Edit Dockerfile to add:
   RUN echo "install YOUR_MODEL_NAME" > /tmp/model.mg5 && \
       $MG5_DIR/bin/mg5_aMC /tmp/model.mg5
   ```

2. **Update process cards**:
   Edit `cards/proc_card_*.dat` to use new model name.

3. **Rebuild image**:
   ```bash
   docker build -t mg5-hnl .
   ```

### Changing Mass Grid

Edit `scripts/run_hnl_scan.py` line 47-50:

```python
MASS_GRID_FULL = [
    # Your custom mass points
    5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80
]
```

### Debugging Inside Container

```bash
# Run container with interactive shell
docker run --rm -it -v "$(pwd)":/work mg5-hnl bash

# Test MadGraph manually
cd /work/production/madgraph
/opt/MG5_aMC_v3_6_6/bin/mg5_aMC

# In MG prompt:
import model SM_HeavyN_CKM_AllMasses_LO
display particles
# Should show n1, n2, n3
```

---

## References

**MadGraph**:
- MadGraph5_aMC@NLO: [arXiv:1405.0301](https://arxiv.org/abs/1405.0301)
- Download: https://launchpad.net/mg5amcnlo

**HeavyN UFO Model**:
- Based on: [arXiv:1802.02537](https://arxiv.org/abs/1802.02537) (HNL phenomenology)
- Model database: https://feynrules.irmp.ucl.ac.be/

**LLP Detector Studies** (MadGraph usage examples):
- MATHUSLA: [arXiv:1811.00927](https://arxiv.org/abs/1811.00927)
- CODEX-b: [arXiv:1911.00481](https://arxiv.org/abs/1911.00481)
- ANUBIS: [arXiv:1909.13022](https://arxiv.org/abs/1909.13022)

**Analysis Pipeline**:
- See `../../analysis_pbc_test/` for downstream analysis
- HNLCalc physics model: [arXiv:2405.07330](https://arxiv.org/abs/2405.07330)

---

## Status

**Pipeline Status**: ✅ Fully functional
**Last Updated**: 2025-11-29
**MadGraph Version**: 3.6.6
**Docker Base**: Ubuntu 22.04
**Python**: 3.10+

For questions or issues, see main project documentation: `../../CLAUDE.md`
