# Moneyplot Generation

This directory contains the script for generating the HNL exclusion "moneyplot" - the final visualization showing exclusion islands for all three lepton flavors.

## Overview

The moneyplot shows:
- **X-axis:** HNL mass (m_N) in GeV (log scale)
- **Y-axis:** Mixing parameter |U_ℓ|² (log scale)
- **Shaded regions:** Excluded parameter space (detector has sensitivity)
- **Red line:** "Too long-lived" boundary (lower |U|², HNL doesn't decay in detector)
- **Blue line:** "Too prompt" boundary (higher |U|², HNL decays before reaching detector)

## Files

```
money_plot/
├── README.md                    # This file
└── plot_money_island.py         # Moneyplot generator
```

## Usage

**Simple execution:**
```bash
cd money_plot
/opt/homebrew/Caskroom/miniconda/base/envs/llpatcolliders/bin/python plot_money_island.py
```

**Or from analysis directory:**
```bash
cd analysis_pbc_test
/opt/homebrew/Caskroom/miniconda/base/envs/llpatcolliders/bin/python ../money_plot/plot_money_island.py
```

## Input

**Required file:** `../output/csv/analysis/HNL_U2_limits_summary.csv`

**Format:**
```csv
mass_GeV,flavour,benchmark,eps2_min,eps2_max,peak_events
0.2,electron,100,1.12e-07,0.01,2.17e10
2.6,muon,010,8.70e-09,2.36e-05,4588
...
```

**Columns:**
- `mass_GeV`: HNL mass point
- `flavour`: electron, muon, or tau
- `benchmark`: Coupling pattern (100/010/001)
- `eps2_min`: Lower exclusion boundary (too long-lived)
- `eps2_max`: Upper exclusion boundary (too prompt)
- `peak_events`: Maximum signal events (at optimal |U|²)

## Output

**File:** `../output/images/HNL_moneyplot_island.png`

**Properties:**
- Format: PNG
- Size: ~100-150 KB
- Resolution: 150 DPI (publication quality)
- Dimensions: 18×5 inches (3 panels)

## How It Works

### 1. Load Data
```python
df = pd.read_csv("../output/csv/analysis/HNL_U2_limits_summary.csv")
```

### 2. Process Each Flavor
For electron, muon, tau:
- Filter by flavor
- Remove NaN values (no sensitivity)
- Deduplicate by mass (keep best sensitivity)

### 3. Handle Regime Transitions
- **Low-mass (<5 GeV):** Meson production - plot with circles
- **High-mass (≥5 GeV):** EW production - plot with squares

### 4. Create Exclusion Islands
- **Shaded region:** Between eps2_min and eps2_max
- **Red line:** Lower boundary (cτ too long)
- **Blue line:** Upper boundary (cτ too short)

### 5. Format and Save
- Log-log axes
- Grid lines
- Legend
- Tight layout
- Save to PNG

## Physics Interpretation

### Island Shape

The characteristic "island" shape arises from:

1. **Left boundary (low mass):**
   - Production cross-section drops → less sensitivity
   - OR kinematic threshold (e.g., m_N < m_kaon)

2. **Lower boundary (red line):**
   - Lifetime too long: λ = βγ × cτ₀ >> detector baseline
   - HNL reaches detector but rarely decays
   - Scales as |U|² decreases → cτ increases

3. **Upper boundary (blue line):**
   - Lifetime too short: λ << detector baseline
   - HNL decays before reaching z=22m
   - Scales as |U|² increases → cτ decreases

4. **Right boundary (island closure):**
   - cτ₀(m) becomes inherently too short
   - Even small |U|² yields prompt decays
   - Boundaries meet → island closes

### Lifetime Sweet Spot

The excluded region represents:
```
10 m < βγ × cτ₀ < 1000 m
```

Where:
- Baseline to detector: ~20-100m
- HNL boost factor: βγ = p/m ~ 3-10 (typical)

## Customization

### Adjust Plot Range

Edit lines 65-68:
```python
if flavour == "tau":
    ax.set_xlim([0.5, 50])   # Tau starts at 0.5 GeV
else:
    ax.set_xlim([0.2, 50])   # Electron/muon from 0.2 GeV
```

### Change Colors

Edit lines 36-38 and 45-47:
```python
ax.fill_between(..., alpha=0.3, color='red')   # Excluded region
ax.plot(..., 'r-', ...)                        # Lower boundary
ax.plot(..., 'b-', ...)                        # Upper boundary
```

### Adjust Resolution

Edit line 78:
```python
plt.savefig(..., dpi=150, ...)  # Increase for higher quality
```

### Modify Markers

Low-mass regime (line 37-38):
```python
marker='o', markersize=3   # Circles for meson production
```

High-mass regime (line 46-47):
```python
marker='s', markersize=4   # Squares for EW production
```

## Data Flow

Complete pipeline:

```
1. Pythia Simulation
   └→ output/csv/simulation_new/HNL_*.csv

2. Geometry Preprocessing (cached)
   └→ output/csv/geometry/HNL_*_geom.csv

3. Limits Calculation
   └→ output/csv/analysis/HNL_U2_limits_summary.csv

4. Moneyplot Generation ← THIS SCRIPT
   └→ output/images/HNL_moneyplot_island.png
```

## Dependencies

**Python packages:**
- `pandas`: Data loading
- `matplotlib`: Plotting
- `numpy`: Numerical operations

**Conda environment:**
```bash
conda activate llpatcolliders
```

## Troubleshooting

**Problem:** "File not found" error
```bash
# Check limits file exists
ls -lh ../output/csv/analysis/HNL_U2_limits_summary.csv

# If missing, run analysis first:
cd ../analysis_pbc_test
/opt/homebrew/Caskroom/miniconda/base/envs/llpatcolliders/bin/python -u limits/run_serial.py
```

**Problem:** Empty or incomplete plot
- Check limits file has data: `wc -l HNL_U2_limits_summary.csv`
- Verify no NaN values: `grep -c "nan" HNL_U2_limits_summary.csv`
- Check mass range: `cut -d, -f1 HNL_U2_limits_summary.csv | sort -n | uniq`

**Problem:** Islands don't close properly
- Need finer mass grid in closure region (4-5.5 GeV)
- Run `production/run_island_closure.sh` for additional mass points
- See `production/ISLAND_CLOSURE_README.md`

## Physics Validation

Expected features of a correct moneyplot:

✓ **Island structure:** Closed contours (not open-ended)
✓ **Smooth boundaries:** No large jumps between mass points
✓ **Flavor ordering:** τ < μ < e in sensitivity (typically)
✓ **Width:** Island spans 2-4 decades in |U|²
✓ **Peak sensitivity:** Around m ~ 0.5-2 GeV (D-meson production)

## References

- **Detector geometry:** `analysis_pbc_test/geometry/per_parent_efficiency.py`
- **Limits calculation:** `analysis_pbc_test/limits/run_serial.py`
- **HNL physics:** `analysis_pbc_test/models/hnl_model_hnlcalc.py`

## Example Output

For a complete dataset (137 mass points):
- **Electrons:** Smooth island from 0.2-5.5 GeV
- **Muons:** Smooth island from 0.2-5.5 GeV
- **Taus:** Smooth island from 0.5-4.5 GeV
- **File size:** ~120 KB
- **Generation time:** < 5 seconds

---

**Last Updated:** November 30, 2024
**Status:** Production-ready
**Dependencies:** Requires completed limits calculation
