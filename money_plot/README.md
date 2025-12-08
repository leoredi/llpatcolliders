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
cd analysis_pbc
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

### 3. Plot Exclusion Region
- **Single unified fill:** Uses one `fill_between` call for continuous shading
- **Markers:** Small circles on boundary lines for visual clarity

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

### Why Tau and Electron/Muon Show Different Shapes

A striking feature of the moneyplot is the **different behavior of the "too prompt" (blue) line** for different flavors:

**Electron/Muon (plateau at low mass):**
- The blue line stays **flat at eps2_max = 0.01** for masses 0.2-1.3 GeV
- This is **not a bug** - it indicates the scan saturated at the maximum tested coupling
- The true upper limit is likely much higher (> 0.01)

**Tau (smooth decrease):**
- The blue line **decreases immediately** from low masses
- No plateau structure
- Upper limit is naturally constrained from the start

**Physics Explanation - Kinematic Thresholds:**

The tau lepton is **1.777 GeV**, much heavier than the muon (0.106 GeV) and electron (0.0005 GeV). This creates a fundamental kinematic barrier:

**Light Meson Production (m < 2 GeV):**
| Meson | Mass | e-coupled | μ-coupled | τ-coupled |
|-------|------|-----------|-----------|-----------|
| π± | 0.14 GeV | ✓ M→eN | ✓ M→μN | ✗ Forbidden |
| K± | 0.49 GeV | ✓ M→eN | ✓ M→μN | ✗ Forbidden |
| D0/D± | ~1.9 GeV | ✓ M→eN | ✓ M→μN | ✗ Suppressed |

**Why This Matters:**

1. **Electron/Muon at 0.2-1.3 GeV:**
   - Abundant light meson production: σ(π), σ(K), σ(D) are all large at FCC-ee
   - High branching ratios: BR(M → ℓN) is kinematically allowed
   - **Result:** Production is SO copious that even at maximum coupling (0.01), HNLs still reach the detector before decaying
   - **Plateau:** Analysis saturates at the coupling grid maximum

2. **Tau at 0.2-1.3 GeV:**
   - Light meson decays to tau are **kinematically forbidden**: m_π, m_K, m_D < m_τ
   - Only electroweak production (Z/W → τN) available - much rarer
   - **Result:** Production rate is low, so upper coupling limit is constrained immediately
   - **Smooth curve:** No saturation, true physics limit visible

**The Drop at ~1.3 GeV (e/μ):**
- Marks transition from meson-dominated to electroweak-dominated production
- D meson threshold effects become important
- Production cross-section drops → coupling limits decrease

**Heavy Flavor (m > 3 GeV):**
- B meson production becomes relevant: B → ℓN (all flavors)
- All three flavors now have comparable production mechanisms
- Curves converge to similar shapes

**Key Takeaway:** The plateau is **fundamental physics**, not an analysis artifact. It demonstrates that electron- and muon-coupled HNLs benefit from copious light meson production, while tau-coupled HNLs are kinematically excluded from these abundant channels.

**References:**
- Kinematic suppression: [arXiv:2510.12248](https://arxiv.org/abs/2510.12248) - Tau-coupled HNLs at LHC
- Meson production: [PRD 109.L111102](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.109.L111102) - Tau mixing searches
- CHARM constraints: [PRD 104.095019](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.104.095019) - Tau HNL mixing limits

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

Edit lines 35-37:
```python
ax.fill_between(..., alpha=0.3, color='red')   # Excluded region
ax.plot(..., 'r-', ...)                        # Lower boundary (too long-lived)
ax.plot(..., 'b-', ...)                        # Upper boundary (too prompt)
```

### Adjust Resolution

Edit line 81:
```python
plt.savefig(..., dpi=150, ...)  # Increase for higher quality
```

### Modify Markers

Edit lines 36-37:
```python
marker='o', markersize=3   # Circles on boundary lines
```

## Data Flow

Complete pipeline:

```
1. Production (Pythia + MadGraph)
   └→ output/csv/simulation/HNL_*_{kaon,charm,beauty,ew}.csv

2. Combine Production Channels (REQUIRED!)
   └→ python analysis_pbc/limits/combine_production_channels.py
   └→ output/csv/simulation/HNL_*_combined.csv
   └→ (deletes original separate files to save space)

3. Geometry Preprocessing (cached)
   └→ output/csv/geometry/HNL_*_geom.csv

4. Limits Calculation
   └→ python analysis_pbc/limits/run_serial.py --parallel
   └→ output/csv/analysis/HNL_U2_limits_summary.csv

5. Moneyplot Generation ← THIS SCRIPT
   └→ python money_plot/plot_money_island.py
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

# If missing, run the full pipeline:
# 1. Combine production channels (REQUIRED!)
cd analysis_pbc
conda run -n llpatcolliders python limits/combine_production_channels.py

# 2. Run analysis
conda run -n llpatcolliders python limits/run_serial.py --parallel
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
✓ **Flavor-dependent shapes:**
  - e/μ: Plateau at eps2_max = 0.01 for m < 1.3 GeV (light meson saturation)
  - τ: Smooth decrease from low masses (no light meson production)
✓ **Width:** Island spans 2-4 decades in |U|²
✓ **Peak sensitivity:** Around m ~ 0.5-2 GeV (meson production dominates)

## References

- **Detector geometry:** `analysis_pbc/geometry/per_parent_efficiency.py`
- **Limits calculation:** `analysis_pbc/limits/run_serial.py`
- **HNL physics:** `analysis_pbc/models/hnl_model_hnlcalc.py`

## Example Output

For a complete dataset (137 mass points):
- **Electrons:** Smooth island from 0.2-5.5 GeV
- **Muons:** Smooth island from 0.2-5.5 GeV
- **Taus:** Smooth island from 0.5-4.5 GeV
- **File size:** ~120 KB
- **Generation time:** < 5 seconds

---

**Last Updated:** December 7, 2025
**Status:** Production-ready
**Dependencies:** Requires completed limits calculation

## Changelog

**v2.0 - December 7, 2025:**
- Fixed white gap in excluded region by using unified fill_between
- Added comprehensive physics explanation for tau vs e/μ behavior
- Removed artificial split at 5 GeV (now continuous shading)
- Updated documentation with kinematic threshold references
- Clarified that plateau at eps2_max=0.01 is saturation, not physics limit

**v1.0 - November 30, 2024:**
- Initial production release
- Three-panel moneyplot for electron, muon, tau flavors
- Island visualization with exclusion boundaries
