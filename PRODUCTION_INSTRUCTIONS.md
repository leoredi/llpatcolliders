# Dense Mass Grid Production Instructions

## Summary

Updated mass grid from **36 → 66 points** with uniform **0.2 GeV spacing** to eliminate sensitivity gaps across all flavor transitions (electron, muon, tau).

## What Changed

### Mass Grid Configuration (`config_mass_grid.py`)

**Before:** 36 mass points with gaps up to 0.8 GeV in critical regions
**After:** 66 mass points with maximum 0.2 GeV spacing throughout

**Key improvements:**
- Low mass (0.2-2.0 GeV): 0.05-0.2 GeV steps through K and D transitions
- Mid mass (2.0-5.0 GeV): 0.2 GeV steps through B transition
- **Critical transition (5.0-8.0 GeV): 0.2 GeV steps through beauty→EW transition** ← This fixes the visible gaps!
- High mass (8.0-20.0 GeV): 0.5-1.0 GeV steps

### Production Requirements

**New mass points to produce:** 36 per flavor

**Total new files needed:**
- Electron: 36 points × 4 channels (kaon/charm/beauty/ew) = ~144 files
- Muon: 36 points × 4 channels = ~144 files
- Tau: ~30 points × 3 channels = ~90 files
- **Grand total: ~378 new simulation files**

## Production Steps

### Option 1: Automated Script (Recommended)

```bash
cd production/pythia_production
./run_gap_filling.sh
```

**Details:**
- Produces 36 new mass points × 3 flavors = 108 jobs
- 200,000 events per mass point
- 10 parallel workers
- **Estimated time: ~12-15 hours** (depending on system)

### Option 2: Manual Production

For specific mass points or testing:

```bash
cd production/pythia_production
conda activate llpatcolliders

# Example: produce 5.4 GeV muon with 200k events
./main_hnl_production 5.4 muon 200000
```

## After Production Completes

### 1. Verify output files

```bash
# Should see ~378 new CSV files
ls -lh ../../output/csv/simulation/HNL_*.csv | wc -l

# Check specific new masses (example)
ls ../../output/csv/simulation/HNL_5p*GeV_muon_*.csv
```

### 2. Re-run analysis with new data

```bash
cd ../../analysis_pbc
conda run -n llpatcolliders python limits/run_serial.py --parallel --workers 10
```

**Time: ~30-40 minutes** with 10 workers

### 3. Regenerate plot

```bash
cd ../money_plot
conda run -n llpatcolliders python plot_money_island.py
```

**Expected result:** Smooth exclusion curves with **no visible gaps** in all three panels!

## New Mass Points List

```
Low mass (D-regime):
  0.20, 1.20, 1.30, 1.50, 1.70, 1.90

Mid mass (B-regime):
  2.20, 2.40, 2.80, 3.20, 3.60, 4.00, 4.40, 4.80

Critical transition (B→EW):
  5.20, 5.40, 5.60, 5.80, 6.20, 6.40, 6.60,
  7.00, 7.20, 7.40, 7.60, 7.80, 8.00

High mass (EW-dominated):
  8.50, 9.00, 9.50, 11.00, 12.00, 14.00, 16.00, 18.00, 19.00
```

## Statistics Increase

Also increased from 100k → 200k events per mass point to improve peak statistics in the transition region.

**Before:** Peak ~600 events at 6.8 GeV
**After:** Expected peak ~1200-1500 events (better sensitivity)

## Troubleshooting

### If production fails for high-mass beauty points (>5 GeV)

This is expected - Pythia beauty production fails at high masses. The EW channel will compensate.

### To check progress during production

```bash
# Monitor running jobs
jobs

# Check most recent output
ls -lht ../../output/csv/simulation/ | head -20

# Check log files
tail -f ../../output/logs/simulation/*.log
```

## Expected Outcome

After completing production → analysis → plotting:

✅ **Smooth, continuous exclusion curves** for all three flavors
✅ **No visible gaps** at ~5-7 GeV transition
✅ **Better statistics** throughout (200k events vs 100k)
✅ **Publication-ready figure**

---

**Status:** Configuration updated ✅
**Next:** Run `./run_gap_filling.sh` to start production
