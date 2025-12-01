# HNL Analysis Pipeline

This directory contains the complete analysis pipeline for HNL exclusion limit calculation.

## Overview

The analysis chain:
1. **Geometry preprocessing** → Ray-trace HNL trajectories to detector (cached)
2. **HNL physics model** → Calculate lifetimes and branching ratios (HNLCalc)
3. **Signal calculation** → Compute expected events for each |U|² value
4. **Limit extraction** → Find exclusion boundaries (eps2_min, eps2_max)

## Quick Start

```bash
# Run complete analysis
/opt/homebrew/Caskroom/miniconda/base/envs/llpatcolliders/bin/python -u limits/run_serial.py

# Output: ../output/csv/analysis/HNL_U2_limits_summary.csv
```

**Runtime:** ~2-3 hours for 137 mass points

## Key Feature: Geometry Caching

Geometry preprocessing is **automatically cached** for efficiency:

```python
if geom_csv.exists():
    geom_df = pd.read_csv(geom_csv)          # Load from cache (<1 sec)
else:
    mesh = build_drainage_gallery_mesh()
    geom_df = preprocess_hnl_csv(sim_csv, mesh)
    geom_df.to_csv(geom_csv, index=False)   # Save to cache (~30 sec)
```

**Cache location:** `../output/csv/geometry/`
**Speedup:** 30x faster on subsequent runs

---

For complete documentation, see `../CLAUDE.md`
