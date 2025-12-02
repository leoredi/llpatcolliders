# MadGraph HNL Production - Quick Start Guide

## Production Complete ✅

**Status:** All EW production runs successful (96/96 mass×flavor combinations)

---

## What Was Generated

### Event Files
- **96 CSV files** in `csv/` directory
- **32 mass points:** 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 25, 28, 30, 32, 35, 38, 40, 45, 50, 55, 60, 65, 70, 75, 80 GeV
- **3 flavors:** electron, muon, tau
- **50,000 events per file**
- **~4.8 million total events**

### File Format
```csv
event,weight,hnl_id,parent_pdg,pt,eta,phi,p,E,mass,prod_x_mm,prod_y_mm,prod_z_mm,boost_gamma
```

Compatible with existing Pythia-based analysis pipeline.

---

## Cross-Section Summary

| Mass (GeV) | σ(pp→W/Z→ℓN) [pb] | Relative to 5 GeV |
|------------|-------------------|-------------------|
| 5          | 25,600            | 100%              |
| 15         | 24,400            | 95%               |
| 50         | 12,100            | 47%               |
| 80         | 610               | 2.4%              |

**Key observation:** Cross-section drops by ~40× from 5→80 GeV due to phase space suppression near W mass threshold (M_W = 80.4 GeV).

**Lepton universality:** All flavors have identical cross-sections (within <1% MC uncertainty).

---

## File Locations

```
csv/
├── HNL_5p0GeV_electron_EW.csv
├── HNL_5p0GeV_muon_EW.csv
├── HNL_5p0GeV_tau_EW.csv
├── ...
├── HNL_80p0GeV_electron_EW.csv
├── HNL_80p0GeV_muon_EW.csv
├── HNL_80p0GeV_tau_EW.csv
└── summary_HNL_EW_production.csv  ← Metadata + cross-sections
```

---

## Next Steps

### 1. Run Analysis Pipeline
```bash
cd ../../analysis_pbc

# Geometry processing
python geometry/per_parent_efficiency.py

# Limit calculation
conda run -n llpatcolliders python limits/u2_limit_calculator.py

# Generate plots
conda run -n llpatcolliders python ../money_plot/plot_money_island.py
```

### 2. Validation Checks
```bash
# All files should have 50,001 lines (header + 50k events)
wc -l csv/HNL_*_EW.csv | grep -v "50001" | grep -v "total"

# Check parent PDG codes
tail -n +2 csv/HNL_15p0GeV_muon_EW.csv | cut -d',' -f4 | sort -u
# Should show: -24, 23, 24 (W+, Z, W-)
```

---

## Documentation

- **Full report:** `EW_PRODUCTION_SUCCESS.md`
- **Docker guide:** `DOCKER_README.md`
- **Physics methodology:** `../../CLAUDE.md`

---

*Last updated: 2025-11-30*
