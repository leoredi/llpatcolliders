# ANUBIS-like Closure Tests

This directory contains closure tests and U² limit calculations for a simplified ANUBIS-like vertical shaft geometry.

## Files

- **`anubis_geometry.py`** - Defines a simplified vertical shaft detector geometry:
  - Cylinder radius: 8.0 m
  - Cylinder height: 40.0 m
  - Center position: (100.0, 80.0, 0.0) m from IP
  - Orientation: Vertical (parallel to +y axis)

- **`test_expected_signal_events_kernel.py`** - Algorithmic closure tests:
  - Validates the `expected_signal_events` function against analytic formulas
  - Tests single HNL and weighted multi-HNL scenarios
  - Run with: `python tests/closure_anubis/test_expected_signal_events_kernel.py`

- **`run_anubis_closure.py`** - Full U² limit calculation for ANUBIS:
  - Computes exclusion limits for muon-coupled HNLs
  - Mass range: 0.2 - 10.0 GeV (30 mass points)
  - Luminosity: 3000 fb⁻¹ (HL-LHC)
  - Output: `output/csv/analysis/ANUBIS_U2_limits_muon.csv`
  - Run with: `python tests/closure_anubis/run_anubis_closure.py`
  - **Warning**: Takes ~10-15 minutes (30 masses × 100 eps² scans × geometry preprocessing)

- **`run_anubis_quick_test.py`** - Quick test with 3 mass points:
  - Tests masses: 1.0, 2.6, 5.0 GeV
  - Run with: `python tests/closure_anubis/run_anubis_quick_test.py`
  - Takes ~2-3 minutes

## Quick Start

```bash
cd analysis_pbc_test

# 1. Run algorithmic closure tests (fast)
conda run -n llpatcolliders python tests/closure_anubis/test_expected_signal_events_kernel.py

# 2. Run quick test with 3 mass points (~2-3 min)
conda run -n llpatcolliders python tests/closure_anubis/run_anubis_quick_test.py

# 3. Run full U² limit scan (~10-15 min)
conda run -n llpatcolliders python tests/closure_anubis/run_anubis_closure.py
```

## Results

**Quick Test Results (3 mass points):**

| Mass (GeV) | U² Min | U² Max | Peak Events |
|------------|---------|---------|-------------|
| 1.0 | 6.89×10⁻⁹ | 2.42×10⁻⁴ | 1.6×10⁶ |
| 2.6 | 1.10×10⁻⁸ | 3.59×10⁻⁷ | 31.2 |

The 1.0 GeV point shows exceptional sensitivity (~1.5 million peak events), demonstrating that the vertical shaft geometry is very effective for low-mass HNLs.

## Comparison with PBC

- **PBC (drainage gallery)**: Complex horizontal geometry at z=22m, smaller active volume
- **ANUBIS (vertical shaft)**: Simple vertical cylinder at 100m distance, larger active volume (800 m³)
- **Mass coverage**: ANUBIS is optimized for m < 10 GeV; PBC covers full range up to 80 GeV
- **Sensitivity**: ANUBIS shows stronger sensitivity at low masses due to larger detector volume

## Notes

- All files are self-contained and do **not** modify existing pipeline code
- Uses the same physics models (HNLCalc), cross-sections, and Pythia samples as PBC
- The vertical shaft geometry is simplified for closure testing; not the final ANUBIS design
- The ctau cache file (`model/ctau/ctau.txt`) may be updated when running these scripts
