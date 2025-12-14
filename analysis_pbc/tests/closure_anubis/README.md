# ANUBIS-like Closure Tests

This directory contains fast algorithmic closure tests for the signal-yield kernel
(`limits/expected_signal.py`) using a simplified ANUBIS-like vertical shaft geometry.

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

## Quick Start

```bash
cd analysis_pbc

# Run algorithmic closure tests (fast)
conda run -n llpatcolliders python tests/closure_anubis/test_expected_signal_events_kernel.py
```

## Notes

- All files are self-contained and do **not** modify existing pipeline code
- Uses the same physics models (HNLCalc), cross-sections, and Pythia samples as PBC
- The vertical shaft geometry is simplified for closure testing; not the final ANUBIS design
