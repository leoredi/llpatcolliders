"""
test_pipeline.py

Minimal smoke tests for the CMS drainage-gallery HNL pipeline.

This file is ONLY for quick end-to-end checks. It must NOT introduce any
new physics models or alternative geometry implementations.

It assumes the following modules already exist and are the SINGLE source of truth:

  - geometry/per_parent_efficiency.py
      • build_drainage_gallery_mesh()
      • preprocess_hnl_csv(csv_file, mesh, origin)

  - models/hnl_model_hnlcalc.py
      • HNLModel(mass_GeV, Ue2, Umu2, Utau2)
          - .ctau0_m           → proper cτ in metres
          - .production_brs()  → {parent_pdg: BR(parent → N + X)}

  - limits/u2_limit_calculator.py
      • expected_signal_events(geom_df, mass_GeV, eps2, benchmark, lumi_fb)
      • run_reach_scan(...)

IMPORTANT CONSTRAINTS FOR ANY LLM EDITING THIS FILE
===================================================

1. Do NOT modify:
     - geometry/per_parent_efficiency.py
     - models/hnl_model_hnlcalc.py
     - config/production_xsecs.py
     - limits/u2_limit_calculator.py

   from this file. You may only IMPORT and USE them.

2. Do NOT reimplement geometry or cross-sections here.
   All physics must go through the modules above.

3. This script is allowed to:
     - run simple test calls,
     - print diagnostic information,
     - exit with an exception if something is clearly broken.

   It is NOT allowed to introduce new physics logic.

Usage
=====

From the repo root:

  cd analysis_pbc
  python test_pipeline.py

The tests will:

  1. Check that HNLModel + HNLCalc work at one mass point.
  2. Build the drainage-gallery mesh and preprocess one Pythia CSV.
  3. Compute N_sig for a single (mass, eps², benchmark) point.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd


# ----------------------------------------------------------------------
# 0. Basic path setup (local relative imports only)
# ----------------------------------------------------------------------

THIS_FILE = Path(__file__).resolve()
TESTS_DIR = THIS_FILE.parent                         # .../analysis_pbc/tests
ANALYSIS_DIR = TESTS_DIR.parent                      # .../analysis_pbc
REPO_ROOT = ANALYSIS_DIR.parent                      # .../llpatcolliders
OUTPUT_DIR = REPO_ROOT / "output" / "csv"
SIM_DIR = OUTPUT_DIR / "simulation_new"
GEOM_DIR = OUTPUT_DIR / "geometry"

# Ensure we can import the local analysis modules
sys.path.insert(0, str(ANALYSIS_DIR))


# ----------------------------------------------------------------------
# 1. Test 1: HNL model wrapper + HNLCalc
# ----------------------------------------------------------------------

def test_model() -> None:
    """
    Test that HNLModel (HNLCalc wrapper) works at a single mass point.

    This is a *sanity check*:
      - ctau0_m is finite and positive
      - production_brs() returns a non-empty dict
    """
    print("=" * 70)
    print("TEST 1: HNL model (HNLCalc wrapper)")
    print("=" * 70)

    from models.hnl_model_hnlcalc import HNLModel

    mass_GeV = 1.0
    # Example: muon-only benchmark (U_mu^2 = 1e-6)
    Ue2, Umu2, Utau2 = 0.0, 1e-6, 0.0

    print(f"Creating HNLModel for m = {mass_GeV} GeV, "
          f"(Ue2, Umu2, Utau2) = ({Ue2:.1e}, {Umu2:.1e}, {Utau2:.1e})")
    model = HNLModel(mass_GeV=mass_GeV, Ue2=Ue2, Umu2=Umu2, Utau2=Utau2)

    ctau0_m = model.ctau0_m
    print(f"  Proper decay length cτ₀ = {ctau0_m:.3e} m")

    brs = model.production_brs()
    print(f"  Number of parent species with non-zero BR: {len(brs)}")

    # Very mild sanity checks
    if not np.isfinite(ctau0_m) or ctau0_m <= 0:
        raise RuntimeError("ctau0_m is not positive and finite")

    if len(brs) == 0:
        raise RuntimeError("production_brs() returned an empty dict")

    # Show a few example parents
    print("\n  Sample production BRs:")
    for i, (pid, br) in enumerate(sorted(brs.items())[:10]):
        print(f"    PDG {pid:6d}: BR = {br:.3e}")
    print("\n✓ TEST 1 passed.\n")


# ----------------------------------------------------------------------
# 2. Test 2: Geometry preprocessing for one CSV
# ----------------------------------------------------------------------

def test_geometry_preprocessing() -> pd.DataFrame:
    """
    Test that geometry/preprocess_hnl_csv produces the expected columns.

    Returns
    -------
    geom_df : pd.DataFrame
        The preprocessed geometry dataframe for reuse in later tests.
    """
    print("=" * 70)
    print("TEST 2: Geometry preprocessing")
    print("=" * 70)

    from geometry.per_parent_efficiency import (
        build_drainage_gallery_mesh,
        preprocess_hnl_csv,
    )

    # Choose one CSV that should exist (adapt flavour/mass if needed)
    csv_path = SIM_DIR / "HNL_1p00GeV_muon_charm.csv"
    if not csv_path.exists():
        raise FileNotFoundError(f"Expected test CSV not found: {csv_path}")

    GEOM_DIR.mkdir(parents=True, exist_ok=True)
    geom_cache = GEOM_DIR / "HNL_1p00GeV_muon_geom.csv"

    print(f"Input CSV:      {csv_path}")
    print(f"Geometry cache: {geom_cache}")

    mesh = build_drainage_gallery_mesh()
    origin = (0.0, 0.0, 0.0)

    if geom_cache.exists():
        print("Loading cached geometry...")
        geom_df = pd.read_csv(geom_cache)
    else:
        print("Computing geometry (first time, may take a while)...")
        geom_df = preprocess_hnl_csv(str(csv_path), mesh, origin=origin)
        geom_df.to_csv(geom_cache, index=False)
        print(f"Saved geometry cache to: {geom_cache}")

    print(f"\nGeometry DataFrame shape: {geom_df.shape}")
    print(f"Columns: {list(geom_df.columns)}")

    required_cols = [
        "parent_id",
        "weight",
        "beta_gamma",
        "hits_tube",
        "entry_distance",
        "path_length",
    ]
    missing = [c for c in required_cols if c not in geom_df.columns]
    if missing:
        raise RuntimeError(
            f"Geometry preprocessing is missing required columns: {missing}"
        )

    n_hits = int(geom_df["hits_tube"].sum())
    frac_hits = float(geom_df["hits_tube"].mean()) if len(geom_df) > 0 else 0.0
    print(f"\nHits in tube: {n_hits} / {len(geom_df)} "
          f"({100.0 * frac_hits:.2f}%)")

    print("\n✓ TEST 2 passed.\n")
    return geom_df


# ----------------------------------------------------------------------
# 3. Test 3: expected_signal_events on one point
# ----------------------------------------------------------------------

def test_expected_signal_events(geom_df: pd.DataFrame) -> None:
    """
    Test expected_signal_events() on a single (mass, eps², benchmark) point.

    This is just a consistency check that:
      - the call works
      - the result is finite (can be very small, that’s fine)
    """
    print("=" * 70)
    print("TEST 3: expected_signal_events()")
    print("=" * 70)

    from limits.u2_limit_calculator import expected_signal_events

    mass_GeV = 1.0
    lumi_fb = 3000.0
    eps2 = 1e-6       # example muon benchmark coupling
    benchmark = "010" # muon-only

    print(f"Computing N_sig for:")
    print(f"  mass_GeV  = {mass_GeV}")
    print(f"  eps2      = {eps2:.1e}")
    print(f"  benchmark = {benchmark} (muon-only)")
    print(f"  lumi_fb   = {lumi_fb:.1f} fb^-1")

    N_sig = expected_signal_events(
        geom_df=geom_df,
        mass_GeV=mass_GeV,
        eps2=eps2,
        benchmark=benchmark,
        lumi_fb=lumi_fb,
    )

    print(f"\nResult: N_sig = {N_sig:.6e} events")

    if not np.isfinite(N_sig):
        raise RuntimeError("expected_signal_events returned non-finite N_sig")

    print("\n✓ TEST 3 passed.\n")


# ----------------------------------------------------------------------
# 4. Optional: small eps² scan (diagnostic only)
# ----------------------------------------------------------------------

def small_eps2_scan(geom_df: pd.DataFrame) -> None:
    """
    Optional diagnostic: scan a few eps² values and print N_sig.
    This is NOT a full reach calculation (run run_reach_scan for that).
    """
    from limits.u2_limit_calculator import expected_signal_events

    print("=" * 70)
    print("OPTIONAL: small eps² scan (diagnostic)")
    print("=" * 70)

    mass_GeV = 1.0
    lumi_fb = 3000.0
    benchmark = "010"

    eps2_values = [1e-8, 1e-7, 1e-6, 1e-5, 1e-4]

    for eps2 in eps2_values:
        N_sig = expected_signal_events(
            geom_df=geom_df,
            mass_GeV=mass_GeV,
            eps2=eps2,
            benchmark=benchmark,
            lumi_fb=lumi_fb,
        )
        print(f"  eps2 = {eps2:.1e} → N_sig = {N_sig:.6e}")

    print("\n(End of diagnostic scan)\n")


# ----------------------------------------------------------------------
# 5. Main entry point
# ----------------------------------------------------------------------

def main() -> None:
    print("\n" + "=" * 70)
    print("HNL DRAINAGE-GALLERY PIPELINE SMOKE TESTS")
    print("=" * 70 + "\n")

    try:
        test_model()
        geom_df = test_geometry_preprocessing()
        test_expected_signal_events(geom_df)

        # Optional quick scan (comment out if you want it silent)
        small_eps2_scan(geom_df)

        print("=" * 70)
        print("ALL SMOKE TESTS COMPLETED ✓")
        print("=" * 70)
        print("\nFor full reach calculation, run:")
        print("  cd analysis_pbc")
        print("  python limits/u2_limit_calculator.py\n")

    except Exception as exc:
        print("\n" + "=" * 70)
        print("SMOKE TEST FAILED ✗")
        print("=" * 70)
        print(f"{type(exc).__name__}: {exc}")
        print("=" * 70)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()