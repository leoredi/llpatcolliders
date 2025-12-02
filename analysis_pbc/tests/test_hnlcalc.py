"""
test_pipeline.py

Minimal smoke tests for the CMS drainage-gallery HNL pipeline.

Constraints:
  - Do NOT modify physics modules from here.
  - Do NOT reimplement geometry or cross-sections here.
  - Only import and call:
      * geometry/per_parent_efficiency.py
      * models/hnl_model_hnlcalc.py
      * limits/u2_limit_calculator.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd


# ----------------------------------------------------------------------
# Path setup: make analysis_pbc importable
# ----------------------------------------------------------------------

THIS_FILE = Path(__file__).resolve()
ANALYSIS_DIR = THIS_FILE.parent              # .../analysis_pbc
REPO_ROOT = ANALYSIS_DIR.parent              # .../llpatcolliders
OUTPUT_DIR = REPO_ROOT / "output" / "csv"
SIM_DIR = OUTPUT_DIR / "simulation"
GEOM_DIR = OUTPUT_DIR / "geometry"

if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))


# ----------------------------------------------------------------------
# 1. Test HNLModel + HNLCalc
# ----------------------------------------------------------------------

def test_model() -> None:
    print("=" * 70)
    print("TEST 1: HNL model (HNLCalc wrapper)")
    print("=" * 70)

    from models.hnl_model_hnlcalc import HNLModel

    mass_GeV = 1.0
    # muon-only benchmark
    Ue2, Umu2, Utau2 = 0.0, 1e-6, 0.0

    print(f"Creating HNLModel for m = {mass_GeV} GeV, "
          f"(Ue2, Umu2, Utau2) = ({Ue2:.1e}, {Umu2:.1e}, {Utau2:.1e})")
    model = HNLModel(mass_GeV=mass_GeV, Ue2=Ue2, Umu2=Umu2, Utau2=Utau2)

    ctau0_m = model.ctau0_m
    print(f"  cτ₀ = {ctau0_m:.3e} m")

    brs = model.production_brs()
    print(f"  parents with non-zero BR: {len(brs)}")

    if not np.isfinite(ctau0_m) or ctau0_m <= 0:
        raise RuntimeError("ctau0_m is not positive and finite")

    if len(brs) == 0:
        raise RuntimeError("production_brs() returned an empty dict")

    print("\n  Sample BRs:")
    for pid, br in list(sorted(brs.items()))[:10]:
        print(f"    PDG {pid:6d}: BR = {br:.3e}")

    print("\n✓ TEST 1 passed.\n")


# ----------------------------------------------------------------------
# 2. Geometry preprocessing for one CSV
# ----------------------------------------------------------------------

def test_geometry_preprocessing() -> pd.DataFrame:
    print("=" * 70)
    print("TEST 2: Geometry preprocessing")
    print("=" * 70)

    from geometry.per_parent_efficiency import (
        build_drainage_gallery_mesh,
        preprocess_hnl_csv,
    )

    csv_path = SIM_DIR / "HNL_mass_1.0_muon_Meson.csv"
    if not csv_path.exists():
        raise FileNotFoundError(f"Missing test CSV: {csv_path}")

    GEOM_DIR.mkdir(parents=True, exist_ok=True)
    geom_cache = GEOM_DIR / "HNL_mass_1.0_muon_geom.csv"

    print(f"Input CSV:      {csv_path}")
    print(f"Geometry cache: {geom_cache}")

    mesh = build_drainage_gallery_mesh()
    origin = (0.0, 0.0, 0.0)

    if geom_cache.exists():
        print("Loading cached geometry...")
        geom_df = pd.read_csv(geom_cache)
    else:
        print("Computing geometry (may take a while)...")
        geom_df = preprocess_hnl_csv(str(csv_path), mesh, origin=origin)
        geom_df.to_csv(geom_cache, index=False)
        print(f"Saved geometry cache to {geom_cache}")

    print(f"\nGeometry shape: {geom_df.shape}")
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
        raise RuntimeError(f"Missing required geometry columns: {missing}")

    n_hits = int(geom_df["hits_tube"].sum())
    frac_hits = float(geom_df["hits_tube"].mean()) if len(geom_df) else 0.0
    print(f"\nHits in tube: {n_hits} / {len(geom_df)} "
          f"({100.0*frac_hits:.2f}%)")

    print("\n✓ TEST 2 passed.\n")
    return geom_df


# ----------------------------------------------------------------------
# 3. expected_signal_events on one point
# ----------------------------------------------------------------------

def test_expected_signal_events(geom_df: pd.DataFrame) -> None:
    print("=" * 70)
    print("TEST 3: expected_signal_events()")
    print("=" * 70)

    from limits.u2_limit_calculator import expected_signal_events

    mass_GeV = 1.0
    lumi_fb = 3000.0
    eps2 = 1e-6
    benchmark = "010"  # muon-only

    print("Computing N_sig for:")
    print(f"  mass_GeV  = {mass_GeV}")
    print(f"  eps2      = {eps2:.1e}")
    print(f"  benchmark = {benchmark}")
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
# 4. Optional: small eps² scan
# ----------------------------------------------------------------------

def small_eps2_scan(geom_df: pd.DataFrame) -> None:
    from limits.u2_limit_calculator import expected_signal_events

    print("=" * 70)
    print("OPTIONAL: small eps² scan")
    print("=" * 70)

    mass_GeV = 1.0
    lumi_fb = 3000.0
    benchmark = "010"

    for eps2 in [1e-8, 1e-7, 1e-6, 1e-5, 1e-4]:
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
# 5. Main
# ----------------------------------------------------------------------

def main() -> None:
    print("\n" + "=" * 70)
    print("HNL DRAINAGE-GALLERY PIPELINE SMOKE TESTS")
    print("=" * 70 + "\n")

    try:
        test_model()
        geom_df = test_geometry_preprocessing()
        test_expected_signal_events(geom_df)
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