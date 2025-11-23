"""
run_anubis_closure.py

Compute U² limits for ANUBIS-like vertical shaft geometry.

- Uses:
    * HNLModel + HNLCalc for physics
    * production_xsecs.get_parent_sigma_pb
    * Pythia HNL CSV samples in output/csv/simulation
    * Simplified vertical shaft geometry from anubis_geometry

It computes U² exclusion limits for muon-coupled HNLs in the mass range
0.2 - 10 GeV assuming HL-LHC luminosity L = 3000 fb^-1.

This script is intended for OFFLINE use only.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Tuple, Optional

import numpy as np
import pandas as pd

# ----------------------------------------------------------------------
# 0. Path setup
# ----------------------------------------------------------------------

THIS_FILE = Path(__file__).resolve()
# .../analysis_pbc_test/tests/closure_anubis/run_anubis_closure.py
ANALYSIS_DIR = THIS_FILE.parents[2]  # .../analysis_pbc_test
REPO_ROOT = ANALYSIS_DIR.parent  # .../llpatcolliders

if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))

from geometry.per_parent_efficiency import preprocess_hnl_csv
from limits.u2_limit_calculator import expected_signal_events, scan_eps2_for_mass

# Import anubis_geometry from the same package
sys.path.insert(0, str(THIS_FILE.parent))
from anubis_geometry import build_anubis_shaft_mesh

SIM_DIR = REPO_ROOT / "output" / "csv" / "simulation"
OUTPUT_DIR = REPO_ROOT / "output" / "csv" / "analysis"


# ----------------------------------------------------------------------
# 1. Helper: find CSV for mass/flavour
# ----------------------------------------------------------------------

def find_hnl_csv(mass_GeV: float, flavour: str = "muon") -> Path:
    """
    Look for a HNL Pythia CSV with a name like:

      HNL_mass_<mass>_<flavour>_Meson.csv (for m < 5 GeV)
      HNL_mass_<mass>_<flavour>_EW.csv (for m >= 5 GeV)

    where <mass> is formatted like f"{mass_GeV:g}". If not found exactly,
    falls back to scanning the directory for files containing the mass
    string and flavour token.
    """
    mass_str = f"{mass_GeV:g}"

    # Try both Meson and EW regimes
    for regime in ["Meson", "EW"]:
        pattern = f"HNL_mass_{mass_str}_{flavour}_{regime}.csv"
        candidate = SIM_DIR / pattern
        if candidate.exists():
            return candidate

    # Fall back to scanning
    for regime in ["Meson", "EW"]:
        for f in SIM_DIR.glob(f"HNL_mass_*_{flavour}_{regime}.csv"):
            if mass_str in f.name:
                return f

    raise FileNotFoundError(
        f"Could not find HNL CSV for mass {mass_GeV} GeV, flavour {flavour} in {SIM_DIR}"
    )


# ----------------------------------------------------------------------
# 2. Compute U² limits for one mass point
# ----------------------------------------------------------------------

def compute_U2_limit_anubis(
    mass_GeV: float,
    flavour: str = "muon",
    benchmark: str = "010",
    lumi_fb: float = 3000.0,
    N_limit: float = 3.0,
) -> Tuple[Optional[float], Optional[float], float]:
    """
    Compute the U² exclusion limits for ANUBIS-like vertical shaft.

    Returns
    -------
    eps2_min : float or None
        Lower boundary of excluded region
    eps2_max : float or None
        Upper boundary of excluded region
    peak_events : float
        Maximum signal events across the eps2 scan
    """
    csv_path = find_hnl_csv(mass_GeV=mass_GeV, flavour=flavour)
    print(f"\n[ANUBIS] Mass {mass_GeV:.2f} GeV: Using CSV {csv_path.name}")

    mesh = build_anubis_shaft_mesh()
    origin = (0.0, 0.0, 0.0)

    # Project HNLs into the ANUBIS shaft geometry
    geom_df = preprocess_hnl_csv(str(csv_path), mesh=mesh, origin=origin)

    # Scan eps2 to find exclusion boundaries
    eps2_grid, Nsig, eps2_min, eps2_max = scan_eps2_for_mass(
        geom_df=geom_df,
        mass_GeV=mass_GeV,
        benchmark=benchmark,
        lumi_fb=lumi_fb,
        N_limit=N_limit,
    )

    peak_events = float(np.max(Nsig))

    if eps2_min is not None:
        print(f"[ANUBIS] Mass {mass_GeV:.2f} GeV: Excluded U² ∈ [{eps2_min:.2e}, {eps2_max:.2e}], peak = {peak_events:.2e} events")
    else:
        print(f"[ANUBIS] Mass {mass_GeV:.2f} GeV: No sensitivity (peak = {peak_events:.2e} events)")

    return eps2_min, eps2_max, peak_events


# ----------------------------------------------------------------------
# 3. Main: scan low masses for muons
# ----------------------------------------------------------------------

def main() -> None:
    print("=" * 70)
    print("ANUBIS-LIKE HNL U² LIMIT CALCULATION")
    print("=" * 70)
    print(f"Simulation dir: {SIM_DIR}")
    print("Geometry: simplified vertical shaft (anubis_geometry.py)")
    print("Luminosity: 3000 fb^-1 (HL-LHC)")
    print("Flavour: muon (benchmark '010')")
    print("Mass range: 0.2 - 10.0 GeV\n")

    # Low-mass muon grid (ANUBIS stops working after ~10 GeV)
    mass_points = [
        0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
        1.2, 1.4, 1.6, 1.8, 2.0, 2.3, 2.6, 3.0, 3.4, 3.8, 4.2, 4.6, 5.0,
        6.0, 7.0, 8.0, 9.0, 10.0
    ]

    results = []

    for mass_GeV in mass_points:
        try:
            eps2_min, eps2_max, peak_events = compute_U2_limit_anubis(
                mass_GeV=mass_GeV,
                flavour="muon",
                benchmark="010",
                lumi_fb=3000.0,
                N_limit=3.0,
            )

            results.append({
                "mass_GeV": mass_GeV,
                "flavour": "muon",
                "benchmark": "010",
                "eps2_min": eps2_min,
                "eps2_max": eps2_max,
                "peak_events": peak_events,
            })

        except FileNotFoundError as exc:
            print(f"[WARN] {exc}")
            continue
        except Exception as exc:
            print(f"[ERROR] Mass {mass_GeV:.2f} GeV failed: {exc}")
            continue

    # Save results to CSV
    if results:
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        df = pd.DataFrame(results)
        output_path = OUTPUT_DIR / "ANUBIS_U2_limits_muon.csv"
        df.to_csv(output_path, index=False)
        print(f"\n{'='*70}")
        print(f"Results saved to: {output_path}")
        print(f"{'='*70}\n")

        # Print summary
        n_excluded = df["eps2_min"].notna().sum()
        print(f"Summary:")
        print(f"  Total mass points: {len(df)}")
        print(f"  Mass points with sensitivity: {n_excluded}")
        print(f"  Mass points without sensitivity: {len(df) - n_excluded}")

        if n_excluded > 0:
            sensitive_df = df[df["eps2_min"].notna()]
            print(f"\nSensitive mass range:")
            print(f"  {sensitive_df['mass_GeV'].min():.2f} - {sensitive_df['mass_GeV'].max():.2f} GeV")
            print(f"\nStrongest sensitivity:")
            idx_max = sensitive_df["peak_events"].idxmax()
            best = sensitive_df.loc[idx_max]
            print(f"  Mass = {best['mass_GeV']:.2f} GeV")
            print(f"  Peak events = {best['peak_events']:.2e}")
            print(f"  U² range = [{best['eps2_min']:.2e}, {best['eps2_max']:.2e}]")
    else:
        print("\n[WARN] No results to save.")

    print("\nDone.\n")


if __name__ == "__main__":
    main()
