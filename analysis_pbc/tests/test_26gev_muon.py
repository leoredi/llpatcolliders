"""
Test U² limit calculation for 2.6 GeV HNL with muon coupling.

This test demonstrates the full PBC methodology for a single mass point.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Path setup
THIS_FILE = Path(__file__).resolve()
TESTS_DIR = THIS_FILE.parent
ANALYSIS_DIR = TESTS_DIR.parent
REPO_ROOT = ANALYSIS_DIR.parent
SIM_DIR = REPO_ROOT / "output" / "csv" / "simulation_new"
GEOM_DIR = REPO_ROOT / "output" / "csv" / "geometry"

sys.path.insert(0, str(ANALYSIS_DIR))

from geometry.per_parent_efficiency import build_drainage_gallery_mesh, preprocess_hnl_csv
from limits.u2_limit_calculator import scan_eps2_for_mass


def main():
    MASS_GEV = 2.6
    FLAVOUR = "muon"
    BENCHMARK = "010"  # Pure muon coupling
    LUMI_FB = 3000.0   # HL-LHC

    print("=" * 70)
    print(f"U² Limit Calculation: m_HNL = {MASS_GEV} GeV, {FLAVOUR} coupling")
    print("=" * 70)

    # Load or compute geometry
    # Note: Simulation files use "p" instead of "." in filenames (e.g., 2p60 not 2.6)
    mass_str = str(MASS_GEV).replace(".", "p") + "0"  # 2.6 → 2p60
    csv_path = SIM_DIR / f"HNL_{mass_str}GeV_{FLAVOUR}_beauty.csv"
    geom_cache = GEOM_DIR / f"HNL_{mass_str}GeV_{FLAVOUR}_geom.csv"

    if not csv_path.exists():
        print(f"ERROR: Simulation file not found: {csv_path}")
        return

    GEOM_DIR.mkdir(parents=True, exist_ok=True)

    print(f"\n[1/2] Geometry Processing")
    print(f"  CSV: {csv_path.name}")

    if geom_cache.exists():
        print(f"  Loading cached geometry...")
        geom_df = pd.read_csv(geom_cache)
    else:
        print(f"  Computing geometry (first run, may take time)...")
        mesh = build_drainage_gallery_mesh()
        geom_df = preprocess_hnl_csv(str(csv_path), mesh, origin=(0.0, 0.0, 0.0))
        geom_df.to_csv(geom_cache, index=False)

    n_total = len(geom_df)
    n_hits = geom_df["hits_tube"].sum()
    print(f"  Total HNLs: {n_total}")
    print(f"  Hit detector: {n_hits} ({100*n_hits/n_total:.2f}%)")

    # Parent distribution
    print(f"\n  Parent meson distribution:")
    parent_counts = geom_df["parent_id"].abs().value_counts()
    for pid, count in parent_counts.head(8).items():
        print(f"    PDG {int(pid):>4}: {count:>6} events ({100*count/n_total:>5.1f}%)")

    # Scan U²
    print(f"\n[2/2] Scanning |U_mu|² for 95% CL limit (N_sig = 3 events)")
    print(f"  Luminosity: {LUMI_FB:.0f} fb⁻¹")

    eps2_grid, Nsig, eps2_min, eps2_max = scan_eps2_for_mass(
        geom_df=geom_df,
        mass_GeV=MASS_GEV,
        benchmark=BENCHMARK,
        lumi_fb=LUMI_FB,
        N_limit=3.0,
    )

    # Results
    print("\n" + "=" * 70)
    print("RESULTS")
    print("=" * 70)
    print(f"Mass:      {MASS_GEV} GeV")
    print(f"Flavour:   {FLAVOUR} (benchmark {BENCHMARK})")
    print(f"Lumi:      {LUMI_FB:.0f} fb⁻¹")
    print(f"Peak N_sig: {Nsig.max():.2e} events")

    if eps2_min is not None:
        print(f"\n95% CL Exclusion Range:")
        print(f"  |U_mu|²_min = {eps2_min:.3e}")
        print(f"  |U_mu|²_max = {eps2_max:.3e}")
        print(f"  Island width: {np.log10(eps2_max/eps2_min):.2f} decades")
    else:
        print(f"\nNO SENSITIVITY (peak events < 3)")

    # Detailed scan table
    print("\n" + "-" * 70)
    print("Detailed Scan Points:")
    print(f"{'|U_mu|²':>12}  {'N_sig':>12}  {'Excluded?':>10}")
    print("-" * 70)

    for i in [0, 20, 40, 50, 60, 80, 99]:
        excluded = "YES" if Nsig[i] >= 3.0 else "NO"
        print(f"{eps2_grid[i]:>12.2e}  {Nsig[i]:>12.2e}  {excluded:>10}")

    print("=" * 70)


if __name__ == "__main__":
    main()
