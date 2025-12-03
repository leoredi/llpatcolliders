"""
run_anubis_quick_test.py

Quick test of ANUBIS U² limit calculation with just 3 mass points.

This is a faster version for testing purposes.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Use the main script's functions
THIS_FILE = Path(__file__).resolve()
sys.path.insert(0, str(THIS_FILE.parent))

from run_anubis_closure import compute_U2_limit_anubis, OUTPUT_DIR, SIM_DIR

import pandas as pd


def main() -> None:
    print("=" * 70)
    print("ANUBIS QUICK TEST - 3 MASS POINTS")
    print("=" * 70)
    print(f"Simulation dir: {SIM_DIR}")
    print("Mass points: 1.0, 2.6, 5.0 GeV\n")

    mass_points = [1.0, 2.6, 5.0]
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

        except Exception as exc:
            print(f"[ERROR] Mass {mass_GeV:.2f} GeV failed: {exc}")
            continue

    # Display results
    if results:
        df = pd.DataFrame(results)
        print(f"\n{'='*70}")
        print("RESULTS:")
        print(f"{'='*70}\n")
        print(df.to_string(index=False))
        print()

        n_excluded = df["eps2_min"].notna().sum()
        print(f"\nSummary: {n_excluded}/{len(df)} mass points have sensitivity")
    else:
        print("\n[WARN] No results.")

    print("\nDone.\n")


if __name__ == "__main__":
    main()
