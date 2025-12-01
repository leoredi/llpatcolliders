#!/usr/bin/env python
"""
Serial version of u2_limit_calculator (no multiprocessing)
Regenerates HNL_U2_limits_summary.csv from today's simulation data
"""

import sys
import re
from pathlib import Path
import numpy as np
import pandas as pd

# Setup paths
THIS_FILE = Path(__file__).resolve()
ANALYSIS_DIR = THIS_FILE.parent
REPO_ROOT = ANALYSIS_DIR.parents[1]
OUTPUT_DIR = REPO_ROOT / "output" / "csv"
SIM_DIR = OUTPUT_DIR / "simulation_new"
GEOM_CACHE_DIR = OUTPUT_DIR / "geometry"
ANALYSIS_OUT_DIR = OUTPUT_DIR / "analysis"

ANALYSIS_ROOT = ANALYSIS_DIR.parent
if str(ANALYSIS_ROOT) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_ROOT))

from geometry.per_parent_efficiency import build_drainage_gallery_mesh, preprocess_hnl_csv
from models.hnl_model_hnlcalc import HNLModel
from config.production_xsecs import get_parent_sigma_pb
from limits.u2_limit_calculator import expected_signal_events

def couplings_from_eps2(eps2, benchmark):
    if benchmark == "100":
        return eps2, 0.0, 0.0
    elif benchmark == "010":
        return 0.0, eps2, 0.0
    elif benchmark == "001":
        return 0.0, 0.0, eps2
    else:
        raise ValueError(f"Unknown benchmark: {benchmark}")

def scan_single_mass(mass_val, mass_str, flavour, benchmark, lumi_fb, sim_csv, geom_csv):
    """Process one mass point"""
    print(f"\n[{flavour} {mass_val} GeV] Processing...")

    # 1. Load/compute geometry
    if geom_csv.exists():
        geom_df = pd.read_csv(geom_csv)
    else:
        print(f"  Computing geometry (caching to {geom_csv.name})...")
        mesh = build_drainage_gallery_mesh()
        geom_df = preprocess_hnl_csv(sim_csv, mesh)
        geom_df.to_csv(geom_csv, index=False)

    n_hits = geom_df['hits_tube'].sum() if 'hits_tube' in geom_df.columns else 0
    print(f"  Loaded {len(geom_df)} HNLs, {n_hits} hit detector")

    if len(geom_df) == 0 or n_hits == 0:
        print(f"  WARNING: No hits, skipping")
        return None

    # 2. Scan |U|²
    eps2_scan = np.logspace(-12, -2, 100)
    N_scan = []

    for eps2 in eps2_scan:
        N = expected_signal_events(geom_df, mass_val, eps2, benchmark, lumi_fb)
        N_scan.append(N)

    N_scan = np.array(N_scan)

    # 3. Find exclusion range (N >= 3)
    mask_excluded = (N_scan >= 3.0)

    if not mask_excluded.any():
        print(f"  No sensitivity (peak = {N_scan.max():.1f})")
        return {
            "mass_GeV": mass_val,
            "flavour": flavour,
            "benchmark": benchmark,
            "eps2_min": np.nan,
            "eps2_max": np.nan,
            "peak_events": N_scan.max()
        }

    indices_excl = np.where(mask_excluded)[0]
    eps2_min = eps2_scan[indices_excl[0]]
    eps2_max = eps2_scan[indices_excl[-1]]
    peak_events = N_scan.max()

    print(f"  ✓ Excluded: |U|² ∈ [{eps2_min:.2e}, {eps2_max:.2e}], peak = {peak_events:.0f}")

    return {
        "mass_GeV": mass_val,
        "flavour": flavour,
        "benchmark": benchmark,
        "eps2_min": eps2_min,
        "eps2_max": eps2_max,
        "peak_events": peak_events
    }

def run_flavour(flavour, benchmark, lumi_fb):
    """Run scan for one flavour"""
    print(f"\n{'='*60}")
    print(f"FLAVOUR: {flavour.upper()} (Benchmark {benchmark})")
    print(f"{'='*60}")

    # Find simulation files
    pattern = re.compile(rf"HNL_([0-9]+p[0-9]+)GeV_{flavour}_(kaon|charm|beauty|ew)(?:_direct|_fromTau)?\.csv")

    files = []
    for f in SIM_DIR.glob(f"*{flavour}*.csv"):
        if f.stat().st_size < 1000:  # Skip empty files
            print(f"[SKIP] Empty file: {f.name}")
            continue

        match = pattern.search(f.name)
        if match and "_fromTau" not in f.name:  # Skip fromTau to avoid double counting
            mass_str = match.group(1)
            mass_val = float(mass_str.replace('p', '.'))
            files.append((mass_val, mass_str, f))

    files.sort(key=lambda x: x[0])
    print(f"Found {len(files)} mass points")

    # Process each file
    results = []
    for mass_val, mass_str, sim_csv in files:
        geom_csv = GEOM_CACHE_DIR / f"HNL_{mass_str}GeV_{flavour}_geom.csv"

        res = scan_single_mass(mass_val, mass_str, flavour, benchmark, lumi_fb, sim_csv, geom_csv)
        if res:
            results.append(res)

    return pd.DataFrame(results)

if __name__ == "__main__":
    print("="*60)
    print("SERIAL U² LIMIT CALCULATOR")
    print("="*60)

    GEOM_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    ANALYSIS_OUT_DIR.mkdir(parents=True, exist_ok=True)

    L_HL_LHC_FB = 3000.0
    results_out = ANALYSIS_OUT_DIR / "HNL_U2_limits_summary.csv"

    all_results = []

    # Process each flavour
    for flavour, benchmark in [("electron", "100"), ("muon", "010"), ("tau", "001")]:
        df = run_flavour(flavour, benchmark, L_HL_LHC_FB)
        all_results.append(df)

    # Combine and save
    final_df = pd.concat(all_results, ignore_index=True)
    final_df.to_csv(results_out, index=False)

    print(f"\n{'='*60}")
    print(f"COMPLETE!")
    print(f"Saved {len(final_df)} mass points to:")
    print(f"  {results_out}")
    print(f"{'='*60}")
