#!/usr/bin/env python
"""
U² limit calculator driver.

Usage:
    python run.py                       # Single-threaded (Majorana, default)
    python run.py --parallel            # Parallel processing (uses all cores)
    python run.py --parallel --workers 8  # Use 8 cores
    python run.py --dirac               # Dirac HNL interpretation (×2 yield)
    python run.py --separation-mm 1.0   # Minimum charged-track separation (mm)
"""

import sys
import re
import os
import io
import tempfile
import contextlib
from pathlib import Path
import numpy as np
import pandas as pd
import argparse
from concurrent.futures import ProcessPoolExecutor
import multiprocessing
from tqdm import tqdm

# Setup paths
THIS_FILE = Path(__file__).resolve()
ANALYSIS_DIR = THIS_FILE.parent
REPO_ROOT = ANALYSIS_DIR.parents[1]
OUTPUT_DIR = REPO_ROOT / "output" / "csv"

# Current simulation output directory (Pythia + MadGraph)
SIM_DIR = OUTPUT_DIR / "simulation"
GEOM_CACHE_DIR = OUTPUT_DIR / "geometry"
ANALYSIS_OUT_DIR = OUTPUT_DIR / "analysis"

ANALYSIS_ROOT = ANALYSIS_DIR.parent
if str(ANALYSIS_ROOT) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_ROOT))

from geometry.per_parent_efficiency import build_drainage_gallery_mesh, preprocess_hnl_csv
from limits.expected_signal import expected_signal_events

def scan_single_mass(
    mass_val,
    mass_str,
    flavour,
    benchmark,
    lumi_fb,
    sim_files,
    dirac=False,
    separation_m=0.001,
    decay_seed=12345,
    quiet=False,
):
    """
    Process one mass point, combining all available production regimes (kaon/charm/beauty/ew).
    sim_files: list of tuples (sim_csv_path, regime)
    dirac: if True, multiply yield by 2 for Dirac HNL interpretation
    quiet: if True, suppress verbose output (for parallel mode)
    """
    # In quiet mode, redirect stdout to suppress all output from this worker
    stdout_ctx = contextlib.redirect_stdout(io.StringIO()) if quiet else contextlib.nullcontext()

    with stdout_ctx:
        return _scan_single_mass_impl(
            mass_val, mass_str, flavour, benchmark, lumi_fb, sim_files,
            dirac, separation_m, decay_seed, quiet
        )


def _scan_single_mass_impl(
    mass_val,
    mass_str,
    flavour,
    benchmark,
    lumi_fb,
    sim_files,
    dirac,
    separation_m,
    decay_seed,
    quiet,
):
    """Internal implementation of scan_single_mass."""
    if not quiet:
        print(f"\n[{flavour} {mass_val} GeV] Processing ({len(sim_files)} production file(s))...")

    # 1. Load/compute geometry for each production file, then concatenate
    geom_dfs = []
    mesh = None  # Build only if needed

    for sim_csv, regime in sim_files:
        # Regime-specific geometry cache to avoid collisions
        geom_cache_name = f"{sim_csv.stem}_geom.csv"
        geom_csv = GEOM_CACHE_DIR / geom_cache_name

        # Backward compatibility: fall back to old cache naming if present
        legacy_geom_csv = GEOM_CACHE_DIR / f"HNL_{mass_str}GeV_{flavour}_geom.csv"
        if not geom_csv.exists() and legacy_geom_csv.exists():
            geom_csv = legacy_geom_csv

        # Delete stale geometry cache (source file was updated)
        if geom_csv.exists() and sim_csv.stat().st_mtime > geom_csv.stat().st_mtime:
            if not quiet:
                print(f"  Stale cache detected, regenerating: {geom_csv.name}")
            geom_csv.unlink()

        if geom_csv.exists():
            geom_df = pd.read_csv(geom_csv)
        else:
            if mesh is None:
                mesh = build_drainage_gallery_mesh()
            if not quiet:
                print(f"  Computing geometry for {sim_csv.name} (caching to {geom_csv.name})...")
            geom_df = preprocess_hnl_csv(sim_csv, mesh)

            # Atomic write to prevent race condition in parallel execution
            with tempfile.NamedTemporaryFile(mode='w', dir=geom_csv.parent,
                                              suffix='.tmp', delete=False) as tmp:
                geom_df.to_csv(tmp.name, index=False)
            os.replace(tmp.name, geom_csv)

        n_hits = geom_df['hits_tube'].sum() if 'hits_tube' in geom_df.columns else 0
        if not quiet:
            print(f"  Loaded {len(geom_df)} HNLs from {regime}, {n_hits} hit detector")
        geom_dfs.append(geom_df)

    if len(geom_dfs) == 0:
        if not quiet:
            print(f"  WARNING: No geometry loaded, skipping")
        return None

    geom_df = pd.concat(geom_dfs, ignore_index=True)
    n_hits_total = geom_df['hits_tube'].sum() if 'hits_tube' in geom_df.columns else 0
    if not quiet:
        print(f"  Total combined: {len(geom_df)} HNLs, {n_hits_total} hit detector")

    if len(geom_df) == 0 or n_hits_total == 0:
        if not quiet:
            print(f"  WARNING: No hits, skipping")
        return None

    # 2. Scan |U|²
    eps2_scan = np.logspace(-12, -2, 100)
    N_scan = []

    from decay.decay_detector import DecaySelection, build_decay_cache

    decay_cache = build_decay_cache(
        geom_df,
        mass_val,
        flavour,
        DecaySelection(separation_m=separation_m, seed=decay_seed),
        verbose=not quiet,
    )

    if not quiet:
        print("  Precomputing decay cache for separation scan...")
    for eps2 in eps2_scan:
        N = expected_signal_events(
            geom_df,
            mass_val,
            eps2,
            benchmark,
            lumi_fb,
            dirac=dirac,
            separation_m=separation_m,
            decay_seed=decay_seed,
            decay_cache=decay_cache,
        )
        N_scan.append(N)

    N_scan = np.array(N_scan)

    # 3. Find exclusion range (N >= 2.996)
    mask_excluded = (N_scan >= 2.996)

    if not mask_excluded.any():
        if not quiet:
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

    if not quiet:
        print(f"  ✓ Excluded: |U|² ∈ [{eps2_min:.2e}, {eps2_max:.2e}], peak = {peak_events:.0f}")

    return {
        "mass_GeV": mass_val,
        "flavour": flavour,
        "benchmark": benchmark,
        "eps2_min": eps2_min,
        "eps2_max": eps2_max,
        "peak_events": peak_events
    }

def run_flavour(
    flavour,
    benchmark,
    lumi_fb,
    use_parallel=False,
    n_workers=None,
    dirac=False,
    separation_m=0.001,
    decay_seed=12345,
):
    """Run scan for one flavour"""
    print(f"\n{'='*60}")
    print(f"FLAVOUR: {flavour.upper()} (Benchmark {benchmark})")
    print(f"{'='*60}")

    # Find simulation files
    # Accept both 1 and 2 decimal formats: 5p0 or 5p00 (transitioning to 2-decimal standard)
    pattern = re.compile(
        rf"^HNL_([0-9]+p[0-9]{{1,2}})GeV_{flavour}_"
        r"((?:kaon|charm|beauty|ew|combined)(?:_ff)?)"
        r"(?:_(direct|fromTau))?"
        r"\.csv$"
    )

    files = []
    empty_files = []  # Track empty files, report only if they would be used
    for f in SIM_DIR.glob(f"*{flavour}*.csv"):
        match = pattern.search(f.name)
        if not match:
            continue

        if f.stat().st_size < 1000:  # Track empty files for later reporting
            empty_files.append(f)
            continue

        mass_str = match.group(1)
        mass_val = float(mass_str.replace("p", "."))
        regime_token = match.group(2)          # e.g. charm, charm_ff, combined
        mode = match.group(3)                 # direct/fromTau (tau only) or None

        is_ff = regime_token.endswith("_ff")
        base_regime = regime_token.replace("_ff", "")

        files.append((mass_val, mass_str, base_regime, mode, is_ff, f))

    # Group by mass and select files to avoid double counting:
    # - If *_combined.csv exists for a mass, use ONLY that file.
    # - Otherwise, include all regimes; for tau, include both direct and fromTau modes.
    # - Prefer *_ff files over base files for the same (regime, mode).
    files_by_mass = {}
    for mass_val, mass_str, base_regime, mode, is_ff, path in files:
        key = (mass_val, mass_str)
        files_by_mass.setdefault(key, []).append((base_regime, mode, is_ff, path))

    def _label(base_regime: str, mode: str | None, is_ff: bool) -> str:
        label = base_regime
        if is_ff:
            label += "_ff"
        if mode:
            label += f"_{mode}"
        return label

    def _sort_key(item):
        base_regime, mode, _, _ = item
        regime_order = {"kaon": 0, "charm": 1, "beauty": 2, "ew": 3, "combined": 4}
        mode_order = {None: 0, "direct": 1, "fromTau": 2}
        return (regime_order.get(base_regime, 99), mode_order.get(mode, 99), base_regime, mode or "")

    selected_by_mass = {}
    for key, items in files_by_mass.items():
        combined = [it for it in items if it[0] == "combined"]
        if combined:
            # Prefer combined_ff if it ever exists, otherwise take first.
            chosen = next((it for it in combined if it[2]), combined[0])

            selected = [(chosen[3], _label(chosen[0], chosen[1], chosen[2]))]

            if len(items) > 1:
                other_names = ", ".join(
                    _label(r, m, ff) for r, m, ff, _ in sorted(items, key=_sort_key) if r != "combined"
                )
                print(
                    f"[WARN] m={key[0]:.2f} {flavour}: found combined + other files ({other_names}); "
                    "using combined only."
                )

            selected_by_mass[key] = selected
            continue

        chosen = {}
        for base_regime, mode, is_ff, path in items:
            k2 = (base_regime, mode)
            if k2 not in chosen or is_ff:
                chosen[k2] = (base_regime, mode, is_ff, path)

        selected = [v for _, v in sorted(chosen.items(), key=lambda kv: _sort_key(kv[1]))]
        selected_by_mass[key] = [(path, _label(base_regime, mode, is_ff)) for base_regime, mode, is_ff, path in selected]

    files_by_mass = selected_by_mass

    # Report empty files only if no valid file exists for that mass at all
    # (i.e., the mass would be completely missing from analysis)
    masses_with_valid_files = {key[0] for key in files_by_mass.keys()}

    for f in empty_files:
        m = pattern.search(f.name)
        if m:
            mass_val = float(m.group(1).replace("p", "."))
            # Only warn if NO valid file exists for this mass
            if mass_val not in masses_with_valid_files:
                print(f"[SKIP] Empty file (no valid alternative): {f.name}")

    mass_points = sorted(files_by_mass.keys(), key=lambda x: x[0])
    print(f"Found {len(mass_points)} mass points")

    # Process each file
    if use_parallel:
        if n_workers is None:
            n_workers = multiprocessing.cpu_count()
        print(f"Using {n_workers} parallel workers")

        # Prepare arguments for parallel processing
        args_list = []
        for (mass_val, mass_str) in mass_points:
            sim_list = files_by_mass[(mass_val, mass_str)]
            args_list.append(
                (mass_val, mass_str, flavour, benchmark, lumi_fb, sim_list, dirac, separation_m, decay_seed)
            )

        # Process in parallel with tqdm progress bar
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            results = list(tqdm(
                executor.map(scan_single_mass_wrapper, args_list),
                total=len(args_list),
                desc=f"  {flavour}",
                unit="mass",
                ncols=80,
            ))

        # Print summary of results
        valid_results = [r for r in results if r is not None]
        excluded = sum(1 for r in valid_results if not np.isnan(r.get("eps2_min", np.nan)))
        print(f"  Completed: {excluded}/{len(valid_results)} mass points have sensitivity")

        # Filter out None results
        results = valid_results
    else:
        results = []
        for (mass_val, mass_str) in mass_points:
            sim_list = files_by_mass[(mass_val, mass_str)]
            res = scan_single_mass(
                mass_val,
                mass_str,
                flavour,
                benchmark,
                lumi_fb,
                sim_list,
                dirac=dirac,
                separation_m=separation_m,
                decay_seed=decay_seed,
            )
            if res:
                results.append(res)

    return pd.DataFrame(results)

def scan_single_mass_wrapper(args):
    """Wrapper for parallel processing - runs in quiet mode"""
    # Unpack args and add quiet=True for parallel execution
    (mass_val, mass_str, flavour, benchmark, lumi_fb, sim_list, dirac, separation_m, decay_seed) = args
    return scan_single_mass(
        mass_val, mass_str, flavour, benchmark, lumi_fb, sim_list,
        dirac=dirac, separation_m=separation_m, decay_seed=decay_seed, quiet=True
    )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate U² limits for HNL search")
    parser.add_argument("--parallel", action="store_true", help="Use parallel processing")
    parser.add_argument("--workers", type=int, default=None, help="Number of workers (default: all CPU cores)")
    parser.add_argument("--dirac", action="store_true", help="Dirac HNL interpretation (×2 yield vs Majorana)")
    parser.add_argument(
        "--separation-mm",
        type=float,
        default=1.0,
        help="Minimum charged-track separation at detector surface in mm (default: 1.0)",
    )
    parser.add_argument(
        "--decay-seed",
        type=int,
        default=12345,
        help="Random seed for decay sampling (default: 12345)",
    )
    args = parser.parse_args()
    if args.separation_mm <= 0:
        raise ValueError("--separation-mm must be positive.")

    n_workers = args.workers if args.workers else multiprocessing.cpu_count()
    mode_str = f"PARALLEL, {n_workers} workers" if args.parallel else "SINGLE-THREADED"
    hnl_type = "DIRAC" if args.dirac else "MAJORANA"
    print("="*60)
    print(f"U² LIMIT CALCULATOR ({mode_str}, {hnl_type})")
    print("="*60)

    GEOM_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    ANALYSIS_OUT_DIR.mkdir(parents=True, exist_ok=True)

    L_HL_LHC_FB = 3000.0
    results_out = ANALYSIS_OUT_DIR / "HNL_U2_limits_summary.csv"
    separation_m = args.separation_mm * 1e-3
    print(f"Decay separation: {args.separation_mm:.3f} mm (seed={args.decay_seed})")

    all_results = []

    # Process each flavour
    for flavour, benchmark in [("electron", "100"), ("muon", "010"), ("tau", "001")]:
        df = run_flavour(
            flavour,
            benchmark,
            L_HL_LHC_FB,
            use_parallel=args.parallel,
            n_workers=args.workers,
            dirac=args.dirac,
            separation_m=separation_m,
            decay_seed=args.decay_seed,
        )
        df["separation_mm"] = args.separation_mm
        all_results.append(df)

    # Combine and save
    final_df = pd.concat(all_results, ignore_index=True)
    final_df.to_csv(results_out, index=False)

    print(f"\n{'='*60}")
    print(f"COMPLETE!")
    print(f"Saved {len(final_df)} mass points to:")
    print(f"  {results_out}")
    print(f"{'='*60}")
