#!/usr/bin/env python
"""
U² limit calculator driver.

Usage:
    python run.py                       # Single-threaded (Majorana, default)
    python run.py --parallel            # Parallel processing (uses all cores)
    python run.py --parallel --workers 8  # Use 8 cores
    python run.py --dirac               # Dirac HNL interpretation (×2 yield)
    python run.py --separation-mm 1.0   # Minimum charged-track separation (mm)
    python run.py --flavour muon --mass 2.6  # Single mass point
    python run.py --hnlcalc-per-eps2    # Legacy: recompute HNLCalc per eps2 (slow)
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

MASS_FILTER_TOL = 5e-4  # GeV tolerance for --mass filtering

from geometry.per_parent_efficiency import build_drainage_gallery_mesh, preprocess_hnl_csv
from limits.expected_signal import expected_signal_events, couplings_from_eps2
from limits.timing_utils import _time_block


def _count(timing: dict | None, key: str, delta: int = 1) -> None:
    if timing is None:
        return
    timing[key] = timing.get(key, 0) + delta


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
    show_progress=None,
    timing_enabled=False,
    hnlcalc_per_eps2=False,
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
            dirac, separation_m, decay_seed, quiet, show_progress, timing_enabled,
            hnlcalc_per_eps2
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
    show_progress,
    timing_enabled,
    hnlcalc_per_eps2,
):
    """Internal implementation of scan_single_mass."""
    timing = {} if timing_enabled else None
    if timing is not None:
        timing["count_geom_files"] = len(sim_files)

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
            _count(timing, "count_geom_cache_hits")
            with _time_block(timing, "time_geom_load_s"):
                geom_df = pd.read_csv(geom_csv)
        else:
            _count(timing, "count_geom_cache_misses")
            if mesh is None:
                with _time_block(timing, "time_mesh_build_s"):
                    mesh = build_drainage_gallery_mesh()
            if not quiet:
                print(f"  Computing geometry for {sim_csv.name} (caching to {geom_csv.name})...")
            with _time_block(timing, "time_geom_compute_s"):
                geom_df = preprocess_hnl_csv(sim_csv, mesh, show_progress=show_progress)

            # Atomic write to prevent race condition in parallel execution
            with _time_block(timing, "time_geom_write_s"):
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

    with _time_block(timing, "time_geom_concat_s"):
        geom_df = pd.concat(geom_dfs, ignore_index=True)
    n_hits_total = geom_df['hits_tube'].sum() if 'hits_tube' in geom_df.columns else 0
    if timing is not None:
        timing["n_geom_rows"] = int(len(geom_df))
        timing["n_hits_total"] = int(n_hits_total)
    if not quiet:
        print(f"  Total combined: {len(geom_df)} HNLs, {n_hits_total} hit detector")

    if len(geom_df) == 0 or n_hits_total == 0:
        if not quiet:
            print(f"  WARNING: No hits, skipping")
        return None

    # 2. Scan |U|²
    eps2_scan = np.logspace(-12, -2, 100)
    if timing is not None:
        timing["n_eps2_points"] = int(len(eps2_scan))
    N_scan = []

    from decay.decay_detector import DecaySelection, build_decay_cache

    eps2_ref = 1e-6
    ctau0_ref = None
    br_ref = None
    if not hnlcalc_per_eps2:
        from models.hnl_model_hnlcalc import HNLModel
        with _time_block(timing, "time_hnlcalc_ref_s"):
            Ue2, Umu2, Utau2 = couplings_from_eps2(eps2_ref, benchmark)
            model = HNLModel(mass_GeV=mass_val, Ue2=Ue2, Umu2=Umu2, Utau2=Utau2)
            ctau0_ref = model.ctau0_m
            br_ref = model.production_brs()
        if timing is not None:
            timing["n_hnlcalc_eps2_ref"] = float(eps2_ref)

    with _time_block(timing, "time_decay_cache_s"):
        decay_cache = build_decay_cache(
            geom_df,
            mass_val,
            flavour,
            DecaySelection(separation_m=separation_m, seed=decay_seed),
            verbose=not quiet,
        )

    if not quiet:
        print("  Precomputing decay cache for separation scan...")
    with _time_block(timing, "time_eps2_scan_s"):
        for eps2 in eps2_scan:
            if hnlcalc_per_eps2:
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
                    timing=timing,
                )
            else:
                ctau0_m = ctau0_ref * (eps2_ref / eps2)
                br_scale = eps2 / eps2_ref
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
                    ctau0_m=ctau0_m,
                    br_per_parent=br_ref,
                    br_scale=br_scale,
                    timing=timing,
                )
            N_scan.append(N)

    N_scan = np.array(N_scan)

    # 3. Find exclusion range (N >= 2.996)
    mask_excluded = (N_scan >= 2.996)

    if not mask_excluded.any():
        if not quiet:
            print(f"  No sensitivity (peak = {N_scan.max():.1f})")
        result = {
            "mass_GeV": mass_val,
            "flavour": flavour,
            "benchmark": benchmark,
            "eps2_min": np.nan,
            "eps2_max": np.nan,
            "peak_events": N_scan.max()
        }
        if timing is not None:
            result.update(timing)
        return result

    indices_excl = np.where(mask_excluded)[0]
    i_lo = indices_excl[0]   # first index above threshold
    i_hi = indices_excl[-1]  # last index above threshold
    N_limit = 2.996

    # Log-linear interpolation at the lower crossing (eps2_min / red curve)
    if i_lo > 0:
        N_below, N_above = N_scan[i_lo - 1], N_scan[i_lo]
        dN = N_above - N_below
        if dN > 0:
            frac = np.clip((N_limit - N_below) / dN, 0.0, 1.0)
            log_lo = np.log10(eps2_scan[i_lo - 1])
            log_hi = np.log10(eps2_scan[i_lo])
            eps2_min = 10.0 ** (log_lo + frac * (log_hi - log_lo))
        else:
            eps2_min = eps2_scan[i_lo]
    else:
        eps2_min = eps2_scan[i_lo]

    # Log-linear interpolation at the upper crossing (eps2_max / blue curve)
    if i_hi < len(eps2_scan) - 1:
        N_above, N_below = N_scan[i_hi], N_scan[i_hi + 1]
        dN = N_above - N_below
        if dN > 0:
            frac = np.clip((N_above - N_limit) / dN, 0.0, 1.0)
            log_lo = np.log10(eps2_scan[i_hi])
            log_hi = np.log10(eps2_scan[i_hi + 1])
            eps2_max = 10.0 ** (log_lo + frac * (log_hi - log_lo))
        else:
            eps2_max = eps2_scan[i_hi]
    else:
        eps2_max = eps2_scan[i_hi]

    peak_events = N_scan.max()

    if not quiet:
        print(f"  ✓ Excluded: |U|² ∈ [{eps2_min:.2e}, {eps2_max:.2e}], peak = {peak_events:.0f}")

    result = {
        "mass_GeV": mass_val,
        "flavour": flavour,
        "benchmark": benchmark,
        "eps2_min": eps2_min,
        "eps2_max": eps2_max,
        "peak_events": peak_events
    }
    if timing is not None:
        result.update(timing)
    return result

def run_flavour(
    flavour,
    benchmark,
    lumi_fb,
    use_parallel=False,
    n_workers=None,
    dirac=False,
    separation_m=0.001,
    decay_seed=12345,
    show_progress=None,
    mass_filter=None,
    timing_enabled=False,
    hnlcalc_per_eps2=False,
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
        warn_for_key = True
        if mass_filter is not None:
            warn_for_key = abs(key[0] - mass_filter) <= MASS_FILTER_TOL
        combined = [it for it in items if it[0] == "combined"]
        if combined:
            # Prefer combined_ff if it ever exists, otherwise take first.
            chosen = next((it for it in combined if it[2]), combined[0])

            selected = [(chosen[3], _label(chosen[0], chosen[1], chosen[2]))]

            if warn_for_key and len(items) > 1:
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
    if mass_filter is not None:
        mass_points = [mp for mp in mass_points if abs(mp[0] - mass_filter) <= MASS_FILTER_TOL]
        if not mass_points:
            all_masses = sorted({m for m, _ in files_by_mass.keys()})
            preview = ", ".join(f"{m:.2f}" for m in all_masses[:10])
            more = "..." if len(all_masses) > 10 else ""
            print(f"[WARN] No mass points found for {flavour} at m={mass_filter:.4f} GeV.")
            print(f"       Available masses (first 10): {preview}{more}")
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
                (mass_val, mass_str, flavour, benchmark, lumi_fb, sim_list, dirac, separation_m, decay_seed, show_progress, timing_enabled, hnlcalc_per_eps2)
            )

        # Process in parallel (with optional progress bar)
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            iterator = executor.map(scan_single_mass_wrapper, args_list)
            if show_progress is None or show_progress:
                iterator = tqdm(
                    iterator,
                    total=len(args_list),
                    desc=f"  {flavour}",
                    unit="mass",
                    ncols=80,
                )
            results = list(iterator)

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
                show_progress=show_progress,
                timing_enabled=timing_enabled,
                hnlcalc_per_eps2=hnlcalc_per_eps2,
            )
            if res:
                results.append(res)

    return pd.DataFrame(results)

def scan_single_mass_wrapper(args):
    """Wrapper for parallel processing - runs in quiet mode"""
    # Unpack args and add quiet=True for parallel execution
    (mass_val, mass_str, flavour, benchmark, lumi_fb, sim_list, dirac, separation_m, decay_seed, show_progress, timing_enabled, hnlcalc_per_eps2) = args
    return scan_single_mass(
        mass_val, mass_str, flavour, benchmark, lumi_fb, sim_list,
        dirac=dirac, separation_m=separation_m, decay_seed=decay_seed, quiet=True,
        show_progress=show_progress, timing_enabled=timing_enabled,
        hnlcalc_per_eps2=hnlcalc_per_eps2,
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
    parser.add_argument(
        "--no-progress",
        action="store_true",
        help="Disable tqdm progress bars (auto-disabled in non-TTY environments)",
    )
    parser.add_argument(
        "--mass",
        type=float,
        default=None,
        help="Only process a single mass point in GeV (e.g., 2.6).",
    )
    parser.add_argument(
        "--flavour",
        type=str,
        choices=["electron", "muon", "tau"],
        default=None,
        help="Only process one flavour: electron, muon, or tau.",
    )
    parser.add_argument(
        "--timing",
        action="store_true",
        help="Record per-mass timing breakdown (adds time_* columns).",
    )
    parser.add_argument(
        "--timing-out",
        type=str,
        default=None,
        help="Optional path for timing CSV (default: output/csv/analysis/HNL_U2_timing.csv).",
    )
    parser.add_argument(
        "--hnlcalc-per-eps2",
        action="store_true",
        help="Recompute HNLCalc for every eps2 point (slow, legacy behavior).",
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

    # Determine progress bar visibility: explicit --no-progress overrides auto-detect
    show_progress = None if not args.no_progress else False
    timing_enabled = args.timing

    # Process each flavour
    flavour_list = [("electron", "100"), ("muon", "010"), ("tau", "001")]
    if args.flavour:
        flavour_list = [fb for fb in flavour_list if fb[0] == args.flavour]
    for flavour, benchmark in flavour_list:
        df = run_flavour(
            flavour,
            benchmark,
            L_HL_LHC_FB,
            use_parallel=args.parallel,
            n_workers=args.workers,
            dirac=args.dirac,
            separation_m=separation_m,
            decay_seed=args.decay_seed,
            show_progress=show_progress,
            mass_filter=args.mass,
            timing_enabled=timing_enabled,
            hnlcalc_per_eps2=args.hnlcalc_per_eps2,
        )
        df["separation_mm"] = args.separation_mm
        all_results.append(df)

    # Combine and save
    final_df = pd.concat(all_results, ignore_index=True)
    final_df.to_csv(results_out, index=False)

    if timing_enabled:
        timing_out = Path(args.timing_out) if args.timing_out else (ANALYSIS_OUT_DIR / "HNL_U2_timing.csv")
        timing_cols = [c for c in final_df.columns if c.startswith("time_") or c.startswith("count_") or c.startswith("n_")]
        timing_df = final_df[["mass_GeV", "flavour", "benchmark"] + timing_cols]
        timing_df.to_csv(timing_out, index=False)

    print(f"\n{'='*60}")
    print(f"COMPLETE!")
    print(f"Saved {len(final_df)} mass points to:")
    print(f"  {results_out}")
    if timing_enabled:
        print(f"Timing breakdown saved to:")
        print(f"  {timing_out}")
    print(f"{'='*60}")
