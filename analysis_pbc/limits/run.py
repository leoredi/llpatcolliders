#!/usr/bin/env python
"""
Limit calculator driver.

Usage:
    python run.py                       # Single-threaded (Majorana, default)
    python run.py --parallel            # Parallel processing (uses all cores)
    python run.py --parallel --workers 8  # Use 8 cores
    python run.py --dirac               # Dirac HNL interpretation (×2 yield)

ALP mode:
    python run.py --particle alp --alp-benchmark BC10
"""

import sys
import re
import os
import tempfile
from pathlib import Path
import numpy as np
import pandas as pd
import argparse
from concurrent.futures import ProcessPoolExecutor
import multiprocessing

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
from limits.expected_signal import expected_signal_events, expected_signal_events_alp

def scan_single_mass(mass_val, mass_str, flavour, benchmark, lumi_fb, sim_files, dirac=False):
    """
    Process one mass point, combining all available production regimes (kaon/charm/beauty/ew).
    sim_files: list of tuples (sim_csv_path, regime)
    dirac: if True, multiply yield by 2 for Dirac HNL interpretation
    """
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

        if geom_csv.exists():
            geom_df = pd.read_csv(geom_csv)
        else:
            if mesh is None:
                mesh = build_drainage_gallery_mesh()
            print(f"  Computing geometry for {sim_csv.name} (caching to {geom_csv.name})...")
            geom_df = preprocess_hnl_csv(sim_csv, mesh)

            # Atomic write to prevent race condition in parallel execution
            with tempfile.NamedTemporaryFile(mode='w', dir=geom_csv.parent,
                                              suffix='.tmp', delete=False) as tmp:
                geom_df.to_csv(tmp.name, index=False)
            os.replace(tmp.name, geom_csv)

        n_hits = geom_df['hits_tube'].sum() if 'hits_tube' in geom_df.columns else 0
        print(f"  Loaded {len(geom_df)} HNLs from {regime}, {n_hits} hit detector")
        geom_dfs.append(geom_df)

    if len(geom_dfs) == 0:
        print(f"  WARNING: No geometry loaded, skipping")
        return None

    geom_df = pd.concat(geom_dfs, ignore_index=True)
    n_hits_total = geom_df['hits_tube'].sum() if 'hits_tube' in geom_df.columns else 0
    print(f"  Total combined: {len(geom_df)} HNLs, {n_hits_total} hit detector")

    if len(geom_df) == 0 or n_hits_total == 0:
        print(f"  WARNING: No hits, skipping")
        return None

    # 2. Scan |U|²
    eps2_scan = np.logspace(-12, -2, 100)
    N_scan = []

    for eps2 in eps2_scan:
        N = expected_signal_events(geom_df, mass_val, eps2, benchmark, lumi_fb, dirac=dirac)
        N_scan.append(N)

    N_scan = np.array(N_scan)

    # 3. Find exclusion range (N >= 2.996)
    mask_excluded = (N_scan >= 2.996)

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

def run_flavour(flavour, benchmark, lumi_fb, use_parallel=False, n_workers=None, dirac=False):
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

            # Transitional support for legacy combined files:
            # If tau fromTau files exist alongside a combined file, include them only
            # if the combined file does NOT already contain tau-parent (PDG 15) rows.
            # This prevents silently dropping τ→NX samples while still avoiding
            # double counting when users run combine with --no-cleanup.
            if flavour == "tau":
                fromtau_candidates = [it for it in items if it[1] == "fromTau" and it[0] != "combined"]
                if fromtau_candidates:
                    import csv

                    def _combined_has_tau_parent(csv_path: Path) -> bool:
                        try:
                            with csv_path.open(newline="") as f:
                                reader = csv.DictReader(f)
                                for row in reader:
                                    try:
                                        pid = int(row["parent_pdg"])
                                    except Exception:
                                        continue
                                    if abs(pid) == 15:
                                        return True
                        except Exception:
                            # If we can't inspect it reliably, be conservative and assume it has tau.
                            return True
                        return False

                    has_tau = _combined_has_tau_parent(chosen[3])
                    if not has_tau:
                        extra = {}
                        for base_regime, mode, is_ff, path in fromtau_candidates:
                            k2 = (base_regime, mode)
                            if k2 not in extra or is_ff:
                                extra[k2] = (base_regime, mode, is_ff, path)
                        extras_sorted = [v for _, v in sorted(extra.items(), key=lambda kv: _sort_key(kv[1]))]
                        selected.extend(
                            [(path, _label(base_regime, mode, is_ff)) for base_regime, mode, is_ff, path in extras_sorted]
                        )
                        extra_names = ", ".join(_label(r, m, ff) for r, m, ff, _ in extras_sorted)
                        print(
                            f"[INFO] m={key[0]:.2f} tau: combined file has no PDG 15 rows; "
                            f"including additional fromTau file(s): {extra_names}"
                        )
                    else:
                        print(
                            f"[WARN] m={key[0]:.2f} tau: found fromTau file(s) but combined already contains PDG 15 rows; "
                            "skipping extra fromTau files to avoid double counting."
                        )

            if len(items) > 1:
                other_names = ", ".join(
                    _label(r, m, ff) for r, m, ff, _ in sorted(items, key=_sort_key) if r != "combined"
                )
                print(
                    f"[WARN] m={key[0]:.2f} {flavour}: found combined + other files ({other_names}); "
                    "using combined (and any needed fromTau) to avoid double counting."
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
            args_list.append((mass_val, mass_str, flavour, benchmark, lumi_fb, sim_list, dirac))

        # Process in parallel
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            results = list(executor.map(scan_single_mass_wrapper, args_list))

        # Filter out None results
        results = [r for r in results if r is not None]
    else:
        results = []
        for (mass_val, mass_str) in mass_points:
            sim_list = files_by_mass[(mass_val, mass_str)]
            res = scan_single_mass(mass_val, mass_str, flavour, benchmark, lumi_fb, sim_list, dirac=dirac)
            if res:
                results.append(res)

    return pd.DataFrame(results)

def scan_single_mass_wrapper(args):
    """Wrapper for parallel processing"""
    return scan_single_mass(*args)


def _alp_infer_channel_from_stem(stem: str) -> str | None:
    s = stem.lower()
    mapping = [
        ("h_to_aa", "h_to_aa"),
        ("htoa a", "h_to_aa"),
        ("htoa", "h_to_aa"),
        ("h_to_za", "h_to_Za"),
        ("htoz a", "h_to_Za"),
        ("htoza", "h_to_Za"),
        ("z_to_gamma_a", "Z_to_gamma_a"),
        ("ztogammaa", "Z_to_gamma_a"),
        ("z_to_gam a", "Z_to_gamma_a"),
        ("b_to_ka", "B_to_Ka"),
        ("btoka", "B_to_Ka"),
    ]
    for needle, mode in mapping:
        if needle.replace(" ", "") in s.replace(" ", ""):
            return mode
    return None


def _alp_collect_files(sim_dir: Path, benchmark: str) -> dict[tuple[float, str], list[tuple[Path, str]]]:
    """
    Return {(mass_val, mass_str): [(path, production_mode), ...]} for ALP samples.

    Expected filename convention (flexible):
        ALP_{mass}GeV_{benchmark}_{channel}.csv
    where {mass} uses p as decimal separator (e.g. 3p00), and channel contains
    one of: h_to_aa, h_to_Za, Z_to_gamma_a, B_to_Ka.
    """
    pattern = re.compile(rf"^ALP_([0-9]+p[0-9]{{1,2}})GeV_({benchmark})_.*\.csv$", re.IGNORECASE)
    files_by_mass: dict[tuple[float, str], list[tuple[Path, str]]] = {}

    for f in sim_dir.glob("ALP_*.csv"):
        m = pattern.search(f.name)
        if not m:
            continue
        mass_str = m.group(1)
        mass_val = float(mass_str.replace("p", "."))
        channel = _alp_infer_channel_from_stem(f.stem)
        if channel is None:
            raise ValueError(
                f"Could not infer ALP production channel from filename: {f.name}. "
                "Include one of: h_to_aa, h_to_Za, Z_to_gamma_a, B_to_Ka."
            )
        files_by_mass.setdefault((mass_val, mass_str), []).append((f, channel))

    return files_by_mass


def scan_single_mass_alp(
    mass_val: float,
    mass_str: str,
    benchmark: str,
    lumi_fb: float,
    sim_files: list[tuple[Path, str]],
    fa_min: float,
    fa_max: float,
    fa_points: int,
    br_h_aa: float | None,
    lambda_aha_GeV: float | None,
    br_h_Za: float | None,
    C_aZh: float | None,
    Lambda_GeV: float,
    C_gamma_Z: float,
    C_WW: float,
) -> dict | None:
    print(f"\n[ALP {benchmark} {mass_val} GeV] Processing ({len(sim_files)} production file(s))...")

    geom_samples: list[tuple[pd.DataFrame, str]] = []
    mesh = None

    for sim_csv, channel in sim_files:
        geom_cache_name = f"{sim_csv.stem}_geom.csv"
        geom_csv = GEOM_CACHE_DIR / geom_cache_name

        if geom_csv.exists():
            geom_df = pd.read_csv(geom_csv)
        else:
            if mesh is None:
                mesh = build_drainage_gallery_mesh()
            print(f"  Computing geometry for {sim_csv.name} (caching to {geom_csv.name})...")
            geom_df = preprocess_hnl_csv(sim_csv, mesh)
            with tempfile.NamedTemporaryFile(mode="w", dir=geom_csv.parent, suffix=".tmp", delete=False) as tmp:
                geom_df.to_csv(tmp.name, index=False)
            os.replace(tmp.name, geom_csv)

        n_hits = geom_df["hits_tube"].sum() if "hits_tube" in geom_df.columns else 0
        print(f"  Loaded {len(geom_df)} ALPs from {channel}, {n_hits} hit detector")
        if len(geom_df) == 0:
            continue
        geom_samples.append((geom_df, channel))

    if not geom_samples:
        print("  WARNING: No geometry loaded, skipping")
        return None

    fa_scan = np.logspace(np.log10(fa_min), np.log10(fa_max), int(fa_points))
    N_scan = np.zeros_like(fa_scan, dtype=float)

    for i, fa in enumerate(fa_scan):
        N_total = 0.0
        for geom_df, channel in geom_samples:
            N_total += expected_signal_events_alp(
                geom_df=geom_df,
                mass_GeV=mass_val,
                fa_GeV=float(fa),
                benchmark=benchmark,
                lumi_fb=lumi_fb,
                production_mode=channel,
                br_h_aa=br_h_aa,
                lambda_aha_GeV=lambda_aha_GeV,
                br_h_Za=br_h_Za,
                C_aZh=C_aZh,
                Lambda_GeV=Lambda_GeV,
                C_gamma_Z=C_gamma_Z,
                C_WW=C_WW,
            )
        N_scan[i] = float(N_total)

    mask_excluded = N_scan >= 2.996
    peak_events = float(np.max(N_scan)) if len(N_scan) else float("nan")

    if not mask_excluded.any():
        print(f"  No sensitivity (peak = {peak_events:.2e})")
        return {
            "mass_GeV": mass_val,
            "benchmark": benchmark,
            "fa_min": np.nan,
            "fa_max": np.nan,
            "peak_events": peak_events,
        }

    idx = np.where(mask_excluded)[0]
    fa_min_excl = float(fa_scan[idx[0]])
    fa_max_excl = float(fa_scan[idx[-1]])
    print(f"  ✓ Excluded: f_a ∈ [{fa_min_excl:.2e}, {fa_max_excl:.2e}], peak = {peak_events:.2e}")

    return {
        "mass_GeV": mass_val,
        "benchmark": benchmark,
        "fa_min": fa_min_excl,
        "fa_max": fa_max_excl,
        "peak_events": peak_events,
    }


def scan_single_mass_alp_wrapper(args):
    return scan_single_mass_alp(*args)


def run_alp_benchmark(
    benchmark: str,
    lumi_fb: float,
    *,
    use_parallel: bool = False,
    n_workers: int | None = None,
    fa_min: float = 1e2,
    fa_max: float = 1e12,
    fa_points: int = 50,
    br_h_aa: float | None = None,
    lambda_aha_GeV: float | None = None,
    br_h_Za: float | None = None,
    C_aZh: float | None = None,
    Lambda_GeV: float = 1000.0,
    C_gamma_Z: float = 1.0,
    C_WW: float = 1.0,
) -> pd.DataFrame:
    print(f"\n{'='*60}")
    print(f"ALP BENCHMARK: {benchmark}")
    print(f"{'='*60}")

    files_by_mass = _alp_collect_files(SIM_DIR, benchmark)
    mass_points = sorted(files_by_mass.keys(), key=lambda x: x[0])
    print(f"Found {len(mass_points)} mass points")

    # Early diagnostics for channels that need explicit BR/coupling parameters.
    channels_present = {ch for items in files_by_mass.values() for _, ch in items}
    if "h_to_aa" in channels_present and br_h_aa is None and lambda_aha_GeV is None:
        print("[WARN] Found ALP h_to_aa samples but neither --br-h-aa nor --lambda-aha is set; h_to_aa will contribute 0.")
    if "h_to_Za" in channels_present and br_h_Za is None and C_aZh is None:
        print("[WARN] Found ALP h_to_Za samples but neither --br-h-Za nor --C-aZh is set; h_to_Za will contribute 0.")

    if use_parallel:
        if n_workers is None:
            n_workers = multiprocessing.cpu_count()
        print(f"Using {n_workers} parallel workers")
        args_list = []
        for (mass_val, mass_str) in mass_points:
            sim_list = files_by_mass[(mass_val, mass_str)]
            args_list.append(
                (
                    mass_val,
                    mass_str,
                    benchmark,
                    lumi_fb,
                    sim_list,
                    fa_min,
                    fa_max,
                    fa_points,
                    br_h_aa,
                    lambda_aha_GeV,
                    br_h_Za,
                    C_aZh,
                    Lambda_GeV,
                    C_gamma_Z,
                    C_WW,
                )
            )
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            results = list(executor.map(scan_single_mass_alp_wrapper, args_list))
        results = [r for r in results if r is not None]
    else:
        results = []
        for (mass_val, mass_str) in mass_points:
            sim_list = files_by_mass[(mass_val, mass_str)]
            res = scan_single_mass_alp(
                mass_val,
                mass_str,
                benchmark,
                lumi_fb,
                sim_list,
                fa_min=fa_min,
                fa_max=fa_max,
                fa_points=fa_points,
                br_h_aa=br_h_aa,
                lambda_aha_GeV=lambda_aha_GeV,
                br_h_Za=br_h_Za,
                C_aZh=C_aZh,
                Lambda_GeV=Lambda_GeV,
                C_gamma_Z=C_gamma_Z,
                C_WW=C_WW,
            )
            if res:
                results.append(res)

    if not results:
        return pd.DataFrame(columns=["mass_GeV", "benchmark", "fa_min", "fa_max", "peak_events"])

    return pd.DataFrame(results)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate LLP limits for HNL/ALP searches")
    parser.add_argument("--particle", choices=["hnl", "alp"], default="hnl", help="Particle hypothesis")
    parser.add_argument("--parallel", action="store_true", help="Use parallel processing")
    parser.add_argument("--workers", type=int, default=None, help="Number of workers (default: all CPU cores)")
    parser.add_argument("--dirac", action="store_true", help="Dirac HNL interpretation (×2 yield vs Majorana)")
    parser.add_argument("--alp-benchmark", default="BC10", help="ALP benchmark (BC9/BC10/BC11/all)")
    parser.add_argument("--fa-min", type=float, default=1e2, help="ALP scan: minimum f_a [GeV]")
    parser.add_argument("--fa-max", type=float, default=1e12, help="ALP scan: maximum f_a [GeV]")
    parser.add_argument("--fa-points", type=int, default=50, help="ALP scan: number of f_a grid points")
    parser.add_argument("--br-h-aa", type=float, default=None, help="Override BR(h→aa) (used for h_to_aa samples)")
    parser.add_argument("--lambda-aha", type=float, default=None, help="Higgs portal coupling λ_aha [GeV] (used for h_to_aa samples)")
    parser.add_argument("--br-h-Za", type=float, default=None, help="Override BR(h→Za) (used for h_to_Za samples)")
    parser.add_argument("--C-aZh", type=float, default=None, help="Effective coupling C_aZh (used for h_to_Za samples)")
    parser.add_argument("--Lambda", type=float, default=1000.0, help="Scale Λ [GeV] (used for h_to_Za samples)")
    parser.add_argument("--C-gamma-Z", type=float, default=1.0, help="Wilson coefficient C_γZ (used for Z_to_gamma_a samples)")
    parser.add_argument("--C-WW", type=float, default=1.0, help="Wilson coefficient C_WW (used for B_to_Ka samples)")
    args = parser.parse_args()

    n_workers = args.workers if args.workers else multiprocessing.cpu_count()
    mode_str = f"PARALLEL, {n_workers} workers" if args.parallel else "SINGLE-THREADED"
    print("="*60)
    if args.particle == "hnl":
        hnl_type = "DIRAC" if args.dirac else "MAJORANA"
        print(f"U² LIMIT CALCULATOR ({mode_str}, {hnl_type})")
    else:
        print(f"ALP LIMIT CALCULATOR ({mode_str})")
    print("="*60)

    GEOM_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    ANALYSIS_OUT_DIR.mkdir(parents=True, exist_ok=True)

    L_HL_LHC_FB = 3000.0

    if args.particle == "hnl":
        results_out = ANALYSIS_OUT_DIR / "HNL_U2_limits_summary.csv"
        all_results = []
        for flavour, benchmark in [("electron", "100"), ("muon", "010"), ("tau", "001")]:
            df = run_flavour(
                flavour,
                benchmark,
                L_HL_LHC_FB,
                use_parallel=args.parallel,
                n_workers=args.workers,
                dirac=args.dirac,
            )
            all_results.append(df)
        final_df = pd.concat(all_results, ignore_index=True)
        final_df.to_csv(results_out, index=False)
    else:
        results_out = ANALYSIS_OUT_DIR / "ALP_fa_limits_summary.csv"
        bench = str(args.alp_benchmark).upper()
        if bench == "ALL":
            parts = []
            for b in ["BC9", "BC10", "BC11"]:
                df_b = run_alp_benchmark(
                    b,
                    L_HL_LHC_FB,
                    use_parallel=args.parallel,
                    n_workers=args.workers,
                    fa_min=args.fa_min,
                    fa_max=args.fa_max,
                    fa_points=args.fa_points,
                    br_h_aa=args.br_h_aa,
                    lambda_aha_GeV=args.lambda_aha,
                    br_h_Za=args.br_h_Za,
                    C_aZh=args.C_aZh,
                    Lambda_GeV=args.Lambda,
                    C_gamma_Z=args.C_gamma_Z,
                    C_WW=args.C_WW,
                )
                parts.append(df_b)
            final_df = pd.concat(parts, ignore_index=True)
        else:
            final_df = run_alp_benchmark(
                bench,
                L_HL_LHC_FB,
                use_parallel=args.parallel,
                n_workers=args.workers,
                fa_min=args.fa_min,
                fa_max=args.fa_max,
                fa_points=args.fa_points,
                br_h_aa=args.br_h_aa,
                lambda_aha_GeV=args.lambda_aha,
                br_h_Za=args.br_h_Za,
                C_aZh=args.C_aZh,
                Lambda_GeV=args.Lambda,
                C_gamma_Z=args.C_gamma_Z,
                C_WW=args.C_WW,
            )
        final_df.to_csv(results_out, index=False)

    print(f"\n{'='*60}")
    print(f"COMPLETE!")
    print(f"Saved {len(final_df)} mass points to:")
    print(f"  {results_out}")
    print(f"{'='*60}")
