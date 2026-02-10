#!/usr/bin/env python

import sys
import re
import os
import io
import json
import tempfile
import contextlib
from pathlib import Path
import numpy as np
import pandas as pd
import argparse
from concurrent.futures import ProcessPoolExecutor
import multiprocessing
from tqdm import tqdm

THIS_FILE = Path(__file__).resolve()
ANALYSIS_DIR = THIS_FILE.parent
REPO_ROOT = ANALYSIS_DIR.parents[1]
OUTPUT_DIR = REPO_ROOT / "output" / "csv"

SIM_DIR = OUTPUT_DIR / "simulation"
GEOM_CACHE_DIR = OUTPUT_DIR / "geometry"
ANALYSIS_OUT_DIR = OUTPUT_DIR / "analysis"

ANALYSIS_ROOT = ANALYSIS_DIR.parent
if str(ANALYSIS_ROOT) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_ROOT))

MASS_FILTER_TOL = 5e-4

from geometry.per_parent_efficiency import build_drainage_gallery_mesh, preprocess_hnl_csv
from limits.expected_signal import expected_signal_events, couplings_from_eps2
from limits.timing_utils import _time_block


def _count(timing: dict | None, key: str, delta: int = 1) -> None:
    if timing is None:
        return
    timing[key] = timing.get(key, 0) + delta


def _load_sim_metadata(sim_csv: Path) -> dict:
    meta_path = Path(f"{sim_csv}.meta.json")
    if not meta_path.exists():
        return {}
    try:
        with meta_path.open("r", encoding="utf-8") as f:
            data = json.load(f)
            if isinstance(data, dict):
                return data
    except Exception as exc:
        print(f"[WARN] Failed to read metadata sidecar {meta_path.name}: {exc}")
    return {}


def _attach_sim_metadata(geom_df: pd.DataFrame, sim_csv: Path) -> pd.DataFrame:
    if all(col in geom_df.columns for col in ("qcd_mode", "sigma_gen_pb", "pthat_min_gev")):
        return geom_df

    meta = _load_sim_metadata(sim_csv)
    if not meta:
        return geom_df

    if "qcd_mode" not in geom_df.columns:
        geom_df["qcd_mode"] = str(meta.get("qcd_mode", "auto"))
    if "sigma_gen_pb" not in geom_df.columns:
        sigma_pb = pd.to_numeric(pd.Series([meta.get("sigma_gen_pb")]), errors="coerce").iloc[0]
        geom_df["sigma_gen_pb"] = float(sigma_pb) if np.isfinite(sigma_pb) else np.nan
    if "pthat_min_gev" not in geom_df.columns:
        pthat = pd.to_numeric(pd.Series([meta.get("pthat_min_gev")]), errors="coerce").iloc[0]
        geom_df["pthat_min_gev"] = float(pthat) if np.isfinite(pthat) else np.nan
    return geom_df


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
    p_min_GeV=0.5,
    reco_efficiency=1.0,
    quiet=False,
    show_progress=None,
    timing_enabled=False,
    hnlcalc_per_eps2=False,
):
    stdout_ctx = contextlib.redirect_stdout(io.StringIO()) if quiet else contextlib.nullcontext()

    with stdout_ctx:
        return _scan_single_mass_impl(
            mass_val, mass_str, flavour, benchmark, lumi_fb, sim_files,
            dirac, separation_m, decay_seed, p_min_GeV, reco_efficiency,
            quiet, show_progress, timing_enabled, hnlcalc_per_eps2
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
    p_min_GeV,
    reco_efficiency,
    quiet,
    show_progress,
    timing_enabled,
    hnlcalc_per_eps2,
):
    timing = {} if timing_enabled else None
    if timing is not None:
        timing["count_geom_files"] = len(sim_files)

    if not quiet:
        print(f"\n[{flavour} {mass_val} GeV] Processing ({len(sim_files)} production file(s))...")

    geom_dfs = []
    mesh = None

    for sim_csv, regime in sim_files:
        geom_cache_name = f"{sim_csv.stem}_geom.csv"
        geom_csv = GEOM_CACHE_DIR / geom_cache_name

        legacy_geom_csv = GEOM_CACHE_DIR / f"HNL_{mass_str}GeV_{flavour}_geom.csv"
        if not geom_csv.exists() and legacy_geom_csv.exists():
            geom_csv = legacy_geom_csv

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

            with _time_block(timing, "time_geom_write_s"):
                with tempfile.NamedTemporaryFile(mode='w', dir=geom_csv.parent,
                                                  suffix='.tmp', delete=False) as tmp:
                    geom_df.to_csv(tmp.name, index=False)
                os.replace(tmp.name, geom_csv)

        geom_df = _attach_sim_metadata(geom_df, sim_csv)

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

    eps2_scan = np.logspace(-12, -2, 100)
    if timing is not None:
        timing["n_eps2_points"] = int(len(eps2_scan))
    N_scan = []

    from decay.decay_detector import DecaySelection, build_decay_cache

    eps2_ref = 1e-6
    ctau0_ref = None
    br_ref = None
    if not hnlcalc_per_eps2:
        from hnl_models.hnl_model_hnlcalc import HNLModel
        with _time_block(timing, "time_hnlcalc_ref_s"):
            Ue2, Umu2, Utau2 = couplings_from_eps2(eps2_ref, benchmark)
            model = HNLModel(mass_GeV=mass_val, Ue2=Ue2, Umu2=Umu2, Utau2=Utau2)
            ctau0_ref = model.ctau0_m
            br_ref = model.production_brs()
            br_vis = model.visible_branching_ratio()
        if not quiet:
            print(f"  BR_vis(≥2 charged) = {br_vis:.3f} [implicit in MC sampling, not applied as factor]")
        if timing is not None:
            timing["n_hnlcalc_eps2_ref"] = float(eps2_ref)
            timing["br_visible"] = br_vis

    with _time_block(timing, "time_decay_cache_s"):
        decay_cache = build_decay_cache(
            geom_df,
            mass_val,
            flavour,
            DecaySelection(separation_m=separation_m, seed=decay_seed, p_min_GeV=p_min_GeV),
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
                    p_min_GeV=p_min_GeV,
                    reco_efficiency=reco_efficiency,
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
                    p_min_GeV=p_min_GeV,
                    reco_efficiency=reco_efficiency,
                    decay_cache=decay_cache,
                    ctau0_m=ctau0_m,
                    br_per_parent=br_ref,
                    br_scale=br_scale,
                    timing=timing,
                )
            N_scan.append(N)

    N_scan = np.array(N_scan)

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
    i_lo = indices_excl[0]
    i_hi = indices_excl[-1]
    N_limit = 2.996

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
    p_min_GeV=0.5,
    reco_efficiency=1.0,
    show_progress=None,
    mass_filter=None,
    timing_enabled=False,
    hnlcalc_per_eps2=False,
    allow_variant_drop=False,
    max_mass=None,
):
    print(f"\n{'='*60}")
    print(f"FLAVOUR: {flavour.upper()} (Benchmark {benchmark})")
    print(f"{'='*60}")

    pattern = re.compile(
        rf"^HNL_([0-9]+p[0-9]{{1,2}})GeV_{flavour}_"
        r"((?:kaon|charm|beauty|Bc|ew|combined)(?:_ff)?)"
        r"(?:_(direct|fromTau))?"
        r"(?:_(hardBc|hardccbar|hardbbbar)(?:_pTHat([0-9]+(?:p[0-9]+)?))?)?"
        r"\.csv$"
    )

    files = []
    empty_files = []
    for f in SIM_DIR.glob(f"*{flavour}*.csv"):
        match = pattern.search(f.name)
        if not match:
            continue

        if f.stat().st_size < 1000:
            empty_files.append(f)
            continue

        mass_str = match.group(1)
        mass_val = float(mass_str.replace("p", "."))
        regime_token = match.group(2)
        mode = match.group(3)
        qcd_mode = match.group(4) or "auto"
        pthat_token = match.group(5)
        pthat_min = float(pthat_token.replace("p", ".")) if pthat_token else None

        is_ff = regime_token.endswith("_ff")
        base_regime = regime_token.replace("_ff", "")

        files.append((mass_val, mass_str, base_regime, mode, is_ff, qcd_mode, pthat_min, f))

    files_by_mass = {}
    for mass_val, mass_str, base_regime, mode, is_ff, qcd_mode, pthat_min, path in files:
        key = (mass_val, mass_str)
        files_by_mass.setdefault(key, []).append((base_regime, mode, is_ff, qcd_mode, pthat_min, path))

    def _label(
        base_regime: str,
        mode: str | None,
        is_ff: bool,
        qcd_mode: str,
        pthat_min: float | None,
    ) -> str:
        label = base_regime
        if is_ff:
            label += "_ff"
        if mode:
            label += f"_{mode}"
        if qcd_mode != "auto":
            label += f"_{qcd_mode}"
            if pthat_min is not None:
                pthat_text = f"{pthat_min:g}".replace(".", "p")
                label += f"_pTHat{pthat_text}"
        return label

    def _variant_priority(
        base_regime: str,
        is_ff: bool,
        qcd_mode: str,
        pthat_min: float | None,
    ) -> tuple[int, int, float]:
        # Prefer hard-QCD sliced samples for heavy-flavor regimes when available.
        if base_regime == "charm" and qcd_mode == "hardccbar":
            qcd_priority = 3
        elif base_regime in {"beauty", "Bc"} and qcd_mode in {"hardbbbar", "hardBc"}:
            qcd_priority = 3
        elif qcd_mode != "auto":
            qcd_priority = 2
        else:
            qcd_priority = 1
        ff_priority = 1 if is_ff else 0
        pthat_priority = float(pthat_min) if pthat_min is not None else -1.0
        return (qcd_priority, ff_priority, pthat_priority)

    def _sort_key(item):
        base_regime, mode, is_ff, qcd_mode, pthat_min, _ = item
        regime_order = {"kaon": 0, "charm": 1, "beauty": 2, "Bc": 3, "ew": 4, "combined": 5}
        mode_order = {None: 0, "direct": 1, "fromTau": 2}
        variant_order = _variant_priority(base_regime, is_ff, qcd_mode, pthat_min)
        return (regime_order.get(base_regime, 99), mode_order.get(mode, 99), base_regime, mode or "", variant_order)

    selected_by_mass = {}
    combined_coexist_count = 0
    for key, items in files_by_mass.items():
        warn_for_key = True
        if mass_filter is not None:
            warn_for_key = abs(key[0] - mass_filter) <= MASS_FILTER_TOL
        combined = [it for it in items if it[0] == "combined"]
        if combined:
            chosen = max(
                combined,
                key=lambda it: _variant_priority(it[0], it[2], it[3], it[4]),
            )

            selected = [(chosen[5], _label(chosen[0], chosen[1], chosen[2], chosen[3], chosen[4]))]

            if warn_for_key and len(items) > 1:
                n_other = sum(1 for it in items if it[0] != "combined")
                combined_coexist_count = combined_coexist_count + n_other

            selected_by_mass[key] = selected
            continue

        chosen = {}
        all_candidates_for_key = {}
        for base_regime, mode, is_ff, qcd_mode, pthat_min, path in items:
            k2 = (base_regime, mode)
            if k2 not in all_candidates_for_key:
                all_candidates_for_key[k2] = []
            all_candidates_for_key[k2].append((base_regime, mode, is_ff, qcd_mode, pthat_min, path))
            if k2 not in chosen:
                chosen[k2] = (base_regime, mode, is_ff, qcd_mode, pthat_min, path)
                continue
            current = chosen[k2]
            new_priority = _variant_priority(base_regime, is_ff, qcd_mode, pthat_min)
            old_priority = _variant_priority(current[0], current[2], current[3], current[4])
            if new_priority > old_priority:
                chosen[k2] = (base_regime, mode, is_ff, qcd_mode, pthat_min, path)

        # Check for dropped variants (compare by path, not identity).
        # Only enforce for masses matching --mass filter (warn_for_key);
        # unrelated masses should not block the run.
        if warn_for_key:
            for k2, candidates in all_candidates_for_key.items():
                if len(candidates) > 1:
                    kept_path = str(chosen[k2][5])
                    dropped = [c for c in candidates if str(c[5]) != kept_path]
                    if dropped:
                        msg = (
                            f"Multiple variants for {k2} at m={key[0]:.2f}: "
                            f"keeping {_label(*chosen[k2][:5])}, would drop "
                            f"{[_label(*d[:5]) for d in dropped]}. "
                            f"Pass --allow-variant-drop to override."
                        )
                        if allow_variant_drop:
                            print(f"[WARN] {msg}")
                        else:
                            raise ValueError(msg)

        selected = [v for _, v in sorted(chosen.items(), key=lambda kv: _sort_key(kv[1]))]
        selected_by_mass[key] = [
            (path, _label(base_regime, mode, is_ff, qcd_mode, pthat_min))
            for base_regime, mode, is_ff, qcd_mode, pthat_min, path in selected
        ]

    files_by_mass = selected_by_mass

    if combined_coexist_count > 0:
        print(
            f"[INFO] {combined_coexist_count} individual regime files coexist with combined files; "
            "using combined only."
        )

    masses_with_valid_files = {key[0] for key in files_by_mass.keys()}

    for f in empty_files:
        m = pattern.search(f.name)
        if m:
            mass_val = float(m.group(1).replace("p", "."))
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
    if max_mass is not None:
        mass_points = [mp for mp in mass_points if mp[0] <= max_mass]
    print(f"Found {len(mass_points)} mass points")

    if use_parallel:
        if n_workers is None:
            n_workers = multiprocessing.cpu_count()
        print(f"Using {n_workers} parallel workers")

        args_list = []
        for (mass_val, mass_str) in mass_points:
            sim_list = files_by_mass[(mass_val, mass_str)]
            args_list.append(
                (mass_val, mass_str, flavour, benchmark, lumi_fb, sim_list, dirac, separation_m, decay_seed, p_min_GeV, reco_efficiency, show_progress, timing_enabled, hnlcalc_per_eps2)
            )

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

        valid_results = [r for r in results if r is not None]
        excluded = sum(1 for r in valid_results if not np.isnan(r.get("eps2_min", np.nan)))
        print(f"  Completed: {excluded}/{len(valid_results)} mass points have sensitivity")

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
                p_min_GeV=p_min_GeV,
                reco_efficiency=reco_efficiency,
                show_progress=show_progress,
                timing_enabled=timing_enabled,
                hnlcalc_per_eps2=hnlcalc_per_eps2,
            )
            if res:
                results.append(res)

    return pd.DataFrame(results)

def scan_single_mass_wrapper(args):
    (mass_val, mass_str, flavour, benchmark, lumi_fb, sim_list, dirac, separation_m, decay_seed, p_min_GeV, reco_efficiency, show_progress, timing_enabled, hnlcalc_per_eps2) = args
    return scan_single_mass(
        mass_val, mass_str, flavour, benchmark, lumi_fb, sim_list,
        dirac=dirac, separation_m=separation_m, decay_seed=decay_seed,
        p_min_GeV=p_min_GeV, reco_efficiency=reco_efficiency, quiet=True,
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
        "--max-mass",
        type=float,
        default=None,
        help="Only process mass points up to this value in GeV (e.g., 5.0).",
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
        "--p-min-gev",
        type=float,
        default=0.5,
        help="Minimum charged-track momentum in GeV/c (default: 0.5, MATHUSLA CDR threshold).",
    )
    parser.add_argument(
        "--reco-efficiency",
        type=float,
        default=None,
        help="Flat reconstruction efficiency factor (recommended: 0.5 per MATHUSLA/ANUBIS). "
             "If omitted, defaults to 1.0 with a notice.",
    )
    parser.add_argument(
        "--hnlcalc-per-eps2",
        action="store_true",
        help="Recompute HNLCalc for every eps2 point (slow, legacy behavior).",
    )
    parser.add_argument(
        "--allow-variant-drop",
        action="store_true",
        help="Allow silently dropping lower-priority pTHat/QCD variants (default: error).",
    )
    args = parser.parse_args()
    if args.separation_mm <= 0:
        raise ValueError("--separation-mm must be positive.")
    if args.p_min_gev < 0:
        raise ValueError("--p-min-gev must be non-negative.")
    if args.reco_efficiency is None:
        args.reco_efficiency = 1.0
        print("[NOTICE] --reco-efficiency not set, defaulting to 1.0 (no efficiency loss).")
        print("         For realistic projections, use --reco-efficiency 0.5 (MATHUSLA/ANUBIS).")
    if not 0.0 < args.reco_efficiency <= 1.0:
        raise ValueError("--reco-efficiency must be in (0, 1].")

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
    p_min_GeV = args.p_min_gev
    reco_efficiency = args.reco_efficiency
    print(f"Decay separation: {args.separation_mm:.3f} mm (seed={args.decay_seed})")
    print(f"Track p_min: {p_min_GeV:.2f} GeV/c | Reco efficiency: {reco_efficiency:.2f}")

    all_results = []

    show_progress = None if not args.no_progress else False
    timing_enabled = args.timing

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
            p_min_GeV=p_min_GeV,
            reco_efficiency=reco_efficiency,
            show_progress=show_progress,
            mass_filter=args.mass,
            timing_enabled=timing_enabled,
            hnlcalc_per_eps2=args.hnlcalc_per_eps2,
            allow_variant_drop=args.allow_variant_drop,
            max_mass=args.max_mass,
        )
        df["separation_mm"] = args.separation_mm
        all_results.append(df)

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
