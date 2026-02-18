#!/usr/bin/env python

import sys
import re
import os
import io
import json
import tempfile
import contextlib
from datetime import datetime, timezone
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
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

MASS_FILTER_TOL = 5e-4

from decay.brvis_kappa import DEFAULT_KAPPA_TABLE, lookup_kappa, resolve_kappa_table_path
from geometry.per_parent_efficiency import (
    DEFAULT_DETECTOR_THICKNESS_M,
    GeometryConfig,
    build_drainage_gallery_mesh,
    geometry_metadata,
    geometry_tag,
    is_default_geometry_config,
    normalize_geometry_config,
    preprocess_hnl_csv,
)
from limits.expected_signal import expected_signal_events, couplings_from_eps2
from limits.overlap_resolution import (
    OverlapSample,
    filter_dataframe_by_norm_keys,
    resolve_parent_overlap,
)
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


def _save_geom_cache(geom_df: pd.DataFrame, geom_csv: Path) -> None:
    """Save only hit rows (row_idx, entry_distance, path_length) to cache."""
    hit_mask = geom_df["hits_tube"].values.astype(bool)
    hit_indices = np.flatnonzero(hit_mask)
    slim = pd.DataFrame({
        "row_idx": hit_indices,
        "entry_distance": geom_df["entry_distance"].values[hit_mask],
        "path_length": geom_df["path_length"].values[hit_mask],
    })
    with tempfile.NamedTemporaryFile(mode='w', dir=geom_csv.parent,
                                      suffix='.tmp', delete=False) as tmp:
        slim.to_csv(tmp.name, index=False)
    os.replace(tmp.name, geom_csv)


def _load_geom_cached(sim_csv: Path, geom_csv: Path) -> pd.DataFrame:
    """Load sim CSV and apply cached geometry columns."""
    cache_df = pd.read_csv(geom_csv)

    # Detect legacy full-copy cache (has 'hits_tube' column) vs slim cache (has 'row_idx')
    if "hits_tube" in cache_df.columns and "row_idx" not in cache_df.columns:
        return cache_df

    geom_df = pd.read_csv(sim_csv)
    n = len(geom_df)
    hits_tube = np.zeros(n, dtype=bool)
    entry_distance = np.full(n, np.nan, dtype=float)
    path_length = np.full(n, np.nan, dtype=float)

    if len(cache_df) > 0:
        idx = cache_df["row_idx"].values.astype(int)
        hits_tube[idx] = True
        entry_distance[idx] = cache_df["entry_distance"].values
        path_length[idx] = cache_df["path_length"].values

    geom_df["hits_tube"] = hits_tube
    geom_df["entry_distance"] = entry_distance
    geom_df["path_length"] = path_length

    # Reproduce derived columns that preprocess_hnl_csv would compute
    if "parent_pdg" in geom_df.columns and "parent_id" not in geom_df.columns:
        geom_df["parent_id"] = geom_df["parent_pdg"].abs()
    if "p" in geom_df.columns and "momentum" not in geom_df.columns:
        geom_df["momentum"] = geom_df["p"]
    if "beta_gamma" not in geom_df.columns:
        geom_df["beta_gamma"] = geom_df["momentum"] / geom_df["mass"]
    if "weight" not in geom_df.columns:
        geom_df["weight"] = 1.0

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
    max_separation_m=None,
    separation_policy="all-pairs-min",
    decay_seed=12345,
    p_min_GeV=0.6,
    geometry_config=None,
    reco_efficiency=1.0,
    quiet=False,
    show_progress=None,
    timing_enabled=False,
    hnlcalc_per_eps2=False,
    decay_mode="library",
    kappa_table_path=None,
):
    stdout_ctx = contextlib.redirect_stdout(io.StringIO()) if quiet else contextlib.nullcontext()

    with stdout_ctx:
        return _scan_single_mass_impl(
            mass_val, mass_str, flavour, benchmark, lumi_fb, sim_files,
            dirac, separation_m, max_separation_m, separation_policy,
            decay_seed, p_min_GeV, geometry_config, reco_efficiency,
            quiet, show_progress, timing_enabled, hnlcalc_per_eps2,
            decay_mode, kappa_table_path,
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
    max_separation_m,
    separation_policy,
    decay_seed,
    p_min_GeV,
    geometry_config,
    reco_efficiency,
    quiet,
    show_progress,
    timing_enabled,
    hnlcalc_per_eps2,
    decay_mode,
    kappa_table_path,
):
    timing = {} if timing_enabled else None
    if timing is not None:
        timing["count_geom_files"] = len(sim_files)

    geometry_cfg = normalize_geometry_config(geometry_config)
    geom_tag = geometry_tag(geometry_cfg)
    default_geometry = is_default_geometry_config(geometry_cfg)

    if not quiet:
        print(f"\n[{flavour} {mass_val} GeV] Processing ({len(sim_files)} production file(s))...")

    geom_dfs = []
    mesh = None

    for entry in sim_files:
        if len(entry) == 2:
            sim_csv, regime = entry
            owned_keys = ()
        elif len(entry) == 3:
            sim_csv, regime, owned_keys = entry
        else:
            raise ValueError(f"Unsupported sim file entry: {entry}")
        geom_cache_name = f"{sim_csv.stem}_geom_{geom_tag}.csv"
        geom_csv = GEOM_CACHE_DIR / geom_cache_name

        # Check legacy cache names only for default geometry config.
        if default_geometry and not geom_csv.exists():
            for legacy_stem in (
                sim_csv.stem + "_geom",
                f"HNL_{mass_str}GeV_{flavour}_combined_geom",
                f"HNL_{mass_str}GeV_{flavour}_geom",
            ):
                legacy_path = GEOM_CACHE_DIR / f"{legacy_stem}.csv"
                if legacy_path.exists():
                    geom_csv = legacy_path
                    break

        if geom_csv.exists() and sim_csv.stat().st_mtime > geom_csv.stat().st_mtime:
            if not quiet:
                print(f"  Stale cache detected, regenerating: {geom_csv.name}")
            geom_csv.unlink()

        if geom_csv.exists():
            _count(timing, "count_geom_cache_hits")
            with _time_block(timing, "time_geom_load_s"):
                geom_df = _load_geom_cached(sim_csv, geom_csv)
        else:
            _count(timing, "count_geom_cache_misses")
            if mesh is None:
                with _time_block(timing, "time_mesh_build_s"):
                    mesh = build_drainage_gallery_mesh(geometry_cfg)
            if not quiet:
                print(f"  Computing geometry for {sim_csv.name} (caching to {geom_csv.name})...")
            with _time_block(timing, "time_geom_compute_s"):
                geom_df = preprocess_hnl_csv(sim_csv, mesh, show_progress=show_progress)

            with _time_block(timing, "time_geom_write_s"):
                _save_geom_cache(geom_df, geom_csv)

        geom_df = _attach_sim_metadata(geom_df, sim_csv)

        if len(owned_keys) > 0:
            n_before = len(geom_df)
            geom_df = filter_dataframe_by_norm_keys(geom_df, owned_keys)
            if len(geom_df) == 0:
                if not quiet:
                    print(f"  Skipping {sim_csv.name}: overlap resolver kept 0 events")
                continue
            if not quiet and len(geom_df) != n_before:
                print(
                    f"  Overlap ownership filter: kept {len(geom_df)}/{n_before} rows from {sim_csv.name}"
                )

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

    decay_mode_norm = str(decay_mode).strip().lower().replace("-", "_")
    if decay_mode_norm not in {"library", "brvis_kappa"}:
        raise ValueError(
            f"Unsupported decay_mode='{decay_mode}'. Use 'library' or 'brvis-kappa'."
        )

    eps2_ref = 1e-6
    ctau0_ref = None
    br_ref = None
    br_vis = None
    need_ref_model = (not hnlcalc_per_eps2) or (decay_mode_norm == "brvis_kappa")
    if need_ref_model:
        from hnl_models.hnl_model_hnlcalc import HNLModel

        with _time_block(timing, "time_hnlcalc_ref_s"):
            Ue2, Umu2, Utau2 = couplings_from_eps2(eps2_ref, benchmark)
            model = HNLModel(mass_GeV=mass_val, Ue2=Ue2, Umu2=Umu2, Utau2=Utau2)
            ctau0_ref = model.ctau0_m
            br_ref = model.production_brs()
            br_vis = model.visible_branching_ratio()
        if timing is not None:
            timing["n_hnlcalc_eps2_ref"] = float(eps2_ref)
            timing["br_visible"] = float(br_vis)

    kappa_eff = np.nan
    kappa_path_resolved = ""
    if decay_mode_norm == "brvis_kappa":
        kappa_path = resolve_kappa_table_path(kappa_table_path)
        kappa_path_resolved = str(kappa_path)
        kappa_eff = lookup_kappa(
            flavour=flavour,
            mass_GeV=mass_val,
            p_min_GeV=p_min_GeV,
            separation_mm=separation_m * 1e3,
            table_path=kappa_path,
            geometry_config=geometry_cfg,
        )
        if not quiet:
            print(
                f"  BR_vis(>=2 charged)={float(br_vis):.4f}, "
                f"kappa={float(kappa_eff):.4f} [{kappa_path.name}]"
            )

    decay_cache = None
    if decay_mode_norm == "library":
        from decay.decay_detector import DecaySelection, build_decay_cache

        with _time_block(timing, "time_decay_cache_s"):
            decay_cache = build_decay_cache(
                geom_df,
                mass_val,
                flavour,
                DecaySelection(
                    separation_m=separation_m,
                    seed=decay_seed,
                    p_min_GeV=p_min_GeV,
                    max_separation_m=max_separation_m,
                    separation_policy=separation_policy,
                ),
                verbose=not quiet,
            )
        if not quiet:
            print("  Precomputing decay cache for separation scan...")

    with _time_block(timing, "time_eps2_scan_s"):
        for eps2 in eps2_scan:
            base_kwargs = dict(
                geom_df=geom_df,
                mass_GeV=mass_val,
                eps2=eps2,
                benchmark=benchmark,
                lumi_fb=lumi_fb,
                dirac=dirac,
                separation_m=separation_m,
                max_separation_m=max_separation_m,
                separation_policy=separation_policy,
                decay_seed=decay_seed,
                p_min_GeV=p_min_GeV,
                geometry_config=geometry_cfg,
                reco_efficiency=reco_efficiency,
                decay_mode=decay_mode_norm,
                br_vis=br_vis if decay_mode_norm == "brvis_kappa" else None,
                kappa_eff=float(kappa_eff) if decay_mode_norm == "brvis_kappa" else None,
                timing=timing,
            )
            if hnlcalc_per_eps2:
                N = expected_signal_events(
                    **base_kwargs,
                    decay_cache=decay_cache,
                )
            else:
                ctau0_m = ctau0_ref * (eps2_ref / eps2)
                br_scale = eps2 / eps2_ref
                N = expected_signal_events(
                    **base_kwargs,
                    decay_cache=decay_cache,
                    ctau0_m=ctau0_m,
                    br_per_parent=br_ref,
                    br_scale=br_scale,
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
            "peak_events": N_scan.max(),
            "decay_mode": decay_mode_norm,
            "br_vis": float(br_vis) if br_vis is not None else np.nan,
            "kappa_eff": float(kappa_eff) if np.isfinite(kappa_eff) else np.nan,
            "kappa_table_path": kappa_path_resolved,
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
        "peak_events": peak_events,
        "decay_mode": decay_mode_norm,
        "br_vis": float(br_vis) if br_vis is not None else np.nan,
        "kappa_eff": float(kappa_eff) if np.isfinite(kappa_eff) else np.nan,
        "kappa_table_path": kappa_path_resolved,
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
    max_separation_m=None,
    separation_policy="all-pairs-min",
    decay_seed=12345,
    p_min_GeV=0.6,
    geometry_config: GeometryConfig | None = None,
    reco_efficiency=1.0,
    show_progress=None,
    mass_filter=None,
    timing_enabled=False,
    hnlcalc_per_eps2=False,
    allow_variant_drop=False,
    overlap_min_events=0,
    strict_overlap_min_events=False,
    max_mass=None,
    decay_mode="library",
    kappa_table_path=None,
    allow_legacy_tau_all=False,
):
    print(f"\n{'='*60}")
    print(f"FLAVOUR: {flavour.upper()} (Benchmark {benchmark})")
    print(f"{'='*60}")

    geometry_cfg = normalize_geometry_config(geometry_config)

    pattern = re.compile(
        rf"^HNL_([0-9]+p[0-9]{{1,2}})GeV_{flavour}_"
        r"((?:kaon|charm|beauty|Bc|ew|all|combined)(?:_ff)?)"
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

        if qcd_mode in {"hardccbar", "hardbbbar"}:
            print(f"[WARN] Skipping hard-sliced sample in nominal analysis: {f.name}")
            continue

        is_ff = regime_token.endswith("_ff")
        base_regime = regime_token.replace("_ff", "")

        files.append((mass_val, mass_str, base_regime, mode, is_ff, qcd_mode, pthat_min, f))

    files_by_mass = {}
    for mass_val, mass_str, base_regime, mode, is_ff, qcd_mode, pthat_min, path in files:
        key = (mass_val, mass_str)
        files_by_mass.setdefault(key, []).append((base_regime, mode, is_ff, qcd_mode, pthat_min, path))

    if flavour == "tau" and not allow_legacy_tau_all:
        filtered_by_mass = {}
        legacy_only_masses = []
        for key, items in files_by_mass.items():
            legacy_items = [it for it in items if it[0] in {"all", "combined"}]
            component_items = [it for it in items if it[0] not in {"all", "combined"}]
            if component_items:
                if legacy_items:
                    dropped = ", ".join(sorted({it[5].name for it in legacy_items}))
                    print(
                        f"[WARN] m={key[0]:.2f} (tau): ignoring legacy tau_all/tau_combined inputs: {dropped}"
                    )
                filtered_by_mass[key] = component_items
                continue
            if legacy_items:
                legacy_only_masses.append((key[0], sorted(it[5].name for it in legacy_items)))

        if legacy_only_masses:
            legacy_only_masses = sorted(legacy_only_masses, key=lambda x: x[0])
            preview = ", ".join(f"{m:.2f}" for m, _ in legacy_only_masses[:10])
            more = "..." if len(legacy_only_masses) > 10 else ""
            raise ValueError(
                "Tau analysis requires explicit component files (direct/fromTau/ew); "
                "legacy *_tau_all.csv or *_tau_combined.csv are not accepted by default. "
                f"Missing component masses: {preview}{more}. "
                "Regenerate tau production and rerun, or pass --allow-legacy-tau-all to bypass."
            )
        files_by_mass = filtered_by_mass

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
        regime_order = {"kaon": 0, "charm": 1, "beauty": 2, "Bc": 3, "ew": 4, "all": 5, "combined": 5}
        mode_order = {None: 0, "direct": 1, "fromTau": 2}
        variant_order = _variant_priority(base_regime, is_ff, qcd_mode, pthat_min)
        return (regime_order.get(base_regime, 99), mode_order.get(mode, 99), base_regime, mode or "", variant_order)

    def _resolve_overlap_selection(
        mass_key: tuple[float, str],
        variants: list[tuple[str, str | None, bool, str, float | None, Path]],
        *,
        warn_for_key: bool,
    ) -> list[tuple[Path, str, tuple[tuple[str, int], ...]]]:
        overlap_samples = [
            OverlapSample(
                base_regime=base_regime,
                mode=mode,
                is_ff=is_ff,
                qcd_mode=qcd_mode,
                pthat_min=pthat_min,
                path=path,
            )
            for base_regime, mode, is_ff, qcd_mode, pthat_min, path in variants
        ]
        context = f"m={mass_key[0]:.2f} ({flavour})"
        overlap = resolve_parent_overlap(
            overlap_samples,
            context=context,
            min_events_per_mass=int(overlap_min_events),
            strict_min_events=bool(strict_overlap_min_events),
        )
        if warn_for_key:
            for msg in overlap.warnings:
                print(msg)

        resolved_entries: list[tuple[Path, str, tuple[tuple[str, int], ...]]] = []
        for resolved in overlap.samples:
            sample = resolved.sample
            resolved_entries.append(
                (
                    sample.path,
                    _label(sample.base_regime, sample.mode, sample.is_ff, sample.qcd_mode, sample.pthat_min),
                    resolved.owned_keys,
                )
            )
        return resolved_entries

    selected_by_mass = {}
    for key, items in files_by_mass.items():
        warn_for_key = True
        if mass_filter is not None:
            warn_for_key = abs(key[0] - mass_filter) <= MASS_FILTER_TOL

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
        selected_by_mass[key] = _resolve_overlap_selection(
            key,
            selected,
            warn_for_key=warn_for_key,
        )

    files_by_mass = selected_by_mass

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
                (
                    mass_val,
                    mass_str,
                    flavour,
                    benchmark,
                    lumi_fb,
                    sim_list,
                    dirac,
                    separation_m,
                    max_separation_m,
                    separation_policy,
                    decay_seed,
                    p_min_GeV,
                    geometry_cfg,
                    reco_efficiency,
                    show_progress,
                    timing_enabled,
                    hnlcalc_per_eps2,
                    decay_mode,
                    kappa_table_path,
                )
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
                max_separation_m=max_separation_m,
                separation_policy=separation_policy,
                decay_seed=decay_seed,
                p_min_GeV=p_min_GeV,
                geometry_config=geometry_cfg,
                reco_efficiency=reco_efficiency,
                show_progress=show_progress,
                timing_enabled=timing_enabled,
                hnlcalc_per_eps2=hnlcalc_per_eps2,
                decay_mode=decay_mode,
                kappa_table_path=kappa_table_path,
            )
            if res:
                results.append(res)

    return pd.DataFrame(results)

def scan_single_mass_wrapper(args):
    (
        mass_val,
        mass_str,
        flavour,
        benchmark,
        lumi_fb,
        sim_list,
        dirac,
        separation_m,
        max_separation_m,
        separation_policy,
        decay_seed,
        p_min_GeV,
        geometry_config,
        reco_efficiency,
        show_progress,
        timing_enabled,
        hnlcalc_per_eps2,
        decay_mode,
        kappa_table_path,
    ) = args
    return scan_single_mass(
        mass_val, mass_str, flavour, benchmark, lumi_fb, sim_list,
        dirac=dirac, separation_m=separation_m, max_separation_m=max_separation_m,
        separation_policy=separation_policy, decay_seed=decay_seed,
        p_min_GeV=p_min_GeV, geometry_config=geometry_config,
        reco_efficiency=reco_efficiency, quiet=True,
        show_progress=show_progress, timing_enabled=timing_enabled,
        hnlcalc_per_eps2=hnlcalc_per_eps2,
        decay_mode=decay_mode,
        kappa_table_path=kappa_table_path,
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
        "--max-separation-mm",
        type=float,
        default=None,
        help="Optional exploratory upper bound on charged-track separation at detector surface in mm.",
    )
    parser.add_argument(
        "--separation-policy",
        type=str,
        choices=["all-pairs-min", "any-pair-window"],
        default="all-pairs-min",
        help=(
            "Pairwise separation policy: all-pairs-min (baseline) requires all pair distances in window; "
            "any-pair-window passes if any pair is in window."
        ),
    )
    parser.add_argument(
        "--geometry-model",
        type=str,
        choices=["tube", "profile"],
        default="tube",
        help="Geometry model: tube (baseline) or exploratory profile.",
    )
    parser.add_argument(
        "--detector-thickness-m",
        type=float,
        default=DEFAULT_DETECTOR_THICKNESS_M,
        help=(
            "Detector thickness in meters for geometry-model=profile (default: "
            f"{DEFAULT_DETECTOR_THICKNESS_M})."
        ),
    )
    parser.add_argument(
        "--profile-inset-floor",
        action="store_true",
        help="For geometry-model=profile, inset the floor by detector thickness (exploratory).",
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
        default=0.6,
        help="Minimum charged-track momentum in GeV/c (default: 0.6).",
    )
    parser.add_argument(
        "--decay-mode",
        type=str,
        choices=["library", "brvis-kappa"],
        default="library",
        help="Decay acceptance mode: full library sampling or calibrated BR_vis*kappa surrogate.",
    )
    parser.add_argument(
        "--kappa-table",
        type=str,
        default=str(DEFAULT_KAPPA_TABLE),
        help="Path to calibrated kappa table (used by --decay-mode brvis-kappa).",
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
    parser.add_argument(
        "--overlap-min-events",
        type=int,
        default=0,
        help=(
            "Optional minimum owned simulation events required per mass point after overlap resolution. "
            "Set 0 to disable."
        ),
    )
    parser.add_argument(
        "--strict-overlap-min-events",
        action="store_true",
        help="Fail when --overlap-min-events is not met (default: warn only).",
    )
    parser.add_argument(
        "--allow-legacy-tau-all",
        action="store_true",
        help=(
            "Allow legacy tau _all/_combined files. "
            "Default behavior requires explicit tau components (direct/fromTau/ew)."
        ),
    )
    args = parser.parse_args()
    if args.separation_mm <= 0:
        raise ValueError("--separation-mm must be positive.")
    if args.max_separation_mm is not None and args.max_separation_mm <= 0:
        raise ValueError("--max-separation-mm must be positive when provided.")
    if args.max_separation_mm is not None and args.max_separation_mm <= args.separation_mm:
        raise ValueError("--max-separation-mm must be strictly greater than --separation-mm.")
    if args.p_min_gev < 0:
        raise ValueError("--p-min-gev must be non-negative.")
    if args.overlap_min_events < 0:
        raise ValueError("--overlap-min-events must be non-negative.")
    if args.reco_efficiency is None:
        args.reco_efficiency = 1.0
        print("[NOTICE] --reco-efficiency not set, defaulting to 1.0 (no efficiency loss).")
        print("         For realistic projections, use --reco-efficiency 0.5 (MATHUSLA/ANUBIS).")
    if not 0.0 < args.reco_efficiency <= 1.0:
        raise ValueError("--reco-efficiency must be in (0, 1].")

    geometry_cfg = normalize_geometry_config(
        GeometryConfig(
            model=args.geometry_model,
            detector_thickness_m=args.detector_thickness_m,
            profile_inset_floor=args.profile_inset_floor,
        )
    )

    if args.decay_mode == "brvis-kappa" and args.max_separation_mm is not None:
        raise ValueError(
            "--max-separation-mm is not supported with --decay-mode brvis-kappa. "
            "Treat max separation as exploratory library-only cut."
        )
    if args.decay_mode == "brvis-kappa" and args.separation_policy != "all-pairs-min":
        raise ValueError(
            "--separation-policy must be 'all-pairs-min' with --decay-mode brvis-kappa."
        )

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
    max_separation_m = args.max_separation_mm * 1e-3 if args.max_separation_mm is not None else None
    p_min_GeV = args.p_min_gev
    reco_efficiency = args.reco_efficiency
    print(f"Decay separation: {args.separation_mm:.3f} mm (seed={args.decay_seed})")
    if args.max_separation_mm is not None:
        print(f"Max separation (exploratory): {args.max_separation_mm:.3f} mm")
    print(f"Separation policy: {args.separation_policy}")
    print(f"Track p_min: {p_min_GeV:.2f} GeV/c | Reco efficiency: {reco_efficiency:.2f}")
    if args.overlap_min_events > 0:
        mode = "strict" if args.strict_overlap_min_events else "warn"
        print(f"Overlap minimum events: {args.overlap_min_events} ({mode})")
    if args.allow_legacy_tau_all:
        print("[WARN] Legacy tau _all/_combined inputs are enabled.")
    print(f"Decay mode: {args.decay_mode}")
    print(f"Geometry: {geometry_cfg.model} (tag={geometry_tag(geometry_cfg)})")
    if args.decay_mode == "brvis-kappa":
        print(f"Kappa table: {resolve_kappa_table_path(args.kappa_table)}")

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
            max_separation_m=max_separation_m,
            separation_policy=args.separation_policy,
            decay_seed=args.decay_seed,
            p_min_GeV=p_min_GeV,
            geometry_config=geometry_cfg,
            reco_efficiency=reco_efficiency,
            show_progress=show_progress,
            mass_filter=args.mass,
            timing_enabled=timing_enabled,
            hnlcalc_per_eps2=args.hnlcalc_per_eps2,
            allow_variant_drop=args.allow_variant_drop,
            overlap_min_events=args.overlap_min_events,
            strict_overlap_min_events=args.strict_overlap_min_events,
            max_mass=args.max_mass,
            decay_mode=args.decay_mode,
            kappa_table_path=args.kappa_table,
            allow_legacy_tau_all=args.allow_legacy_tau_all,
        )
        df["separation_mm"] = args.separation_mm
        df["p_min_gev"] = p_min_GeV
        df["decay_mode"] = args.decay_mode
        df["kappa_table_path"] = (
            str(resolve_kappa_table_path(args.kappa_table))
            if args.decay_mode == "brvis-kappa"
            else ""
        )
        all_results.append(df)

    final_df = pd.concat(all_results, ignore_index=True)
    final_df.to_csv(results_out, index=False)

    summary_meta = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "separation_mm": float(args.separation_mm),
        "max_separation_mm": None if args.max_separation_mm is None else float(args.max_separation_mm),
        "separation_policy": args.separation_policy,
        "p_min_gev": float(p_min_GeV),
        "decay_mode": args.decay_mode,
        "overlap_min_events": int(args.overlap_min_events),
        "strict_overlap_min_events": bool(args.strict_overlap_min_events),
        "allow_legacy_tau_all": bool(args.allow_legacy_tau_all),
        "kappa_table_path": (
            str(resolve_kappa_table_path(args.kappa_table)) if args.decay_mode == "brvis-kappa" else ""
        ),
        "geometry": geometry_metadata(geometry_cfg),
    }
    summary_meta_path = ANALYSIS_OUT_DIR / "HNL_U2_limits_summary.meta.json"
    with summary_meta_path.open("w", encoding="utf-8") as f:
        json.dump(summary_meta, f, indent=2, sort_keys=True)

    if timing_enabled:
        timing_out = Path(args.timing_out) if args.timing_out else (ANALYSIS_OUT_DIR / "HNL_U2_timing.csv")
        timing_cols = [c for c in final_df.columns if c.startswith("time_") or c.startswith("count_") or c.startswith("n_")]
        timing_df = final_df[["mass_GeV", "flavour", "benchmark"] + timing_cols]
        timing_df.to_csv(timing_out, index=False)

    print(f"\n{'='*60}")
    print(f"COMPLETE!")
    print(f"Saved {len(final_df)} mass points to:")
    print(f"  {results_out}")
    print("Run metadata:")
    print(f"  {summary_meta_path}")
    if timing_enabled:
        print(f"Timing breakdown saved to:")
        print(f"  {timing_out}")
    print(f"{'='*60}")
