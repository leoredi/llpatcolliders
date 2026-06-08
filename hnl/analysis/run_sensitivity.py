#!/usr/bin/env python3
"""
hnl/analysis/run_sensitivity.py

Main driver: compute HNL sensitivity from the combined 4-vector CSVs that
``hnl/production/`` writes, and produce the (m_N, U^2) money plot.

Default flavor is Umu only (HNL pattern 010 -- the muon-mixing scenario).

Usage examples:

    python -m analysis.run_sensitivity                       # Umu, full grid
    python -m analysis.run_sensitivity --flavor Ue Umu       # multi-flavor
    python -m analysis.run_sensitivity --mass 1.0 2.0        # custom masses
    python -m analysis.run_sensitivity --workers 3           # parallel workers
    python -m analysis.run_sensitivity --plot-only           # re-plot only

Inputs (read via ``production.paths``):

    {LLP_VECTORS_DIR}/{flavor}/combined/mN_{mass}.csv    -- 4-vectors per point
    {HNL_ROOT}/data/ctau/ctau_{flavor}.dat               -- ctau(m_N) at U^2=1
    {HNL_ROOT}/data/ctau/br_vis_{flavor}.dat             -- BR_vis(m_N), optional
    {HNL_ROOT}/geometry/grendel_geometry.py              -- detector mesh

Outputs:

    {ANALYSIS_DIR}/hnl_sensitivity.csv
    {ANALYSIS_DIR}/hnl_exclusion.{pdf,png}
    {ANALYSIS_DIR}/geometry_cache/{flavor}/geom_{mass}.npz
    {ANALYSIS_DIR}/diagnostics/                          -- with --diagnostics
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

HNL_ROOT = Path(__file__).resolve().parent.parent
if str(HNL_ROOT) not in sys.path:
    sys.path.insert(0, str(HNL_ROOT))

from config_mass_grid import MASS_GRID, format_mass_for_filename
from production.paths import ANALYSIS_DIR, LLP_VECTORS_DIR

from analysis.constants import (
    CMS_ORIGIN,
    DEFAULT_FLAVORS,
    FLAVORS,
    ANALYSIS_MASS_MAX,
    L_INT_PB,
    LOG_U2_MAX,
    LOG_U2_MIN,
    M_ELECTRON,
    N_THRESHOLD,
    N_U2_POINTS,
)
from analysis.exclusion import find_exclusion_band
from analysis.format_bridge import load_combined_csv
from analysis.plot_exclusion import plot_exclusion, plot_nsignal_vs_u2
from analysis.sensitivity import scan_u2

CTAU_DIR = HNL_ROOT / "data" / "ctau"
GEOMETRY_DIR = HNL_ROOT / "geometry"
GEOM_CACHE_DIR = ANALYSIS_DIR / "geometry_cache"


# =========================================================================
# Geometry: ray-cast through the GRENDEL mesh (lazy import)
# =========================================================================

def _get_mesh():
    """Lazy-import the GRENDEL mesh from ``hnl/geometry/``."""
    if not GEOMETRY_DIR.exists():
        raise FileNotFoundError(
            f"geometry/ directory not found at {GEOMETRY_DIR}.\n"
            "The checkout is incomplete; hnl/geometry/grendel_geometry.py is required."
        )
    if str(GEOMETRY_DIR) not in sys.path:
        sys.path.insert(0, str(GEOMETRY_DIR))
    from grendel_geometry import mesh_fiducial
    return mesh_fiducial


def _eta_phi_to_directions_batch(eta, phi):
    """(eta, phi) arrays -> (N, 3) unit direction vectors."""
    theta = 2.0 * np.arctan(np.exp(-eta))
    dx = np.sin(theta) * np.cos(phi)
    dy = np.sin(theta) * np.sin(phi)
    dz = np.cos(theta)
    return np.column_stack([dx, dy, dz])


def compute_geometry(eta, phi, mesh, origin=CMS_ORIGIN, batch_label=""):
    """Batch ray-cast (eta, phi) directions against the GRENDEL mesh."""
    n = len(eta)
    origin_arr = np.array(origin, dtype=np.float64)
    hits = np.zeros(n, dtype=bool)
    entry_d = np.full(n, np.nan)
    exit_d = np.full(n, np.nan)

    directions = _eta_phi_to_directions_batch(eta, phi)
    candidates = np.where(directions[:, 1] > 0.01)[0]
    n_cand = len(candidates)
    if n_cand == 0:
        return hits, entry_d, exit_d

    print(f"  Batch ray-casting {n_cand}/{n} candidates {batch_label}...",
          flush=True)

    cand_dirs = directions[candidates]
    origins = np.tile(origin_arr, (n_cand, 1))
    locations, ray_ids, _ = mesh.ray.intersects_location(
        ray_origins=origins, ray_directions=cand_dirs)

    if len(locations) == 0:
        print(f"  0/{n} events hit detector", flush=True)
        return hits, entry_d, exit_d

    dists = np.linalg.norm(locations - origin_arr, axis=1)
    order = np.argsort(ray_ids)
    sorted_ray_ids = ray_ids[order]
    sorted_dists = dists[order]

    unique_rays, start_idx, counts = np.unique(
        sorted_ray_ids, return_index=True, return_counts=True)

    valid = counts >= 2
    valid_rays = unique_rays[valid]
    valid_starts = start_idx[valid]
    valid_counts = counts[valid]

    for i in range(len(valid_rays)):
        ray_local = valid_rays[i]
        orig_idx = candidates[ray_local]
        s = valid_starts[i]
        e = s + valid_counts[i]
        ray_dists = sorted_dists[s:e]
        ray_dists.sort()
        hits[orig_idx] = True
        entry_d[orig_idx] = ray_dists[0]
        exit_d[orig_idx] = ray_dists[1]

    n_hits = int(hits.sum())
    print(f"  {n_hits}/{n} events hit detector ({n_hits / n * 100:.2f}%)",
          flush=True)
    if n_hits > 0:
        path_lens = exit_d[hits] - entry_d[hits]
        print(f"  Mean path length: {path_lens.mean():.2f} m", flush=True)
    return hits, entry_d, exit_d


def load_or_compute_geometry(flavor, mass_label, eta, phi, mesh, force=False):
    """Load cached geometry or compute + cache as NPZ."""
    cache_dir = GEOM_CACHE_DIR / flavor
    cache_path = cache_dir / f"geom_{mass_label}.npz"

    if cache_path.exists() and not force:
        data = np.load(cache_path)
        return data["hits"].astype(bool), data["entry_d"], data["exit_d"]

    hits, entry_d, exit_d = compute_geometry(
        eta, phi, mesh, CMS_ORIGIN, batch_label=f"[{flavor}/{mass_label}]")

    cache_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(cache_path, hits=hits, entry_d=entry_d, exit_d=exit_d)
    return hits, entry_d, exit_d


# =========================================================================
# c tau / BR_vis tables
# =========================================================================

def _load_two_column_table(path):
    table = {}
    for line in path.read_text().strip().split("\n"):
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        table[float(parts[0])] = float(parts[1])
    return table


def load_ctau_table(flavor):
    """ctau(m_N) at U^2=1 in meters."""
    path = CTAU_DIR / f"ctau_{flavor}.dat"
    if not path.exists():
        raise FileNotFoundError(
            f"ctau table not found at {path}. "
            "The checkout is incomplete; hnl/data/ctau is required."
        )
    return _load_two_column_table(path)


def load_br_vis_table(flavor):
    """BR_vis(m_N) (visible decay fraction). None -> use 1.0 downstream."""
    path = CTAU_DIR / f"br_vis_{flavor}.dat"
    if not path.exists():
        print(f"  NOTE: BR_vis table not found at {path}; using BR_vis = 1.0")
        return None
    return _load_two_column_table(path)


def _lookup(table, mass):
    if table is None:
        return 1.0
    if mass in table:
        return table[mass]
    masses = np.array(list(table.keys()))
    idx = int(np.argmin(np.abs(masses - mass)))
    return table[masses[idx]]


# =========================================================================
# Single-point processing
# =========================================================================

def process_mass_point(flavor, mass, ctau_table, mesh,
                       br_vis_table=None,
                       force_geometry=False, save_diagnostic=False):
    """Process one (flavor, mass) point and return its exclusion band."""
    mass_label = format_mass_for_filename(mass)
    csv_path = LLP_VECTORS_DIR / flavor / "combined" / f"mN_{mass_label}.csv"

    if not csv_path.exists() or csv_path.stat().st_size == 0:
        return None

    data = load_combined_csv(csv_path, mass)
    n_events = len(data["weight"])
    if n_events == 0:
        return None

    hits, entry_d, exit_d = load_or_compute_geometry(
        flavor, mass_label, data["eta"], data["phi"], mesh,
        force=force_geometry)

    n_hits = int(hits.sum())
    if n_hits == 0:
        return {
            "mass_GeV": mass, "flavor": flavor,
            "u2_min": np.nan, "u2_max": np.nan,
            "u2_min_open": False, "u2_max_open": False,
            "peak_N": 0.0, "peak_u2": np.nan,
            "has_sensitivity": False, "n_events": n_events, "n_hits": 0,
        }

    ctau_u2_1 = _lookup(ctau_table, mass)
    if ctau_u2_1 <= 0:
        return None

    br_vis = _lookup(br_vis_table, mass)

    u2_grid, N_grid = scan_u2(
        data["weight"], data["beta_gamma"], data["gamma"], data["beta"],
        hits, entry_d, exit_d,
        ctau_u2_1, mass, L_INT_PB,
        m_daughter=M_ELECTRON, use_acceptance=True,
        log_u2_min=LOG_U2_MIN, log_u2_max=LOG_U2_MAX,
        n_points=N_U2_POINTS,
        br_vis=br_vis)

    result = find_exclusion_band(u2_grid, N_grid, N_THRESHOLD)
    result["mass_GeV"] = mass
    result["flavor"] = flavor
    result["n_events"] = n_events
    result["n_hits"] = n_hits

    if save_diagnostic:
        diag_dir = ANALYSIS_DIR / "diagnostics"
        plot_nsignal_vs_u2(u2_grid, N_grid, mass, flavor, diag_dir)

    return result


# =========================================================================
# Parallel worker plumbing
# =========================================================================

_WORKER_MESH = None


def _worker_init():
    global _WORKER_MESH
    _WORKER_MESH = _get_mesh()


def _worker_process_point(args):
    flavor, mass, ctau_table, br_vis_table, force_geom, save_diag = args
    t0 = time.time()
    result = process_mass_point(
        flavor, mass, ctau_table, _WORKER_MESH,
        br_vis_table=br_vis_table,
        force_geometry=force_geom, save_diagnostic=save_diag)
    elapsed = time.time() - t0
    mass_label = format_mass_for_filename(mass)
    tag = f"{flavor}/mN_{mass_label}"
    if result and result.get("has_sensitivity"):
        print(f"  {tag}: sensitivity found ({elapsed:.1f}s)", flush=True)
    elif result:
        print(f"  {tag}: done, no sensitivity ({elapsed:.1f}s)", flush=True)
    else:
        print(f"  {tag}: skipped ({elapsed:.1f}s)", flush=True)
    return result


# =========================================================================
# Driver
# =========================================================================

def run(flavors, masses, n_workers=1, force_geometry=False,
        save_diagnostics=False):
    ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)
    GEOM_CACHE_DIR.mkdir(parents=True, exist_ok=True)

    work_items = []
    for flavor in flavors:
        ctau = load_ctau_table(flavor)
        br_vis_tab = load_br_vis_table(flavor)
        for mass in masses:
            work_items.append(
                (flavor, mass, ctau, br_vis_tab,
                 force_geometry, save_diagnostics))

    n_total = len(work_items)
    print(f"\nProcessing {n_total} mass points "
          f"with {n_workers} worker(s)...\n", flush=True)

    t_start = time.time()
    results = []
    status_path = ANALYSIS_DIR / "scan_status.json"

    def _write_status(done, n_sens):
        elapsed = time.time() - t_start
        eta = (elapsed / max(done, 1)) * (n_total - done) if done > 0 else 0
        status = {
            "ts": datetime.now().isoformat(),
            "done": done, "total": n_total,
            "n_sensitive": n_sens,
            "elapsed_s": round(elapsed, 1),
            "eta_s": round(eta, 1),
        }
        status_path.write_text(json.dumps(status, indent=2))

    if n_workers <= 1:
        mesh = _get_mesh()
        for i, item in enumerate(work_items):
            flavor, mass, ctau, br_vis_tab, fg, sd = item
            t0 = time.time()
            r = process_mass_point(flavor, mass, ctau, mesh,
                                   br_vis_table=br_vis_tab,
                                   force_geometry=fg, save_diagnostic=sd)
            elapsed = time.time() - t0
            mass_label = format_mass_for_filename(mass)
            if r is not None:
                results.append(r)
            n_sens = sum(1 for r in results if r.get("has_sensitivity"))
            _write_status(i + 1, n_sens)
            print(f"  [{i+1}/{n_total}] {flavor}/mN_{mass_label} "
                  f"({elapsed:.1f}s)", flush=True)
    else:
        with ProcessPoolExecutor(
            max_workers=n_workers, initializer=_worker_init
        ) as pool:
            futures = {
                pool.submit(_worker_process_point, item): i
                for i, item in enumerate(work_items)
            }
            done_count = 0
            for future in as_completed(futures):
                done_count += 1
                r = future.result()
                if r is not None:
                    results.append(r)
                n_sens = sum(1 for r in results if r.get("has_sensitivity"))
                _write_status(done_count, n_sens)
                if done_count % 10 == 0 or done_count == n_total:
                    elapsed = time.time() - t_start
                    print(f"  [{done_count}/{n_total}] "
                          f"{n_sens} sensitive, "
                          f"{elapsed:.0f}s elapsed", flush=True)

    total_time = time.time() - t_start

    if not results:
        print("No results produced.")
        return

    results.sort(key=lambda r: (r["flavor"], r["mass_GeV"]))
    df = pd.DataFrame(results)
    out_csv = ANALYSIS_DIR / "hnl_sensitivity.csv"
    df.to_csv(out_csv, index=False)
    print(f"\nResults saved: {out_csv}")

    n_sens = int(df["has_sensitivity"].sum())
    print(f"Sensitivity found at {n_sens}/{len(df)} mass points")
    print(f"Total time: {total_time:.0f}s ({total_time/60:.1f} min)")

    meta = {
        "timestamp": datetime.now().isoformat(),
        "flavors": flavors,
        "n_masses": len(masses),
        "n_workers": n_workers,
        "n_results": len(results),
        "n_sensitive": n_sens,
        "mass_range": [float(min(masses)), float(max(masses))],
        "total_time_s": round(total_time, 1),
        "llp_vectors_dir": str(LLP_VECTORS_DIR),
    }
    (ANALYSIS_DIR / "run_metadata.json").write_text(json.dumps(meta, indent=2))

    plot_exclusion(out_csv, ANALYSIS_DIR)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="HNL sensitivity + (m_N, U^2) money plot")
    parser.add_argument(
        "--flavor", nargs="+", default=None, choices=FLAVORS,
        help=f"Flavors to process (default: {DEFAULT_FLAVORS})")
    parser.add_argument(
        "--mass", nargs="+", type=float, default=None,
        help=f"Specific masses in GeV (default: full grid <= {ANALYSIS_MASS_MAX} GeV)")
    parser.add_argument(
        "--plot-only", action="store_true",
        help="Re-plot from existing results CSV")
    parser.add_argument(
        "--force-geometry", action="store_true",
        help="Recompute geometry cache (ignore existing NPZ)")
    parser.add_argument(
        "--diagnostics", action="store_true",
        help="Save N_signal vs U^2 plots for each mass point")
    parser.add_argument(
        "--workers", type=int, default=3,
        help="Number of parallel workers (default: 3)")
    args = parser.parse_args(argv)

    if args.plot_only:
        csv_path = ANALYSIS_DIR / "hnl_sensitivity.csv"
        if not csv_path.exists():
            print(f"No results file found: {csv_path}")
            return 1
        plot_exclusion(csv_path, ANALYSIS_DIR)
        return 0

    flavors = args.flavor or DEFAULT_FLAVORS
    masses = args.mass or [m for m in MASS_GRID if m <= ANALYSIS_MASS_MAX]

    print("HNL Sensitivity Analysis")
    print(f"  Flavors: {flavors}")
    print(f"  Masses: {len(masses)} points "
          f"({min(masses):.2f} - {max(masses):.2f} GeV)")
    print(f"  Workers: {args.workers}")
    print(f"  Luminosity: {L_INT_PB:.0f} pb^-1 ({L_INT_PB/1e3:.0f} fb^-1)")
    print(f"  Threshold: N_signal >= {N_THRESHOLD}")
    print(f"  Inputs:    {LLP_VECTORS_DIR}")
    print(f"  Outputs:   {ANALYSIS_DIR}")

    run(flavors, masses,
        n_workers=args.workers,
        force_geometry=args.force_geometry,
        save_diagnostics=args.diagnostics)
    return 0


if __name__ == "__main__":
    sys.exit(main())
