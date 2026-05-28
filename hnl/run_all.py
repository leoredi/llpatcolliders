#!/usr/bin/env python3
"""
hnl/run_all.py

Parallel driver for the complete HNL production pipeline.

Submits one process per (flavor, channel) job:

  (Ue, Bmeson)  (Ue, Dmeson)  (Ue, Bc)  (Ue, tau)
  (Umu, ...)    ...
  (Utau, ...)   ...

= 12 independent jobs. Each one builds its own meson/tau pool (small cost
relative to BR computation) and processes all 116 mass points for that
(flavor, channel). At the end the combine step runs serially.

Usage:
    python -m hnl.run_all                       # full grid, all cores
    python -m hnl.run_all --workers 4           # cap concurrency
    python -m hnl.run_all --n-pool 50000        # smaller pools (faster)
    python -m hnl.run_all --masses 0.5 1.0 2.0  # subset of masses

Run from the hnl/ directory:
    cd hnl && python run_all.py
or as a module from the repo root with PYTHONPATH=hnl.
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

HNL_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(HNL_ROOT))
sys.path.insert(0, str(HNL_ROOT / "vendored" / "HNLCalc"))

from config_mass_grid import MASS_GRID
from production.constants import M_TAU


# Module-level worker functions (must be picklable for ProcessPoolExecutor).

def _worker_meson(flavor: str, channel: str, masses: list, n_pool: int, seed: int) -> str:
    """Process one (flavor, meson_channel) for all masses. channel ∈ {bottom, charm, bc}."""
    from production.decay_engine.generate_meson_csvs import (
        generate_pool, generate_bc_pool, process_channel,
    )
    rng = np.random.default_rng(seed)
    if channel == "bc":
        pool = generate_bc_pool(n_pool, rng)
        sigma = 0.0
    else:
        pool, sigma = generate_pool(channel, n_pool, rng)
    process_channel(flavor, channel, pool, sigma, masses, rng)
    return f"meson {flavor:>5s} {channel}"


def _worker_tau(flavor: str, masses: list, n_pool: int, seed: int) -> str:
    """Process induced-tau channel for one flavor across all (sub-mtau) masses."""
    from production.decay_engine.generate_induced_tau import build_tau_pool, process_flavor
    rng = np.random.default_rng(seed)
    tau_4v, tau_w = build_tau_pool(n_pool, rng)
    sub_masses = [m for m in masses if m < M_TAU]
    process_flavor(flavor, tau_4v, tau_w, sub_masses, rng)
    return f"tau   {flavor:>5s}"


def main():
    ap = argparse.ArgumentParser(description="Parallel HNL production driver")
    ap.add_argument("--flavor", choices=["Ue", "Umu", "Utau"], nargs="+",
                    default=["Ue", "Umu", "Utau"])
    ap.add_argument("--n-pool", type=int, default=100_000,
                    help="Sampled mesons per (flavor, channel) job (default: 100k)")
    ap.add_argument("--workers", type=int, default=None,
                    help="Max concurrent processes (default: os.cpu_count())")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--masses", type=float, nargs="+", default=None,
                    help="Custom mass list (default: full grid)")
    ap.add_argument("--skip-combine", action="store_true",
                    help="Skip the final combine_channels step")
    args = ap.parse_args()

    workers = args.workers or os.cpu_count() or 1
    masses = args.masses if args.masses else MASS_GRID
    meson_channels = ["bottom", "charm", "bc"]

    jobs = []   # list of (fn, args_tuple, label) used by submit loop
    seed = args.seed
    for flavor in args.flavor:
        for ch in meson_channels:
            jobs.append((_worker_meson, (flavor, ch, masses, args.n_pool, seed)))
            seed += 1
        jobs.append((_worker_tau, (flavor, masses, args.n_pool, seed)))
        seed += 1

    print(f"Submitting {len(jobs)} jobs across {workers} workers "
          f"(flavors={args.flavor}, masses={len(masses)}, n_pool={args.n_pool})")
    t0 = time.time()
    failures = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futures = {ex.submit(fn, *a): a for fn, a in jobs}
        for f in as_completed(futures):
            try:
                label = f.result()
                print(f"  [{time.time()-t0:6.1f}s]  done: {label}")
            except Exception as exc:
                a = futures[f]
                failures.append((a, exc))
                print(f"  FAILED job args={a}: {exc}")

    if failures:
        print(f"\n{len(failures)} job(s) failed.")
        return 1

    if not args.skip_combine:
        print("\nCombining channels...")
        from production.combine_channels import combine_for_point
        n_combined = 0
        for flavor in args.flavor:
            for m in masses:
                if combine_for_point(flavor, m) > 0:
                    n_combined += 1
        print(f"  combined files written for {n_combined} (flavor, mass) points")

    print(f"\nTotal wall time: {time.time()-t0:.1f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
