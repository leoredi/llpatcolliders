#!/usr/bin/env python3
"""
hnl/run_all.py

Parallel driver for the complete HNL production pipeline.

Per-flavor channels processed:

  meson sources    bottom, charm, bc       (FONLL + HNLCalc)
  light-meson      kaon                    (parametric + HNLCalc)
  induced tau      Ds, B+ -> tau -> N+X    (FONLL + HNLCalc, output: induced_tau/)
  prompt tau       W -> tau nu, Z -> tau tau -> tau -> N+X
                                           (MG5 SM + HNLCalc, output: tau/)
  electroweak      W/Z -> l N              (MG5 + HeavyN UFO, output: WZ/)

Each (flavor, channel) is one independent process. Prompt-tau is structurally
different: it runs *Stage 1* (one shared MG5 job that produces
vendored/tau_pool.csv) before the parallel pool, then fans out a *Stage 2*
job per flavor that decays the cached pool. The W/Z driver is default-on
since the December refactor; opt out with --no-wz.

Usage:
    python -m hnl.run_all                       # full grid, all cores, all channels
    python -m hnl.run_all --workers 4           # cap concurrency
    python -m hnl.run_all --n-pool 50000        # smaller meson/tau pools
    python -m hnl.run_all --masses 0.5 1.0 2.0  # subset of masses
    python -m hnl.run_all --no-wz --no-prompt-tau   # skip the MG5-heavy paths

Run from the hnl/ directory:
    cd hnl && python run_all.py
or as a module from the repo root with PYTHONPATH=hnl.
"""

from __future__ import annotations

import argparse
import os
import random
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

HNL_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(HNL_ROOT))
sys.path.insert(0, str(HNL_ROOT / "vendored" / "HNLCalc"))

from config_mass_grid import MASS_GRID


# Module-level worker functions (must be picklable for ProcessPoolExecutor).

def _worker_meson(flavor: str, channel: str, masses: list, n_pool: int, seed: int) -> str:
    """Process one (flavor, meson_channel) for all masses. channel ∈ {bottom, charm, bc}."""
    from production.decay_engine.generate_meson_csvs import (
        generate_pool, generate_bc_pool, process_channel,
    )
    rng = np.random.default_rng(seed)
    random.seed(seed)  # HNLCalc's 3-body BR integrator uses stdlib random, not numpy
    if channel == "bc":
        pool = generate_bc_pool(n_pool, rng)
        sigma = 0.0
    else:
        pool, sigma = generate_pool(channel, n_pool, rng)
    process_channel(flavor, channel, pool, sigma, masses, rng)
    return f"meson {flavor:>5s} {channel}"


def _worker_induced_tau(flavor: str, masses: list, n_pool: int, seed: int) -> str:
    """Induced-tau channel (Ds, B+ -> tau -> N+X). Output: induced_tau/."""
    from production.decay_engine.generate_induced_tau import build_tau_pool, process_flavor
    rng = np.random.default_rng(seed)
    random.seed(seed)
    tau_4v, tau_w = build_tau_pool(n_pool, rng)
    # Pass the full mass list, not a pre-filtered subset. process_flavor writes
    # an empty sentinel CSV for m_N >= m_tau so combine_channels sees a fresh
    # "channel closed" marker every run instead of inheriting stale rows from
    # a previous run on a wider grid.
    process_flavor(flavor, tau_4v, tau_w, masses, rng)
    return f"itau  {flavor:>5s}"


def _worker_kaon(flavor: str, masses: list, n_pool: int, seed: int) -> str:
    """Kaon channel (K+ -> l N) for one flavor across all masses."""
    from production.decay_engine.generate_kaon_csvs import sample_kaon_4vectors, process_flavor
    rng = np.random.default_rng(seed)
    random.seed(seed)
    pool = sample_kaon_4vectors(n_pool, rng)
    process_flavor(flavor, pool, masses, rng)
    return f"kaon  {flavor:>5s}"


def _worker_prompt_tau_stage2(flavor: str, masses: list, seed: int) -> str:
    """Stage 2 of prompt tau: decay vendored/tau_pool.csv → output/.../tau/.

    Expects vendored/tau_pool.csv to already exist (Stage 1 has run, or
    the cached pool was vendored). Output: tau/.
    """
    from production.madgraph.run_tau_production import _load_pool, process_flavor
    rng = np.random.default_rng(seed)
    random.seed(seed)
    pool_w, pool_E, pool_px, pool_py, pool_pz, pool_origin = _load_pool()
    # Pass the full mass list; process_flavor writes the empty sentinel for
    # m_N >= m_tau (see _worker_induced_tau for the rationale).
    process_flavor(flavor, pool_w, pool_E, pool_px, pool_py, pool_pz, pool_origin,
                   masses, rng)
    return f"ptau  {flavor:>5s}"


def _worker_wz(flavor: str, masses: list, n_events: int, nb_core: int) -> str:
    """W/Z → l N MadGraph pipeline for one flavor across all masses.

    Default-on since the December refactor. Each mass spawns a MadGraph
    generate_events run (process dir cached per flavor inside the worker).
    We iterate through every mass so the failure report lists *all* bad
    points, then raise once at the end. Raising surfaces the failure to the
    parallel driver's ``failures`` list (the previous "count successes
    silently" path let combine_channels happily skip absent WZ CSVs and the
    driver still exited 0).
    """
    from production.madgraph.run_wz_production import run_single_point
    failed = []
    for m in masses:
        if not run_single_point(flavor, m, n_events, nb_core=nb_core):
            failed.append(m)
    if failed:
        raise RuntimeError(f"W/Z failed for {flavor} masses: {failed}")
    return f"wz    {flavor:>5s} ({len(masses)}/{len(masses)} points)"


def main():
    ap = argparse.ArgumentParser(description="Parallel HNL production driver")
    ap.add_argument("--flavor", choices=["Ue", "Umu", "Utau"], nargs="+",
                    default=["Ue", "Umu", "Utau"])
    ap.add_argument("--n-pool", type=int, default=100_000,
                    help="Sampled mesons/taus per (flavor, channel) job (default: 100k)")
    ap.add_argument("--workers", type=int, default=None,
                    help="Max concurrent processes (default: os.cpu_count())")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--masses", type=float, nargs="+", default=None,
                    help="Custom mass list (default: full grid)")
    ap.add_argument("--skip-combine", action="store_true",
                    help="Skip the final combine_channels step")
    ap.add_argument("--no-wz", action="store_true",
                    help="Skip the W/Z -> l N MadGraph channel (default on)")
    ap.add_argument("--wz-nevents", type=int, default=None,
                    help="Events per (flavor, mass) MadGraph point "
                         "(default: config N_EVENTS_DEFAULT)")
    ap.add_argument("--wz-nb-core", type=int, default=1,
                    help="CPU cores per MadGraph generate_events job")
    ap.add_argument("--no-prompt-tau", action="store_true",
                    help="Skip the prompt-tau channel (Stage 1 MG5 + Stage 2 decay)")
    ap.add_argument("--prompt-tau-nevents", type=int, default=None,
                    help="Stage 1 tau-pool size (default: config N_EVENTS_DEFAULT). "
                         "Ignored if vendored/tau_pool.csv already exists.")
    ap.add_argument("--prompt-tau-nb-core", type=int, default=1,
                    help="CPU cores for Stage 1 MG5 generate_events")
    args = ap.parse_args()

    workers = args.workers or os.cpu_count() or 1
    masses = args.masses if args.masses else MASS_GRID
    meson_channels = ["bottom", "charm", "bc"]

    # Run prompt-tau Stage 1 up-front (single shared MG5 job, gated on cache).
    # Stage 2 is a per-flavor numpy job that fans out alongside the meson workers.
    if not args.no_prompt_tau:
        from production.madgraph.run_tau_production import POOL_CSV, generate_tau_pool
        from config_mass_grid import N_EVENTS_DEFAULT
        pt_nevents = args.prompt_tau_nevents or N_EVENTS_DEFAULT

        regen_reason = None
        if not POOL_CSV.exists() or POOL_CSV.stat().st_size == 0:
            regen_reason = "no cached pool"
        else:
            # Guard against silently driving production with a smoke-sized
            # pool (the LHE typically yields ~1.5–2 taus per MG5 event after
            # the W → τν + DY → ττ split, but we conservatively require
            # row_count >= 0.5 × nevents).
            with open(POOL_CSV) as f:
                n_rows = sum(1 for _ in f)
            threshold = max(int(0.5 * pt_nevents), 1)
            if n_rows < threshold:
                regen_reason = (f"cached pool has {n_rows} rows < threshold "
                                f"{threshold} for --prompt-tau-nevents={pt_nevents}")
            else:
                print(f"Prompt-tau Stage 1: reusing cached pool at {POOL_CSV} "
                      f"({n_rows} rows)")

        if regen_reason is not None:
            print(f"Prompt-tau Stage 1: regenerating ({regen_reason}); "
                  f"building pool ({pt_nevents} events) ...")
            ok = generate_tau_pool(pt_nevents, nb_core=args.prompt_tau_nb_core)
            if not ok:
                # Stage 1 is the gating step for the entire prompt-tau channel.
                # Continuing here would let combine_channels silently treat the
                # missing tau/ folder as zero and the driver exit 0; abort
                # instead so the failure is visible.
                print("Prompt-tau Stage 1 failed; aborting.")
                return 1

    jobs = []   # list of (fn, args_tuple) used by submit loop
    seed = args.seed
    for flavor in args.flavor:
        for ch in meson_channels:
            jobs.append((_worker_meson, (flavor, ch, masses, args.n_pool, seed)))
            seed += 1
        jobs.append((_worker_induced_tau, (flavor, masses, args.n_pool, seed)))
        seed += 1
        jobs.append((_worker_kaon, (flavor, masses, args.n_pool, seed)))
        seed += 1
        if not args.no_prompt_tau:
            jobs.append((_worker_prompt_tau_stage2, (flavor, masses, seed)))
            seed += 1

    if not args.no_wz:
        from config_mass_grid import N_EVENTS_DEFAULT
        wz_nevents = args.wz_nevents or N_EVENTS_DEFAULT
        for flavor in args.flavor:
            jobs.append((_worker_wz, (flavor, masses, wz_nevents, args.wz_nb_core)))

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
        # active_channels reflects what this invocation actually scheduled, so
        # combine's strict mode does not flag --no-wz / --no-prompt-tau output
        # trees as missing the WZ/tau directories.
        active_channels = ["Bmeson", "Dmeson", "Bc", "induced_tau", "Kmeson"]
        if not args.no_prompt_tau:
            active_channels.append("tau")
        if not args.no_wz:
            active_channels.append("WZ")
        n_combined = 0
        for flavor in args.flavor:
            for m in masses:
                if combine_for_point(flavor, m, channels=active_channels,
                                     strict=True) > 0:
                    n_combined += 1
        print(f"  combined files written for {n_combined} (flavor, mass) points")

    print(f"\nTotal wall time: {time.time()-t0:.1f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
