#!/usr/bin/env python3
"""BC4 driver: produce B -> K S four-vectors, scan the coupling, plot the island.

Coupling-controlled, model-complete (like the HNL, not the model-agnostic
(BR, c*tau) scan the higgs/ BC5 analysis uses): for each (m_S, sin^2 theta) the
*single* coupling sets the production yield, the lifetime c*tau and the visible
BR simultaneously, and we require N_signal >= 3 at 3000 fb^-1 background-free.
The result is a CLOSED ISLAND in (m_S, sin^2 theta): a lower edge (too little
production) and an upper edge (decays before reaching GRENDEL).

Usage:
    python -m scalar.run_sensitivity                 # full grid; produces only missing
                                                     # CSVs, reuses cached ray-casts
    python -m scalar.run_sensitivity --force-produce # regenerate the four-vector CSVs
    python -m scalar.run_sensitivity --plot-only     # re-plot from the island CSV
    python -m scalar.run_sensitivity --masses 0.5 1 2 --n-pool 100000
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_SCALAR_ROOT = Path(__file__).resolve().parent
_REPO_ROOT = _SCALAR_ROOT.parent
for _p in (str(_REPO_ROOT), str(_REPO_ROOT / "hnl"), str(_REPO_ROOT / "higgs")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from scalar import production as prod                # noqa: E402
from scalar import acceptance as acc                 # noqa: E402
from scalar.plot_exclusion import plot_island        # noqa: E402
from grendel_geometry import mesh_fiducial            # noqa: E402
from production.fonll.fonll_parser import get_sigma_total       # noqa: E402
from production.fonll.meson_sampler import sample_meson_4vectors  # noqa: E402

# sin^2 theta scan range (log10).  BC4 reach lives ~1e-11 - 1e-6.
LOG_S2T_MIN, LOG_S2T_MAX, N_S2T = -12.0, -2.0, 200

OUT_DIR = _SCALAR_ROOT / "tmp"
VEC_DIR = OUT_DIR / "llp_4vectors"
ISLAND_CSV = OUT_DIR / "bc4_island.csv"


def run(masses, n_pool, seed, n_samples, force_produce=False):
    """Produce-if-missing + cached ray-casts (the BC10/HNL pattern): the
    four-vector CSVs are only (re)generated for masses that have none, or for
    all masses with ``force_produce``; untouched CSVs keep their mtime so the
    geometry cache in ``acceptance._geometry`` stays valid and a re-scan only
    redoes the decay MC + reconstruction + coupling scan."""
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    s2t_grid = np.logspace(LOG_S2T_MIN, LOG_S2T_MAX, N_S2T)
    rng = np.random.default_rng(seed)

    if force_produce:
        to_produce = list(masses)
    else:
        to_produce = [m for m in masses
                      if not (VEC_DIR / f"mS_{prod._mass_label(m)}.csv").exists()]
    if to_produce:
        sigma_bottom = get_sigma_total("bottom")
        pool = sample_meson_4vectors(n_pool, "bottom", rng=rng)
        print(f"sigma_FONLL(bottom) = {sigma_bottom:.3e} pb; "
              f"producing {len(to_produce)}/{len(masses)} masses")
        for m_S in to_produce:
            prod.write_scalar_csv(m_S, VEC_DIR, n_pool, rng,
                                  sigma_bottom=sigma_bottom, pool=pool)
    else:
        print(f"reusing existing four-vector CSVs for all {len(masses)} masses "
              "(--force-produce to regenerate)")

    rows = []
    for i, m_S in enumerate(masses):
        csv = VEC_DIR / f"mS_{prod._mass_label(m_S)}.csv"
        out = acc.process_mass_point(
            m_S, mesh_fiducial, csv, s2t_grid, n_samples=n_samples,
            rng=np.random.default_rng(1000 + i))
        if out is None:
            continue
        band = out[0] if isinstance(out, tuple) else out
        rows.append(band)
        tag = "SENS" if band.get("has_sensitivity") else "----"
        print(f"  [{i+1}/{len(masses)}] m_S={m_S:.3f} {tag} "
              f"peak_N={band['peak_N']:.1f} "
              f"island=[{band['u2_min']:.2e}, {band['u2_max']:.2e}]", flush=True)

    df = pd.DataFrame(rows).sort_values("mass_GeV")
    df.to_csv(ISLAND_CSV, index=False)
    print(f"\nIsland CSV: {ISLAND_CSV}")
    _summarize(df)
    plot_island(ISLAND_CSV, OUT_DIR)
    return df


def _summarize(df):
    sens = df[df["has_sensitivity"].fillna(False)]
    if sens.empty:
        print("No sensitive mass points.")
        return
    print(f"Sensitive masses: {len(sens)}  "
          f"m_S in [{sens['mass_GeV'].min():.3f}, {sens['mass_GeV'].max():.3f}] GeV")
    print(f"Deepest lower edge (sin^2 theta): {np.nanmin(sens['u2_min']):.2e}")
    print(f"Highest upper edge (sin^2 theta): {np.nanmax(sens['u2_max']):.2e}")


def main(argv=None):
    p = argparse.ArgumentParser(description="BC4 dark-scalar GRENDEL sensitivity")
    p.add_argument("--masses", type=float, nargs="+", default=None)
    p.add_argument("--n-pool", type=int, default=prod.N_POOL_DEFAULT)
    p.add_argument("--n-samples", type=int, default=100)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--plot-only", action="store_true")
    p.add_argument("--force-produce", action="store_true",
                   help="regenerate the four-vector CSVs even where they exist "
                        "(default: produce only missing masses and reuse the "
                        "geometry cache)")
    args = p.parse_args(argv)

    if args.plot_only:
        plot_island(ISLAND_CSV, OUT_DIR)
        return 0

    masses = args.masses if args.masses else prod.MASS_GRID
    run(masses, args.n_pool, args.seed, args.n_samples,
        force_produce=args.force_produce)
    return 0


if __name__ == "__main__":
    sys.exit(main())
