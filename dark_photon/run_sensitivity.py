"""BC1 dark-photon GRENDEL sensitivity driver.

Produces A' four-vector CSVs, runs the shared acceptance + eps^2 scan per
mass, writes the exclusion island CSV, and plots the (m_A', eps^2) contour.
The canonical grid is a 1-MeV scan over 0.02--0.20 GeV; custom high-mass
Drell--Yan controls remain available through ``--masses``.

    python -m dark_photon.run_sensitivity                 # full grid
    python -m dark_photon.run_sensitivity --masses 0.3    # one mass
    python -m dark_photon.run_sensitivity --plot-only

Nothing here regenerates the committed Pythia spectra; four-vector CSVs are
resampled from them (fast, pure Python) only for masses that have none, or for
all masses with ``--force-produce``.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_ROOT = Path(__file__).resolve().parent
_REPO_ROOT = _ROOT.parent
# Repo root (so ``dark_photon`` is an importable package) + shared hnl/higgs;
# never dark_photon/ itself (its production.py would shadow hnl's production/).
for _p in (str(_REPO_ROOT), str(_REPO_ROOT / "hnl"), str(_REPO_ROOT / "higgs")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from dark_photon import model                            # noqa: E402,F401
from dark_photon import production as prod               # noqa: E402
from dark_photon import acceptance as acc                # noqa: E402
from dark_photon.mass_grid import MASS_GRID              # noqa: E402
from dark_photon.plot_exclusion import plot_island       # noqa: E402
from grendel_geometry import mesh_fiducial               # noqa: E402

# eps^2 scan range (log10). Dark-photon reach lives ~1e-17 -- 1e-10.
LOG_EPS2_MIN, LOG_EPS2_MAX, N_EPS2 = -18.0, -8.0, 220

OUT_DIR = _ROOT / "tmp"
VEC_DIR = OUT_DIR / "llp_4vectors"
ISLAND_CSV = OUT_DIR / "bc1_island.csv"
N_PER_CHANNEL_DEFAULT = 200_000


def _reconstruction_seed(m_a):
    """Mass-derived seed, invariant under grid refinement and subset scans."""
    return 100_000 + int(round(float(m_a) * 1000))


def _production_seed(m_a):
    return 500_000 + int(round(float(m_a) * 1000))


def _write_checkpoint(rows, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame(rows).sort_values("mass_GeV")
    tmp = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(tmp, index=False)
    tmp.replace(path)


def _production_is_stale(m_a, csv):
    """Whether a cached vector CSV is absent, empty, or older than its input.

    This matters when a spectrum grid is refined or regenerated: existence
    alone is not a valid cache key.  In particular, a zero-byte CSV left by a
    previously unsupported mass must not suppress production forever.
    """
    csv = Path(csv)
    if not csv.exists() or csv.stat().st_size == 0:
        return True
    inputs = []
    if m_a < model.M_OMEGA - model.M_PI0:
        inputs.append(prod.MESON_SPECTRUM)
    dy = prod.DY_SPECTRUM_DIR / f"dy_{prod._mass_label(m_a)}.npz"
    if dy.exists():
        inputs.append(dy)
    return any(path.exists() and path.stat().st_mtime > csv.stat().st_mtime
               for path in inputs)


def run(masses, n_per_channel, n_samples=100, force_produce=False,
        vector_dir=VEC_DIR, output=ISLAND_CSV):
    vector_dir = Path(vector_dir)
    eps2_grid = np.logspace(LOG_EPS2_MIN, LOG_EPS2_MAX, N_EPS2)

    for m_a in masses:
        csv = vector_dir / f"mA_{prod._mass_label(m_a)}.csv"
        if force_produce or _production_is_stale(m_a, csv):
            rng_p = np.random.default_rng(_production_seed(m_a))
            path, n, sig = prod.write_csv(m_a, vector_dir, n_per_channel, rng_p)
            if n == 0:
                print(f"  m_A={m_a:.3f}: no open production channel, skip")

    rows = []
    for i, m_a in enumerate(masses):
        csv = vector_dir / f"mA_{prod._mass_label(m_a)}.csv"
        out = acc.process_mass_point(
            m_a, mesh_fiducial, csv, eps2_grid, n_samples=n_samples,
            rng=np.random.default_rng(_reconstruction_seed(m_a)))
        if out is None:
            continue
        rows.append(out)
        _write_checkpoint(rows, output)
        tag = "SENS" if out.get("has_sensitivity") else "----"
        emin, emax = out.get("eps2_min", np.nan), out.get("eps2_max", np.nan)
        print(f"  [{i+1}/{len(masses)}] m_A={m_a:.3f} {tag} "
              f"peak_N={out['peak_N']:.1f} island=[{emin:.2e}, {emax:.2e}]",
              flush=True)

    if not rows:
        print("No results.")
        return None
    df = pd.DataFrame(rows).sort_values("mass_GeV")
    _write_checkpoint(rows, output)
    print(f"\nIsland CSV: {output}")
    _summarize(df)
    plot_island(output, output.parent)
    return df


def _summarize(df):
    sens = df[df["has_sensitivity"].fillna(False)]
    if sens.empty:
        print("No sensitive mass points.")
        return
    print(f"Sensitive masses: {len(sens)}  "
          f"m_A in [{sens['mass_GeV'].min():.3f}, {sens['mass_GeV'].max():.3f}] GeV")
    print(f"Deepest lower edge (eps^2): {np.nanmin(sens['eps2_min']):.2e}")
    print(f"Highest upper edge (eps^2): {np.nanmax(sens['eps2_max']):.2e}")


def main(argv=None):
    p = argparse.ArgumentParser(description="BC1 dark-photon GRENDEL sensitivity")
    p.add_argument("--masses", type=float, nargs="+", default=None)
    p.add_argument("--n-per-channel", type=int, default=N_PER_CHANNEL_DEFAULT)
    p.add_argument("--n-samples", type=int, default=100)
    p.add_argument("--vector-dir", type=Path, default=VEC_DIR)
    p.add_argument("--output", type=Path, default=ISLAND_CSV)
    p.add_argument("--force-produce", action="store_true")
    p.add_argument("--plot-only", action="store_true")
    args = p.parse_args(argv)

    if args.plot_only:
        plot_island(args.output, args.output.parent)
        return 0

    masses = args.masses if args.masses else MASS_GRID
    run(masses, args.n_per_channel, n_samples=args.n_samples,
        force_produce=args.force_produce, vector_dir=args.vector_dir,
        output=args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
