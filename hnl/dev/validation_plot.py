#!/usr/bin/env python3
"""
hnl/dev/validation_plot.py

Overlay our combined HNL production yield (sigma at U^2 = 1, in pb) versus the
MATHUSLA RHN reference 4-vector files, per flavor. This is the validation plot
called for in the production-review next-steps: a factor-of-2 disagreement on
the common channels would flag a quark+antiquark double-count. Dev-only — see
hnl/dev/README.md.

Our data:   output/llp_4vectors/{flavor}/{channel}/mN_{label}.csv  (weight col 0)
Reference:  llpatcolliders_FONLL/vendored/MATHUSLA_LLPfiles_RHN_U{e,mu,tau}/
            All_RHN_*/RHN_*_LLPweight4vector{Bmeson,Dmeson,Tau,WZ}list_mN_*.csv

Only the channels common to both pipelines (Bmeson, Dmeson, induced_tau —
matched against the MATHUSLA "Tau" reference, which is induced-only) are
summed for the apples-to-apples overlay; our full total (incl. Bc, Kmeson,
prompt-tau, WZ) is drawn too.

Usage:
    python dev/validation_plot.py
    python dev/validation_plot.py --ref /path/to/MATHUSLA_root_dir
"""

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HNL_ROOT = Path(__file__).resolve().parent.parent  # hnl/dev/ -> hnl/
sys.path.insert(0, str(HNL_ROOT))

from config_mass_grid import MASS_GRID, format_mass_for_filename

OUTPUT_BASE = HNL_ROOT / "output" / "llp_4vectors"
PLOT_DIR = HNL_ROOT / "output" / "debug_plots"

FLAVORS = ["Ue", "Umu", "Utau"]
OUR_CHANNELS = ["Bmeson", "Dmeson", "Bc", "tau", "induced_tau", "Kmeson", "WZ"]
# COMMON_CHANNELS are those present in the MATHUSLA reference. The MATHUSLA
# 'Tau' files are induced-only (no prompt-tau), so we match against
# induced_tau on our side; comparing against our 'tau' would mix in
# prompt-tau contributions the reference doesn't have.
COMMON_CHANNELS = ["Bmeson", "Dmeson", "induced_tau"]

# Map our channel folder -> reference channel token in the filename.
REF_TOKEN = {
    "Bmeson": "Bmeson", "Dmeson": "Dmeson",
    "induced_tau": "Tau", "WZ": "WZ",
}
FLAVOR_SUFFIX = {"Ue": "Ue", "Umu": "Umu", "Utau": "Utau"}


def _sum_weights(csv_path):
    if not csv_path.exists() or csv_path.stat().st_size == 0:
        return 0.0
    try:
        data = np.loadtxt(csv_path, delimiter=",")
    except ValueError:
        return 0.0
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.size == 0:
        return 0.0
    return float(data[:, 0].sum())


def our_yield(flavor, channels):
    """Sum-of-weights vs mass for the given channels of our output."""
    masses, ys = [], []
    for m in MASS_GRID:
        label = format_mass_for_filename(m)
        s = sum(_sum_weights(OUTPUT_BASE / flavor / ch / f"mN_{label}.csv")
                for ch in channels)
        if s > 0:
            masses.append(m)
            ys.append(s)
    return np.array(masses), np.array(ys)


def _resolve_ref_base(ref_arg):
    if ref_arg:
        return Path(ref_arg)
    # Sibling llpatcolliders_FONLL checkout (one or two levels above hnl/).
    rel = Path("llpatcolliders_FONLL") / "vendored"
    for base in (HNL_ROOT.parent, HNL_ROOT.parent.parent):
        if (base / rel).is_dir():
            return base / rel
    return HNL_ROOT.parent.parent / rel


def reference_yield(ref_base, flavor, channels):
    """Sum-of-weights vs mass for the MATHUSLA reference, common channels."""
    suffix = FLAVOR_SUFFIX[flavor]
    flavor_dir = ref_base / f"MATHUSLA_LLPfiles_RHN_{suffix}"
    if not flavor_dir.is_dir():
        return np.array([]), np.array([])
    tokens = {REF_TOKEN[c] for c in channels if c in REF_TOKEN}
    pat = re.compile(rf"RHN_{suffix}_LLPweight4vector(\w+?)list_mN_([0-9.]+)\.csv$")

    per_mass = {}
    for csv in flavor_dir.rglob("RHN_*_LLPweight4vector*list_mN_*.csv"):
        mobj = pat.search(csv.name)
        if not mobj:
            continue
        ch_token, mass_str = mobj.group(1), mobj.group(2)
        if ch_token not in tokens:
            continue
        try:
            mass = float(mass_str)
        except ValueError:
            continue
        per_mass[mass] = per_mass.get(mass, 0.0) + _sum_weights(csv)

    if not per_mass:
        return np.array([]), np.array([])
    masses = np.array(sorted(per_mass))
    ys = np.array([per_mass[m] for m in masses])
    keep = ys > 0
    return masses[keep], ys[keep]


def main():
    ap = argparse.ArgumentParser(description="Validation overlay vs MATHUSLA RHN reference")
    ap.add_argument("--ref", default=None, help="MATHUSLA reference root (the 'vendored' dir)")
    ap.add_argument("--out", default=str(PLOT_DIR / "yield_vs_mass_validation.png"))
    args = ap.parse_args()

    ref_base = _resolve_ref_base(args.ref)
    print(f"Reference base: {ref_base} (exists: {ref_base.is_dir()})")

    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 3, figsize=(16, 9), sharex="col")

    for col, flavor in enumerate(FLAVORS):
        ax, axr = axes[0, col], axes[1, col]
        m_tot, y_tot = our_yield(flavor, OUR_CHANNELS)
        m_com, y_com = our_yield(flavor, COMMON_CHANNELS)
        m_ref, y_ref = reference_yield(ref_base, flavor, COMMON_CHANNELS)

        # Single global rescale: bring the reference onto our B+D+tau curve via
        # the median per-mass ratio (interpolated in log-log over the overlap).
        scale, ratio_spread = np.nan, np.nan
        if len(m_com) and len(m_ref):
            lo, hi = max(m_com.min(), m_ref.min()), min(m_com.max(), m_ref.max())
            mask = (m_ref >= lo) & (m_ref <= hi)
            if mask.sum() >= 3:
                our_interp = np.exp(np.interp(np.log(m_ref[mask]),
                                              np.log(m_com), np.log(y_com)))
                r = y_ref[mask] / our_interp
                scale = float(np.median(r))
                # Spread of ratio/median: ~1 everywhere => pure constant offset.
                ratio_spread = float(np.std(np.log10(r / scale)))

        if len(m_tot):
            ax.plot(m_tot, y_tot, "-", color="C0", lw=1.8, label="ours: all channels")
        if len(m_com):
            ax.plot(m_com, y_com, "--", color="C1", lw=1.8, label="ours: B+D+tau")
        if len(m_ref):
            ax.plot(m_ref, y_ref, "o", color="k", ms=4, mfc="none",
                    label="MATHUSLA ref: B+D+Tau")
        ax.set_yscale("log")
        ax.set_title(flavor)
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize=8)

        # Bottom row: reference rescaled by 1/scale overlaid on our B+D+tau.
        if len(m_com):
            axr.plot(m_com, y_com, "--", color="C1", lw=1.8, label="ours: B+D+tau")
        if len(m_ref) and np.isfinite(scale):
            axr.plot(m_ref, y_ref / scale, "o", color="k", ms=4, mfc="none",
                     label=f"ref / {scale:.2e}")
            axr.text(0.04, 0.08,
                     f"shape scatter (dex): {ratio_spread:.2f}",
                     transform=axr.transAxes, fontsize=9,
                     bbox=dict(boxstyle="round", fc="w", alpha=0.7))
        axr.set_yscale("log")
        axr.set_xlabel(r"$m_N$ [GeV]")
        axr.grid(True, which="both", alpha=0.3)
        axr.legend(fontsize=8)

        print(f"{flavor}: global ref/ours scale = {scale:.3e}, "
              f"shape scatter = {ratio_spread:.3f} dex")

    axes[0, 0].set_ylabel(r"$\sum_i w_i$  [pb at $|U|^2=1$]")
    axes[1, 0].set_ylabel("rescaled (shape comparison)")
    fig.suptitle("HNL production yield vs MATHUSLA RHN reference "
                 "(top: absolute; bottom: reference rescaled by a single constant)")
    fig.tight_layout()
    fig.savefig(args.out, dpi=130)
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
