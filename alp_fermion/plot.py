"""BC10 exclusion-island plot in the PBC (m_a, 1/f) plane.

Draws the GRENDEL closed island (lower edge = too little production, upper edge
= decays before PX56) and overlays competitor / existing-bound curves when a
digitized reference table is available.  The unified-calc competitor
(arXiv:2311.00507) and the existing-bound contours are external data products;
point ``--overlay-csv`` at a long-format CSV with columns

    curve, m_a_GeV, invf_GeV_inv

(one row per point; ``curve`` names the experiment/bound) to overlay them.  See
EXTERNAL_INPUTS_NEEDED.md.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PLOT_INVF_MIN = 1e-9      # GeV^-1
PLOT_INVF_MAX = 1e-2
PLOT_MA_MIN = 0.2
PLOT_MA_MAX = 5.0


def _segments(df):
    """Contiguous runs of sensitive masses (no bridging insensitive points)."""
    s = df["has_sensitivity"].fillna(False).astype(bool)
    grp = s.ne(s.shift(fill_value=False)).cumsum()
    for _, seg in df[s].groupby(grp[s]):
        yield seg.sort_values("mass_GeV")


def plot_island(results_csv, output_dir=None, basename="bc10_island",
                overlay_csv=None):
    results_csv = Path(results_csv)
    df = pd.read_csv(results_csv).sort_values("mass_GeV")
    output_dir = Path(output_dir) if output_dir else results_csv.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(7.5, 6))
    plotted = False
    for seg in _segments(df):
        m = seg["mass_GeV"].to_numpy(float)
        lo = seg["invf_min"].to_numpy(float)
        hi = seg["invf_max"].to_numpy(float)
        lo_open = seg["invf_min_open"].to_numpy(bool)
        hi_open = seg["invf_max_open"].to_numpy(bool)
        lo_fill = np.where(lo_open | ~np.isfinite(lo), PLOT_INVF_MIN, lo)
        hi_fill = np.where(hi_open | ~np.isfinite(hi), PLOT_INVF_MAX, hi)
        ax.fill_between(m, lo_fill, hi_fill, color="red", alpha=0.25,
                        zorder=5, label="GRENDEL (PX56), 3000 fb$^{-1}$"
                        if not plotted else None)
        ax.plot(m, np.where(lo_open, np.nan, lo), "r-", lw=1.8, zorder=6)
        ax.plot(m, np.where(hi_open, np.nan, hi), "r-", lw=1.8, zorder=6)
        for end in (0, -1):
            if not lo_open[end] and not hi_open[end] and np.isfinite(lo[end]) and np.isfinite(hi[end]):
                ax.plot([m[end], m[end]], [lo[end], hi[end]], "r-", lw=1.8, zorder=6)
        plotted = True

    if overlay_csv and Path(overlay_csv).exists():
        ov = pd.read_csv(overlay_csv)
        for name, g in ov.groupby("curve"):
            g = g.sort_values("m_a_GeV")
            ax.plot(g["m_a_GeV"], g["invf_GeV_inv"], lw=1.5, ls="--",
                    label=name, zorder=4)

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(PLOT_MA_MIN, PLOT_MA_MAX)
    ax.set_ylim(PLOT_INVF_MIN, PLOT_INVF_MAX)
    ax.set_xlabel(r"$m_a$ [GeV]", fontsize=13)
    ax.set_ylabel(r"$1/f$ [GeV$^{-1}$]  ($c_f=1$, fermiophilic ALP)", fontsize=13)
    ax.set_title("BC10: fermiophilic ALP, $B\\to K\\,a$, $a\\to$ charged tracks",
                 fontsize=12)
    ax.grid(True, which="both", alpha=0.2, lw=0.5)
    if plotted or overlay_csv:
        ax.legend(fontsize=10, loc="lower right")
    else:
        ax.text(0.5, 0.5, "no closed island (rate-starved)",
                ha="center", va="center", transform=ax.transAxes, fontsize=12)

    out_png = output_dir / f"{basename}.png"
    out_pdf = output_dir / f"{basename}.pdf"
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out_png}")
    print(f"Saved {out_pdf}")
    return out_png


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("results_csv")
    ap.add_argument("--overlay-csv", default=None)
    ap.add_argument("--out-dir", default=None)
    a = ap.parse_args()
    plot_island(a.results_csv, a.out_dir, overlay_csv=a.overlay_csv)
