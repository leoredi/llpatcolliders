"""BC1 dark-photon exclusion plot with the established SensCalc envelope."""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PLOT_EPS2_MIN, PLOT_EPS2_MAX = 1e-16, 1e-8
PLOT_M_MIN, PLOT_M_MAX = 0.02, 0.20
SIGNAL_THRESHOLD = 3.0
EXISTING_CONSTRAINTS = (
    Path(__file__).resolve().parent / "data" / "constraints"
    / "existing_exclusion_senscalc.npz"
)


def _segments(df):
    """Yield resolved sensitive intervals, not zero-width singleton bins."""
    sens = df["has_sensitivity"].fillna(False).astype(bool)
    grp = sens.ne(sens.shift(fill_value=False)).cumsum()
    for _, seg in df[sens].groupby(grp[sens]):
        if len(seg) >= 2:
            yield seg.sort_values("mass_GeV")


def _plot_existing_exclusions(ax, path=EXISTING_CONSTRAINTS):
    """Draw the SensCalc BC1 baseline envelope, converting eps to eps^2."""
    path = Path(path)
    if not path.exists():
        return False
    data = np.load(path, allow_pickle=False)
    n_polygons = int(data["n_polygons"])
    labelled = False
    for index in range(n_polygons):
        points = np.asarray(data[f"polygon_{index}"], float)
        mass = points[:, 0]
        eps2 = points[:, 1] ** 2
        label = "Existing exclusions (SensCalc)" if not labelled else None
        if bool(data[f"closed_{index}"]):
            ax.fill(mass, eps2, color="0.72", edgecolor="0.48", lw=0.8,
                    alpha=0.75, label=label, zorder=1)
        else:
            ax.fill_between(mass, eps2, PLOT_EPS2_MAX, color="0.72",
                            edgecolor="0.48", lw=0.8, alpha=0.75,
                            label=label, zorder=1)
        labelled = True
    return labelled


def plot_island(island_csv, output_dir=None, basename="bc1_exclusion"):
    island_csv = Path(island_csv)
    df = pd.read_csv(island_csv).sort_values("mass_GeV")
    output_dir = Path(output_dir) if output_dir else island_csv.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 6))
    existing_plotted = _plot_existing_exclusions(ax)
    plotted = False
    for seg in _segments(df):
        mass = seg["mass_GeV"].to_numpy(float)
        e_min = seg["eps2_min"].to_numpy(float)
        e_max = seg["eps2_max"].to_numpy(float)
        lo_open = seg.get("eps2_min_open", pd.Series(False, index=seg.index)).fillna(False).to_numpy(bool)
        hi_open = seg.get("eps2_max_open", pd.Series(False, index=seg.index)).fillna(False).to_numpy(bool)
        lower = np.where(lo_open, PLOT_EPS2_MIN, e_min)
        upper = np.where(hi_open, PLOT_EPS2_MAX, e_max)
        ax.fill_between(mass, lower, upper, color="red", alpha=0.30,
                        label="GRENDEL (3000 fb$^{-1}$)" if not plotted else None,
                        zorder=5)
        ax.plot(mass, np.where(lo_open, np.nan, e_min), "r-", lw=1.8, zorder=6)
        ax.plot(mass, np.where(hi_open, np.nan, e_max), "r-", lw=1.8, zorder=6)
        for end in (0, -1):
            if not lo_open[end] and not hi_open[end]:
                ax.plot([mass[end]] * 2, [e_min[end], e_max[end]], "r-", lw=1.8, zorder=6)
        plotted = True

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(PLOT_M_MIN, PLOT_M_MAX)
    ax.set_ylim(PLOT_EPS2_MIN, PLOT_EPS2_MAX)
    ax.set_xlabel(r"$m_{A'}$ [GeV]", fontsize=14)
    ax.set_ylabel(r"$\varepsilon^2$", fontsize=14)
    ax.set_title("BC1 dark photon — GRENDEL (PX56, CMS)", fontsize=13)
    ax.grid(True, which="both", alpha=0.2, lw=0.5)
    if plotted or existing_plotted:
        ax.legend(fontsize=9, loc="lower left")

    out_png = output_dir / f"{basename}.png"
    out_pdf = output_dir / f"{basename}.pdf"
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_png}\nSaved: {out_pdf}")
    return out_png


if __name__ == "__main__":
    import sys
    plot_island(sys.argv[1] if len(sys.argv) > 1 else "tmp/bc1_island.csv")
