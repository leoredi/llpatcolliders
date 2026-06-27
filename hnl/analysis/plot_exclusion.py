"""Publication-quality (m_N, U^2) sensitivity plot for GRENDEL."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path

PLOT_U2_MIN = 1e-12
PLOT_U2_MAX = 1e-1

_FLAVOR_LABEL = {
    "Ue":   r"$|U_e|^2$",
    "Umu":  r"$|U_\mu|^2$",
    "Utau": r"$|U_\tau|^2$",
}


def plot_exclusion(results_csv, output_dir=None, basename="hnl_exclusion",
                   band_csv=None):
    """Create the (m_N, U^2) exclusion figure from a sensitivity CSV.

    If ``band_csv`` (from ``analysis/combine_band.py``) is given, the FONLL
    theory-uncertainty ribbons are overlaid on each exclusion boundary while the
    central excluded-region fill is preserved.
    """
    results_csv = Path(results_csv)
    df = pd.read_csv(results_csv)
    if output_dir is None:
        output_dir = results_csv.parent
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    band_df = None
    if band_csv is not None and Path(band_csv).exists():
        band_df = pd.read_csv(band_csv)

    flavors_present = [f for f in ["Ue", "Umu", "Utau"]
                       if f in df["flavor"].unique()]
    if not flavors_present:
        print("No flavors in results CSV; nothing to plot.")
        return

    n = len(flavors_present)
    fig, axes = plt.subplots(1, n, figsize=(6 * n, 6), sharey=True,
                             squeeze=False)
    axes = axes[0]

    for idx, flavor in enumerate(flavors_present):
        _plot_single_panel(axes[idx], df, flavor, is_leftmost=(idx == 0),
                           band_df=band_df)

    fig.tight_layout(w_pad=2.0)
    out_pdf = output_dir / f"{basename}.pdf"
    out_png = output_dir / f"{basename}.png"
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_pdf}")
    print(f"Saved: {out_png}")


def _open_flags(valid, column, boundary_column, scan_limit):
    """Read open-edge flags, inferring them for pre-flag result CSVs."""
    if column in valid:
        return valid[column].fillna(False).astype(bool).to_numpy()
    boundary = valid[boundary_column].to_numpy(dtype=float)
    return np.isclose(boundary, scan_limit, rtol=1e-10, atol=0.0)


def _sensitive_segments(sel):
    """Yield contiguous sensitive rows without bridging insensitive masses."""
    is_sensitive = sel["has_sensitivity"].fillna(False).astype(bool)
    group = is_sensitive.ne(is_sensitive.shift(fill_value=False)).cumsum()
    for _, segment in sel[is_sensitive].groupby(group[is_sensitive]):
        yield segment


def _closure_vertex(mass, u2_min, u2_max):
    """High-mass pinch where the rising lower edge meets the falling upper edge.

    Linearly extrapolates log10(u2_min) and log10(u2_max) from the last two
    sensitive points to their crossing. Returns ``(m_star, u2_star)`` or
    ``None`` when the edges are not converging to a crossing just past the grid.
    """
    if len(mass) < 2:
        return None
    dm = mass[-1] - mass[-2]
    if dm <= 0:
        return None
    l1, l2 = np.log10(u2_min[-2]), np.log10(u2_min[-1])
    h1, h2 = np.log10(u2_max[-2]), np.log10(u2_max[-1])
    s_lo = (l2 - l1) / dm          # lower edge: rising  -> s_lo > 0
    s_hi = (h2 - h1) / dm          # upper edge: falling -> s_hi < 0
    if not (s_lo > 0 > s_hi):
        return None
    m_star = mass[-1] + (h2 - l2) / (s_lo - s_hi)
    if m_star <= mass[-1]:
        return None
    return m_star, 10.0 ** (l2 + s_lo * (m_star - mass[-1]))


def _closure_arc(m_last, u2min_last, u2max_last, m_star, u2_star, n=16):
    """Rounded closing nose instead of a sharp linear V.

    The two edges follow quarter-ellipses that meet at the pinch ``(m_star,
    u2_star)`` with a vertical tangent, so the island closes like the physical
    ~sqrt(m_close - m) narrowing rather than an unphysical arrowhead. Returns
    ``(mass, u2_min, u2_max)`` arrays for the arc (excluding the start point,
    which is already the last data point)."""
    yL, yU, yS = np.log10(u2min_last), np.log10(u2max_last), np.log10(u2_star)
    phi = np.linspace(0.0, np.pi / 2.0, n)[1:]
    x = m_last + (m_star - m_last) * np.sin(phi)
    lo = 10.0 ** (yS + (yL - yS) * np.cos(phi))
    hi = 10.0 ** (yS + (yU - yS) * np.cos(phi))
    return x, lo, hi


def _thin_marker_mask(open_mask, max_markers=12):
    """Sparse, evenly-spaced subset of an open-edge mask.

    When an exclusion boundary runs off the scan ceiling/floor across many dense
    grid points, drawing an open-arrow at every point piles them into an
    overlapping sawtooth. Keep at most ``max_markers`` evenly spaced ones so the
    row still reads as "boundary open beyond the axis" without the clutter.
    """
    idx = np.flatnonzero(open_mask)
    if idx.size <= max_markers:
        return open_mask
    keep = idx[np.linspace(0, idx.size - 1, max_markers).round().astype(int)]
    thinned = np.zeros_like(open_mask)
    thinned[keep] = True
    return thinned


def _band_ribbons(ax, flavor, mass, band_df, labelled):
    """Overlay theory-uncertainty ribbons on the two boundaries for one segment."""
    if band_df is None:
        return labelled
    bsel = band_df[band_df["flavor"] == flavor].set_index("mass_GeV")
    for col in ("u2_min", "u2_max"):
        lo = np.array([bsel[f"{col}_band_lo"].get(m, np.nan) for m in mass], dtype=float)
        hi = np.array([bsel[f"{col}_band_hi"].get(m, np.nan) for m in mass], dtype=float)
        good = np.isfinite(lo) & np.isfinite(hi)
        if not good.any():
            continue
        ax.fill_between(
            mass, np.where(good, lo, np.nan), np.where(good, hi, np.nan),
            alpha=0.35, color="orange", linewidth=0, zorder=4,
            label="FONLL theory band" if not labelled else None)
        labelled = True
    return labelled


def _plot_single_panel(ax, df, flavor, is_leftmost=True, band_df=None,
                       close_island=False):
    sel = df[df["flavor"] == flavor].sort_values("mass_GeV")
    plotted = False
    band_labelled = False

    for valid in _sensitive_segments(sel):
        mass = valid["mass_GeV"].to_numpy(dtype=float)
        u2_min = valid["u2_min"].to_numpy(dtype=float)
        u2_max = valid["u2_max"].to_numpy(dtype=float)
        min_open = _open_flags(
            valid, "u2_min_open", "u2_min", PLOT_U2_MIN)
        max_open = _open_flags(
            valid, "u2_max_open", "u2_max", PLOT_U2_MAX)

        # Close the high-mass dome to its interpolated pinch when the next grid
        # mass is insensitive (the island physically shut between the two), so
        # the contour ends in a point instead of a blunt residual-gap wall.
        if (close_island and len(mass) >= 2
                and not min_open[-1] and not max_open[-1]):
            later = sel[sel["mass_GeV"] > mass[-1]]
            if not later.empty and not bool(later.iloc[0]["has_sensitivity"]):
                cv = _closure_vertex(mass, u2_min, u2_max)
                if cv is not None and mass[-1] < cv[0] <= float(later.iloc[0]["mass_GeV"]):
                    am, alo, ahi = _closure_arc(mass[-1], u2_min[-1], u2_max[-1], cv[0], cv[1])
                    mass = np.append(mass, am)
                    u2_min = np.append(u2_min, alo)
                    u2_max = np.append(u2_max, ahi)
                    min_open = np.append(min_open, np.zeros(len(am), dtype=bool))
                    max_open = np.append(max_open, np.zeros(len(am), dtype=bool))

        # FONLL theory ribbons sit under the central fill so the excluded
        # region stays legible while the boundary uncertainty shows through.
        band_labelled = _band_ribbons(ax, flavor, mass, band_df, band_labelled)

        lower_fill = np.where(min_open, PLOT_U2_MIN, u2_min)
        upper_fill = np.where(max_open, PLOT_U2_MAX, u2_max)
        ax.fill_between(
            mass, lower_fill, upper_fill,
            alpha=0.25, color="red",
            label="GRENDEL" if not plotted else None, zorder=5)
        plotted = True

        lower_line = np.where(min_open, np.nan, u2_min)
        upper_line = np.where(max_open, np.nan, u2_max)
        ax.plot(mass, lower_line, "r-", linewidth=1.8, zorder=6)
        ax.plot(mass, upper_line, "r-", linewidth=1.8, zorder=6)

        # Close the red outline where the island physically shuts: join the lower
        # and upper edges with a vertical cap at any segment end where both edges
        # are sensitive (not open). This matches the filled region (whose polygon
        # already closes there) -- no extrapolation, just the real closing edge.
        for end in (0, -1):
            if not min_open[end] and not max_open[end]:
                ax.plot([mass[end], mass[end]], [u2_min[end], u2_max[end]],
                        "r-", linewidth=1.8, zorder=6)

        min_mark = _thin_marker_mask(min_open)
        if np.any(min_mark):
            ax.scatter(
                mass[min_mark], np.full(min_mark.sum(), PLOT_U2_MIN),
                marker="v", s=18, facecolors="none", edgecolors="red",
                linewidths=0.8, clip_on=False, zorder=7)
        max_mark = _thin_marker_mask(max_open)
        if np.any(max_mark):
            ax.scatter(
                mass[max_mark], np.full(max_mark.sum(), PLOT_U2_MAX),
                marker="^", s=18, facecolors="none", edgecolors="red",
                linewidths=0.8, clip_on=False, zorder=7)

    ax.set_xlabel(r"$m_N$ [GeV]", fontsize=14)
    if is_leftmost:
        ax.set_ylabel(r"$|U|^2$", fontsize=14)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim([0.15, 10.0])
    ax.set_ylim([PLOT_U2_MIN, PLOT_U2_MAX])
    ax.grid(True, which="both", alpha=0.2, linewidth=0.5)
    if plotted:
        ax.legend(fontsize=10, loc="upper right")
    ax.set_title(f"HNL {_FLAVOR_LABEL.get(flavor, flavor)}", fontsize=14)


def plot_nsignal_vs_u2(u2_grid, N_grid, mass, flavor, output_dir):
    """Diagnostic: N_signal vs U^2 for a single mass point."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(u2_grid, N_grid, "b-", linewidth=1.5)
    ax.axhline(3.0, color="r", ls="--", alpha=0.7, label=r"$N_\mathrm{thr} = 3$")
    ax.set_xlabel(r"$U^2$", fontsize=13)
    ax.set_ylabel(r"$N_\mathrm{signal}$", fontsize=13)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_title(f"{flavor}, $m_N$ = {mass:.2f} GeV", fontsize=13)
    ax.legend(fontsize=11)
    ax.grid(True, which="both", alpha=0.2)

    out = output_dir / f"nsignal_{flavor}_mN_{mass:.2f}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
