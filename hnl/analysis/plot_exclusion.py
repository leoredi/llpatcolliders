"""
hnl/analysis/plot_exclusion.py

Publication-quality (m_N, U^2) "money plot" for HNL sensitivity.

Adapted from llpatcolliders_FONLL/analysis/plot_exclusion.py:
- Panels are built only for flavors present in the results CSV (so a
  Umu-only scan produces a single-panel figure instead of 1x3 with two
  empty panels).
- Output filename is "hnl_exclusion" (project-neutral name).
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.signal import savgol_filter

from analysis.reference_curves import load_all_references

_REF_STYLE = {
    "MATHUSLA":  {"color": "#2166ac", "ls": "--",  "lw": 1.5},
    "ANUBIS":    {"color": "#4dac26", "ls": "-.",  "lw": 1.5},
    "CODEX-b":   {"color": "#e08214", "ls": ":",   "lw": 1.8},
    "SHiP":      {"color": "#7b3294", "ls": "--",  "lw": 1.2},
}

_FLAVOR_LABEL = {
    "Ue":   r"$|U_e|^2$",
    "Umu":  r"$|U_\mu|^2$",
    "Utau": r"$|U_\tau|^2$",
}


def plot_exclusion(results_csv, output_dir=None, basename="hnl_exclusion"):
    """Create the (m_N, U^2) exclusion figure from a sensitivity CSV."""
    results_csv = Path(results_csv)
    df = pd.read_csv(results_csv)
    if output_dir is None:
        output_dir = results_csv.parent
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    flavors_present = [f for f in ["Ue", "Umu", "Utau"]
                       if f in df["flavor"].unique()]
    if not flavors_present:
        print("No flavors in results CSV; nothing to plot.")
        return

    ref_curves = load_all_references(flavors=flavors_present)

    n = len(flavors_present)
    fig, axes = plt.subplots(1, n, figsize=(6 * n, 6), sharey=True,
                             squeeze=False)
    axes = axes[0]

    for idx, flavor in enumerate(flavors_present):
        _plot_single_panel(axes[idx], df, flavor, ref_curves,
                           is_leftmost=(idx == 0))

    fig.tight_layout(w_pad=2.0)
    out_pdf = output_dir / f"{basename}.pdf"
    out_png = output_dir / f"{basename}.png"
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_pdf}")
    print(f"Saved: {out_png}")


def _smooth_log(values, window=7):
    if len(values) < window:
        return values
    if window % 2 == 0:
        window += 1
    log_vals = np.log10(values)
    smoothed = savgol_filter(log_vals, window, polyorder=3)
    return 10.0 ** smoothed


def _build_island_polygon(sel):
    """Build a closed polygon for the GRENDEL exclusion island."""
    valid = sel[sel["has_sensitivity"] == True].copy()
    if len(valid) == 0:
        return None

    mass = valid["mass_GeV"].values
    u2_lo = _smooth_log(valid["u2_min"].values)
    u2_hi = _smooth_log(valid["u2_max"].values)

    all_sorted = sel.sort_values("mass_GeV")
    last_sens_mass = mass[-1]
    after = all_sorted[
        (all_sorted["mass_GeV"] > last_sens_mass) &
        (all_sorted["has_sensitivity"] == False)
    ]

    m_upper = list(mass)
    u2_upper = list(u2_hi)
    m_lower = list(mass)
    u2_lower = list(u2_lo)

    if len(after) > 0:
        row = after.iloc[0]
        peak_N_last = valid.iloc[-1]["peak_N"]
        peak_N_next = row["peak_N"]
        if peak_N_next < peak_N_last and peak_N_last > 0:
            frac = (3.0 - peak_N_last) / (peak_N_next - peak_N_last)
            frac = np.clip(frac, 0.01, 0.99)
            m_tip = last_sens_mass + frac * (row["mass_GeV"] - last_sens_mass)
        else:
            m_tip = last_sens_mass
        u2_tip = row["peak_u2"]
        m_upper.append(m_tip)
        u2_upper.append(u2_tip)
        m_lower.append(m_tip)
        u2_lower.append(u2_tip)
    else:
        m_upper.append(mass[-1])
        u2_upper.append(np.sqrt(u2_lo[-1] * u2_hi[-1]))
        m_lower.append(mass[-1])
        u2_lower.append(np.sqrt(u2_lo[-1] * u2_hi[-1]))

    first_sens_mass = mass[0]
    before = all_sorted[
        (all_sorted["mass_GeV"] < first_sens_mass) &
        (all_sorted["has_sensitivity"] == False)
    ]
    if len(before) > 0:
        row = before.iloc[-1]
        peak_N_first = valid.iloc[0]["peak_N"]
        peak_N_prev = row["peak_N"]
        if peak_N_prev < peak_N_first and peak_N_first > 0:
            frac = (3.0 - peak_N_first) / (peak_N_prev - peak_N_first)
            frac = np.clip(frac, 0.01, 0.99)
            m_tip = first_sens_mass + frac * (row["mass_GeV"] - first_sens_mass)
        else:
            m_tip = first_sens_mass
        u2_tip = row["peak_u2"]
        m_upper.insert(0, m_tip)
        u2_upper.insert(0, u2_tip)
        m_lower.insert(0, m_tip)
        u2_lower.insert(0, u2_tip)

    m_poly = np.array(m_upper + m_lower[::-1])
    u2_poly = np.array(u2_upper + u2_lower[::-1])
    return m_poly, u2_poly


def _plot_single_panel(ax, df, flavor, ref_curves, is_leftmost=True):
    sel = df[df["flavor"] == flavor].sort_values("mass_GeV")
    valid = sel[sel["has_sensitivity"] == True]

    if len(valid) > 0:
        poly = _build_island_polygon(sel)
        if poly is not None:
            m_poly, u2_poly = poly
            ax.fill(m_poly, u2_poly, alpha=0.25, color="red",
                    label="GRENDEL", zorder=5)
            ax.plot(m_poly, u2_poly, "r-", linewidth=1.8, zorder=6)

    for exp, curves in ref_curves.items():
        if flavor not in curves:
            continue
        c = curves[flavor]
        style = _REF_STYLE.get(exp, {"color": "gray", "ls": "-", "lw": 1.0})
        ax.fill_between(c["mass"], c["u2_min"], c["u2_max"],
                        alpha=0.08, color=style["color"])
        ax.plot(c["mass"], c["u2_min"],
                color=style["color"], ls=style["ls"], lw=style["lw"],
                label=exp)
        ax.plot(c["mass"], c["u2_max"],
                color=style["color"], ls=style["ls"], lw=style["lw"])

    ax.set_xlabel(r"$m_N$ [GeV]", fontsize=14)
    if is_leftmost:
        ax.set_ylabel(r"$|U|^2$", fontsize=14)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim([0.15, 10.0])
    ax.set_ylim([1e-12, 1e-1])
    ax.grid(True, which="both", alpha=0.2, linewidth=0.5)
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
