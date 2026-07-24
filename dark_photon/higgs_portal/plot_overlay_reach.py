"""Overlay the GRENDEL h->A'A' (m_A' = 2 GeV) reach on the PBC reference curves.

Pure re-plot from saved CSVs (GRENDEL curve + higgs/external/*.csv). No ray-cast,
no recompute -> safe/fast. Reproduces the PBC 'Higgs decay to long-lived dark
photons' plane and shows where GRENDEL sits vs CODEX-b / MATHUSLA / ANUBIS and
the existing ATLAS/CMS displaced limits.

All external curves are (ctau_m, BR) just as decayProbPerEvent_2body.py loads them.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_HERE = Path(__file__).resolve().parent
_EXT = _HERE.parents[1] / "higgs" / "external"   # llpatcolliders_BC1/higgs/external

# Projected HL-LHC 95% CL upper limit on BR(h->invisible) ~ 2.5% (a hard cap on
# any exotic-Higgs BR regardless of lifetime).
H_INV_HLLHC = 0.025


def _load(name):
    d = np.loadtxt(_EXT / name, delimiter=",")
    return d[:, 0], d[:, 1]


def plot(csv, mass, out_png):
    d = np.loadtxt(csv, delimiter=",", skiprows=1)
    ctau, br_full = d[:, 0], d[:, 1]
    bi = int(np.nanargmin(br_full))

    fig, ax = plt.subplots(figsize=(7.6, 5.8))

    # --- existing constraints (digitized current LHC displaced searches) ---
    xlim = (1e-3, 1e3)
    grid = np.logspace(np.log10(xlim[0]), np.log10(xlim[1]), 400)
    env = np.full_like(grid, np.nan)
    for name in ("ATLAS_current.csv", "CMS_current.csv"):
        try:
            x, y = _load(name)
        except OSError:
            continue
        yi = np.interp(grid, x, y, left=np.nan, right=np.nan)
        env = np.fmin(env, yi) if np.any(np.isfinite(env)) else yi
    # existing constraints also cap at the h->inv HL-LHC plateau
    env = np.fmin(np.where(np.isfinite(env), env, np.inf), H_INV_HLLHC)
    ax.fill_between(grid, env, 1e0, color="0.55", alpha=0.35, lw=0, zorder=0)
    ax.text(3e-2, 3.5e-1, "excluded", color="0.35", fontsize=12, style="italic")

    # --- future comparators (loaded exactly like the machinery does) ---
    for name, color, lab in (
        ("MATHUSLA.csv", "green", "MATHUSLA"),
        ("CODEX.csv", "teal", "CODEX-b"),
        ("ANUBISPBC.csv", "magenta", "ANUBIS (PBC)"),
    ):
        try:
            x, y = _load(name)
            ax.loglog(x, y, color=color, lw=2, ls=":", label=lab, zorder=3)
        except OSError:
            pass

    # --- h -> inv HL-LHC line ---
    ax.axhline(H_INV_HLLHC, color="steelblue", lw=1.8, zorder=2)
    ax.text(1.3e-3, H_INV_HLLHC * 1.15, r"$h\to$inv (HL-LHC)",
            color="steelblue", fontsize=9)

    # --- GRENDEL ---
    ax.loglog(ctau, br_full, color="crimson", lw=2.6,
              label="GRENDEL (3 ab$^{-1}$)", zorder=5)
    ax.plot(ctau[bi], br_full[bi], "o", color="crimson", ms=6, zorder=6)
    ax.annotate(rf"best BR $\approx$ {br_full[bi]:.1e}" + "\n"
                rf"at $c\tau \approx$ {ctau[bi]:.2g} m",
                xy=(ctau[bi], br_full[bi]), xytext=(1.4e-2, 1.6e-5),
                fontsize=9, color="crimson",
                arrowprops=dict(arrowstyle="->", color="crimson", lw=1))

    ax.set_xlim(*xlim)
    ax.set_ylim(1e-5, 1e0)
    ax.set_xlabel(r"proper decay length $c\tau_{A'}$ (m)", fontsize=13)
    ax.set_ylabel(r"branching ratio BR($h\to A'A'$)", fontsize=13)
    ax.set_title(rf"Higgs decay to long-lived dark photons ($m_{{A'}} = {mass:g}$ GeV)",
                 fontsize=12)
    ax.grid(True, which="both", ls="-", alpha=0.15)
    ax.legend(fontsize=9, loc="lower right", framealpha=0.9)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    fig.savefig(str(out_png).replace(".png", ".pdf"), bbox_inches="tight")
    print("wrote", out_png)
    print(f"GRENDEL best BR = {br_full[bi]:.3e} at ctau = {ctau[bi]:.3g} m")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default=str(_HERE / "tmp" / "grendel_hAA_2GeV.csv"))
    ap.add_argument("--mass", type=float, default=2.0)
    ap.add_argument("--out", default=str(_HERE / "tmp" / "grendel_hAA_2GeV_overlay.png"))
    args = ap.parse_args(argv)
    plot(args.csv, args.mass, Path(args.out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
