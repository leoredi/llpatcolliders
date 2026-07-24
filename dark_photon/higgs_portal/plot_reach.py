"""Re-plot the GRENDEL h->A'A' reach from the saved curve CSV (fast, no recompute)."""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_HERE = Path(__file__).resolve().parent


def plot(csv, mass, out_png):
    d = np.loadtxt(csv, delimiter=",", skiprows=1)
    ctau, br_full, br_acc = d[:, 0], d[:, 1], d[:, 2]
    bi = int(np.nanargmin(br_full))

    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    ax.loglog(ctau, br_full, color="crimson", lw=2.4,
              label="GRENDEL (full selection)")
    ax.loglog(ctau, br_acc, color="crimson", lw=1.8, ls="--", alpha=0.45,
              label="GRENDEL (acceptance only)")
    ax.plot(ctau[bi], br_full[bi], "o", color="crimson", ms=6, zorder=6)
    ax.annotate(rf"best BR $\approx$ {br_full[bi]:.1e}" + "\n"
                rf"at $c\tau \approx$ {ctau[bi]:.2g} m",
                xy=(ctau[bi], br_full[bi]), xytext=(1.6, 2.5e-4),
                fontsize=9, color="crimson",
                arrowprops=dict(arrowstyle="->", color="crimson", lw=1))

    ax.set_xlim(1e-3, 1e3)
    ax.set_ylim(1e-5, 1e0)
    ax.set_xlabel(r"proper decay length $c\tau_{A'}$ (m)", fontsize=13)
    ax.set_ylabel(r"branching ratio BR($h\to A'A'$)", fontsize=13)
    ax.set_title(rf"Higgs decay to long-lived dark photons ($m_{{A'}} = {mass:g}$ GeV)",
                 fontsize=12)
    ax.grid(True, which="both", ls="-", alpha=0.2)
    ax.legend(fontsize=10, loc="upper center")
    fig.tight_layout()
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    fig.savefig(str(out_png).replace(".png", ".pdf"), bbox_inches="tight")
    print("wrote", out_png)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default=str(_HERE / "tmp" / "grendel_hAA_2GeV.csv"))
    ap.add_argument("--mass", type=float, default=2.0)
    ap.add_argument("--out", default=str(_HERE / "tmp" / "grendel_hAA_2GeV.png"))
    args = ap.parse_args(argv)
    plot(args.csv, args.mass, Path(args.out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
