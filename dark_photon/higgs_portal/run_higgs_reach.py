"""GRENDEL reach for h -> A' A' (long-lived dark photon) in the BR-vs-ctau plane.

Reuses the collaboration's Higgs-portal machinery (higgs/decayProbPerEvent_2body.py)
unchanged: it ray-casts each A' from the CMS IP into the PX56 tunnel fiducial
volume, MC-samples the two-body decay + full GRENDEL reconstruction selection
once, and reweights that single pass to every proper lifetime. The excluded
branching ratio at each ctau is

    BR(h->A'A')_excl = 3 / ( N_h * P(>=1 A' decays & passes) ),
    N_h = 52 pb (sigma_h) * 3000 fb^-1 (HL-LHC) = 1.56e8 Higgs,

i.e. the y-axis of the PBC h->A'A' plot. This is the identical normalization the
dark-scalar (h->SS) curves use; only the signal kinematics (m_A' = 2 GeV, ggF
Higgs) and the interpretation (dark photon) differ.

Usage:
    python run_higgs_reach.py --csv tmp/llp_hAA_2GeV.csv --mass 2.0
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
_REPO_ROOT = _HERE.parents[1]          # .../llpatcolliders_BC1
_HIGGS = _REPO_ROOT / "higgs"
for _p in (str(_HIGGS),):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import decayProbPerEvent_2body as sig          # noqa: E402
from grendel_geometry import mesh_fiducial       # noqa: E402

# ctau range to match the PBC h->A'A' figure: ctau in [1e-3, 1e3] m.
# tau = ctau / c, so tau in ~[3.3e-12, 3.3e-6] s.
N_TAU = 45
LIFETIMES = np.logspace(-11.5, -5.5, N_TAU)     # seconds
N_SAMPLES_PER_PARTICLE = 200


def run(csv, mass, out_dir):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    geo = sig.cache_geometry(csv, mesh_fiducial, [0, 0, 0])
    n_hit = int(np.sum(geo["hits"]))
    print(f"m_A' = {mass} GeV | {len(geo['hits'])} A', {n_hit} ray-cast into tunnel "
          f"({100*n_hit/len(geo['hits']):.1f}%)")

    # Acceptance-only curve (analytic; full-coverage tracking, kinematic cuts only).
    scan = sig.analyze_decay_vs_lifetime(csv, geo, LIFETIMES)
    # Headline curve: full reconstruction selection via one reweighted MC pass.
    mc = sig.sample_separations(geo, 1e-6, n_samples_per_particle=N_SAMPLES_PER_PARTICLE)
    mc_scan = sig.mc_exclusion_vs_lifetime(mc, LIFETIMES, scan["total_events"])

    ctau_m = LIFETIMES * sig.SPEED_OF_LIGHT
    br_full = np.asarray(mc_scan["exclusion"], float)
    br_acc = np.asarray(scan["exclusion"], float)

    best_i = int(np.nanargmin(br_full))
    print(f"  Best excluded BR (full selection): {br_full[best_i]:.3e} "
          f"at ctau = {ctau_m[best_i]:.3g} m")
    print(f"  Best excluded BR (acceptance only): {np.nanmin(br_acc):.3e}")

    # --- write curve CSV ---
    csv_out = out_dir / f"grendel_hAA_{_mass_tag(mass)}.csv"
    np.savetxt(
        csv_out,
        np.column_stack([ctau_m, br_full, br_acc]),
        delimiter=",", header="ctau_m,BR_full_selection,BR_acceptance_only",
        comments="",
    )
    print(f"  wrote {csv_out}")

    # --- plot (house style: PBC-matching log-log BR vs ctau) ---
    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    ax.loglog(ctau_m, br_full, color="crimson", lw=2.2,
              label="GRENDEL (full selection)")
    ax.loglog(ctau_m, br_acc, color="crimson", lw=1.8, ls="--", alpha=0.5,
              label="GRENDEL (acceptance only)")
    ax.set_xlim(1e-3, 1e3)
    ax.set_ylim(1e-4, 1e0)
    ax.set_xlabel(r"proper decay length $c\tau_{A'}$ (m)", fontsize=13)
    ax.set_ylabel(r"branching ratio BR($h\to A'A'$)", fontsize=13)
    ax.set_title(rf"Higgs decay to long-lived dark photons ($m_{{A'}} = {mass:g}$ GeV)",
                 fontsize=12)
    ax.grid(True, which="both", ls="-", alpha=0.2)
    ax.legend(fontsize=10, loc="lower right")
    fig.tight_layout()
    png = out_dir / f"grendel_hAA_{_mass_tag(mass)}.png"
    pdf = out_dir / f"grendel_hAA_{_mass_tag(mass)}.pdf"
    fig.savefig(png, dpi=150, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {png}\n  wrote {pdf}")
    return csv_out, png


def _mass_tag(mass):
    return f"{mass:g}".replace(".", "p") + "GeV"


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default=str(_HERE / "tmp" / "llp_hAA_2GeV.csv"))
    ap.add_argument("--mass", type=float, default=2.0)
    ap.add_argument("--out-dir", default=str(_HERE / "tmp"))
    args = ap.parse_args(argv)
    run(args.csv, args.mass, args.out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
