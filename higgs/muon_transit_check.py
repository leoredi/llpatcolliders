"""
IP-muon-transit MC and diagnostic plot.

Models the background where a single muon shoots out from the CMS IP,
traverses the tunnel, and leaves 4 detector hits — two on the entry
wall and two on the exit wall. Reconstruction can pair these as
"entry-track" and "exit-track" and fake a vertex inside the fiducial
volume.

This script samples muon directions uniformly over the (η, φ) box that
covers the tunnel acceptance, finds the entry/exit points on the actual
fiducial mesh via ray-casting, places 4 hits with σ_hit smearing along
a transverse basis, and computes the same reconstruction observables
that `sample_separations` produces for signal:

  sep_inner, sep_outer, open_angle, DCA, collinearity

The output plot `muon_transit_observables.png` shows how each of these
distributes for muon transits, with the analysis cut values annotated.

Usage:
    python muon_transit_check.py [--n-dir 200000] [--seed 42] [--interactive]
"""

import os
import sys
import argparse
import numpy as np
import matplotlib
import trimesh
from collections import defaultdict

from grendel_geometry import (
    mesh_fiducial, DETECTOR_THICKNESS, eta_phi_to_direction,
)
from decayProbPerEvent_2body import (
    SEP_MIN, SEP_MAX, DCA_CUT,
    THETA_PARALLEL, SEP_OUT_MAX_PARALLEL,
    COLLIN_MIN, SEP_OUT_COLLIN_GATE,
    HIT_RESOLUTION,
)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--n-dir', type=int, default=200_000,
                   help='Number of muon directions to sample from the IP.')
    p.add_argument('--eta-max', type=float, default=4.0,
                   help='|η| bound for direction sampling.')
    p.add_argument('--seed', type=int, default=42)
    p.add_argument('--interactive', action='store_true',
                   help='Show plots instead of just saving them.')
    p.add_argument('--out', default='muon_transit_observables.png')
    return p.parse_args()


def run_muon_mc(n_dir, eta_max, seed, sigma_hit, L):
    """Sample muons, ray-cast, place hits, return reconstructed observables."""
    rng = np.random.default_rng(seed)
    eta = rng.uniform(-eta_max, eta_max, n_dir)
    phi = rng.uniform(-np.pi, np.pi, n_dir)
    dirs = np.array([eta_phi_to_direction(e, p) for e, p in zip(eta, phi)])

    intersector = trimesh.ray.ray_triangle.RayMeshIntersector(mesh_fiducial)
    print(f"Ray-casting {n_dir} directions against fiducial mesh...")
    locs, ray_idx, _ = intersector.intersects_location(
        np.zeros_like(dirs), dirs, multiple_hits=True)

    ray_hits = defaultdict(list)
    for loc, r in zip(locs, ray_idx):
        ray_hits[int(r)].append((np.linalg.norm(loc), loc))
    valid = [(r, dirs[r], min(h[0] for h in v), max(h[0] for h in v))
             for r, v in ray_hits.items() if len(v) >= 2]
    nv = len(valid)
    print(f"  {nv} / {n_dir} muons traverse the fiducial "
          f"({nv/n_dir*100:.2f}% acceptance).")

    if nv == 0:
        return None

    mu_dir = np.array([v[1] for v in valid])
    d_in   = np.array([v[2] for v in valid])
    d_out  = np.array([v[3] for v in valid])
    chord  = d_out - d_in

    # True 4 hit positions along each muon's straight-line trajectory
    P_eo = (d_in  - L)[:, None] * mu_dir
    P_ei = (d_in     )[:, None] * mu_dir
    P_xi = (d_out    )[:, None] * mu_dir
    P_xo = (d_out + L)[:, None] * mu_dir

    # Per-muon transverse basis for hit smearing
    helper = np.tile(np.array([0., 1., 0.]), (nv, 1))
    parallel = np.abs((mu_dir * helper).sum(axis=1)) > 0.95
    helper[parallel] = np.array([1., 0., 0.])
    e1 = np.cross(mu_dir, helper)
    e1 /= np.linalg.norm(e1, axis=1, keepdims=True)
    e2 = np.cross(mu_dir, e1)

    def smear(P):
        a = rng.normal(0, sigma_hit, (nv, 1))
        b = rng.normal(0, sigma_hit, (nv, 1))
        return P + a*e1 + b*e2

    H_eo, H_ei = smear(P_eo), smear(P_ei)
    H_xi, H_xo = smear(P_xi), smear(P_xo)

    # Track directions (inner-hit → outer-hit, mirroring sample_separations)
    d_t1 = H_eo - H_ei
    d_t2 = H_xo - H_xi
    mag1 = np.linalg.norm(d_t1, axis=1)
    mag2 = np.linalg.norm(d_t2, axis=1)

    sep_in  = np.linalg.norm(H_ei - H_xi, axis=1)
    sep_out = np.linalg.norm(H_eo - H_xo, axis=1)
    cos_op  = (d_t1 * d_t2).sum(axis=1) / (mag1 * mag2)
    open_a  = np.arccos(np.clip(cos_op, -1, 1))

    # DCA between the two reconstructed lines (in 3D)
    nvec = np.cross(d_t1, d_t2)
    nmag = np.linalg.norm(nvec, axis=1)
    w    = H_ei - H_xi
    dca  = np.where(nmag > 1e-20,
                    np.abs((w * nvec).sum(axis=1)) / nmag,
                    np.linalg.norm(np.cross(w, d_t1 / mag1[:, None]), axis=1))

    # Collinearity: SVD perpendicular residual of 4 hits to best-fit line
    pts = np.stack([H_eo, H_ei, H_xi, H_xo], axis=1)
    ctr = pts.mean(axis=1, keepdims=True)
    _, S, _ = np.linalg.svd(pts - ctr, full_matrices=False)
    collin = np.sqrt((S[:, 1]**2 + S[:, 2]**2) / 4)

    return dict(nv=nv, chord=chord, sep_in=sep_in, sep_out=sep_out,
                open_a=open_a, dca=dca, collin=collin)


def make_plot(res, out_path, interactive):
    """4-panel diagnostic figure for IP-muon-transit reco observables."""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # --- collinearity (main discriminator) ---
    ax = axes[0, 0]
    co_mm = res['collin'] * 1000
    bins_lin = np.linspace(0, max(15, COLLIN_MIN*1000*1.5), 80)
    ax.hist(co_mm, bins=bins_lin, color='steelblue', edgecolor='black',
            linewidth=0.3, alpha=0.85)
    ax.axvline(COLLIN_MIN*1000, color='red', linestyle='--', linewidth=2,
               label=f'cut = {COLLIN_MIN*1000:.0f} mm')
    pct99 = np.percentile(co_mm, 99)
    ax.axvline(pct99, color='gray', linestyle=':', linewidth=1.5,
               label=f'99%ile = {pct99:.2f} mm')
    ax.set_xlabel('collinearity (mm)')
    ax.set_ylabel('Counts')
    ax.set_title(f'IP-muon-transit collinearity '
                 f'(N = {res["nv"]}, σ_hit = {HIT_RESOLUTION*1000:.1f} mm)')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # --- collinearity log scale (shows the tail clearly) ---
    ax = axes[0, 1]
    bins_log = np.logspace(np.log10(max(co_mm.min(), 0.05)),
                           np.log10(max(co_mm.max(), 100)), 80)
    ax.hist(co_mm, bins=bins_log, color='steelblue', edgecolor='black',
            linewidth=0.3, alpha=0.85)
    ax.axvline(COLLIN_MIN*1000, color='red', linestyle='--', linewidth=2,
               label=f'cut = {COLLIN_MIN*1000:.0f} mm')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel('collinearity (mm)')
    ax.set_ylabel('Counts')
    ax.set_title(f'Collinearity log-log\n'
                 f'(max observed = {co_mm.max():.2f} mm)')
    ax.legend(fontsize=9)
    ax.grid(True, which='both', alpha=0.3)

    # --- sep_outer vs collinearity, with rejection region shaded ---
    ax = axes[1, 0]
    so_cm = res['sep_out'] * 100
    h = ax.hist2d(so_cm, co_mm,
                  bins=[np.logspace(0, np.log10(max(so_cm.max(), 1e3)), 50),
                        np.logspace(-1, np.log10(max(co_mm.max(), 100)), 50)],
                  cmap='viridis', cmin=1)
    ax.axvline(SEP_OUT_COLLIN_GATE*100, color='red', linestyle='--',
               linewidth=1.5,
               label=f'gate = {SEP_OUT_COLLIN_GATE*100:.0f} cm')
    ax.axhline(COLLIN_MIN*1000, color='red', linestyle='--', linewidth=1.5,
               label=f'cut = {COLLIN_MIN*1000:.0f} mm')
    ax.fill_between([SEP_OUT_COLLIN_GATE*100, 1e4], 1e-1, COLLIN_MIN*1000,
                    color='red', alpha=0.15, label='rejected')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel('sep_outer (cm)')
    ax.set_ylabel('collinearity (mm)')
    ax.set_title('Conditional cut region (gate × collinearity)')
    ax.legend(fontsize=9, loc='lower right')
    plt.colorbar(h[3], ax=ax, label='Counts')

    # --- companion observables: open_angle and DCA ---
    ax = axes[1, 1]
    op_mrad = res['open_a'] * 1000
    dca_cm  = res['dca'] * 100
    h = ax.hist2d(op_mrad, dca_cm,
                  bins=[np.linspace(0, max(op_mrad.max(), 3200), 50),
                        np.linspace(0, max(np.percentile(dca_cm, 99.5),
                                            DCA_CUT*100*1.5), 50)],
                  cmap='magma', cmin=1)
    ax.axhline(DCA_CUT*100, color='cyan', linestyle='--', linewidth=1.5,
               label=f'DCA cut = {DCA_CUT*100:.1f} cm')
    ax.axvline(THETA_PARALLEL*1000, color='lime', linestyle='--', linewidth=1.5,
               label=f'θ_parallel = {THETA_PARALLEL*1000:.0f} mrad')
    ax.set_xlabel('open_angle (mrad)')
    ax.set_ylabel('DCA (cm)')
    ax.set_title('open_angle vs DCA (cut values overlaid)')
    ax.legend(fontsize=9, loc='upper left')
    plt.colorbar(h[3], ax=ax, label='Counts')

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    print(f"Saved {out_path}")
    if interactive:
        plt.show()
    else:
        plt.close('all')


def cutflow_summary(res):
    """Print the survival of each cut step on the muon sample."""
    sep_in, sep_out = res['sep_in'], res['sep_out']
    open_a, dca, collin = res['open_a'], res['dca'], res['collin']
    nv = res['nv']
    print(f"\nCutflow on IP-muon transits (N = {nv}):")
    print(f"  {'cut':<38} {'survival':>14}")
    print("  " + "-"*52)
    m = np.ones(nv, dtype=bool)
    print(f"  {'reco':<38} {m.mean()*100:>10.4f}%")
    m &= (sep_in >= SEP_MIN) & (sep_out >= SEP_MIN)
    print(f"  {'sep_in & sep_out > SEP_MIN':<38} {m.mean()*100:>10.4f}%")
    m &= (sep_in <= SEP_MAX)
    print(f"  {'sep_in < SEP_MAX':<38} {m.mean()*100:>10.4f}%")
    m &= (dca <= DCA_CUT)
    print(f"  {'DCA < DCA_CUT':<38} {m.mean()*100:>10.4f}%")
    par = open_a < THETA_PARALLEL
    m &= (~par) | (sep_out < SEP_OUT_MAX_PARALLEL)
    print(f"  {'conditional max sep_outer':<38} {m.mean()*100:>10.4f}%")
    gated = sep_out > SEP_OUT_COLLIN_GATE
    m &= (~gated) | (collin > COLLIN_MIN)
    print(f"  {'collinearity > COLLIN_MIN (gated)':<38} {m.mean()*100:>10.4f}%")

    print()
    print("Observable distributions (IP-muon transits):")
    def pct(x): return tuple(round(v, 4) for v in np.percentile(x, [50, 90, 99]))
    print(f"  fiducial chord [m] med/90/99 = {pct(res['chord'])}")
    print(f"  sep_inner      [m] med/90/99 = {pct(sep_in)}")
    print(f"  sep_outer      [m] med/90/99 = {pct(sep_out)}")
    print(f"  open_angle    [rad] med/90/99 = {pct(open_a)}")
    print(f"  DCA           [m]  med/90/99 = {pct(dca)}")
    print(f"  collinearity  [m]  med/90/99 = {pct(collin)}, "
          f"max = {res['collin'].max():.4f}")


def main():
    args = parse_args()
    # Batch mode by default
    if not args.interactive and not os.environ.get('MPLBACKEND'):
        matplotlib.use('Agg')

    res = run_muon_mc(args.n_dir, args.eta_max, args.seed,
                      HIT_RESOLUTION, DETECTOR_THICKNESS)
    if res is None:
        print("No muons traversed the fiducial — nothing to plot.")
        return

    cutflow_summary(res)
    make_plot(res, args.out, args.interactive)


if __name__ == '__main__':
    main()
