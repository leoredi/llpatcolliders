"""
Cosmic-muon transit MC and diagnostic plot.

Models the background where an atmospheric muon enters the tunnel from
above with the standard cos²θ_zenith surface angular distribution and
traverses the fiducial volume, leaving 4 collinear hits. Same
reconstruction and cutflow as for IP muons; this script only changes
how the rays are sampled.

Sampling:
  - Origins: uniform on a horizontal plane at altitude Y_THROW above the
    tunnel centreline, with the (X, Z) box padded by Y_THROW * tan(θ_max)
    so any direction within the cos²θ acceptance can still hit the tunnel.
  - Directions: zenith angle drawn from f(θ) = 3 cos²θ sinθ on [0, π/2],
    azimuth uniform on [0, 2π]. Sign chosen so the muon goes downward
    (negative Y, since the tunnel is at +Y = 22 m above the IP).

Outputs:
  - Console cutflow on the cosmic sample.
  - Figure cosmic_transit_observables.png mirroring the IP-muon plot.

Usage:
    python cosmic_muon_check.py [--n-rays 500000] [--y-throw 50] [--seed 42]
"""

import os
import sys
import argparse
import numpy as np
import matplotlib
import trimesh
from collections import defaultdict

from grendel_geometry import mesh_fiducial, DETECTOR_THICKNESS
from decayProbPerEvent_2body import (
    SEP_MIN, SEP_MAX, DCA_CUT,
    THETA_PARALLEL, SEP_OUT_MAX_PARALLEL,
    COLLIN_MIN, SEP_OUT_COLLIN_GATE,
    HIT_RESOLUTION,
)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--n-rays', type=int, default=500_000,
                   help='Number of cosmic rays to sample.')
    p.add_argument('--y-throw', type=float, default=50.0,
                   help='Y altitude (m) above the IP at which to throw cosmics. '
                        'The tunnel centreline sits at Y = 22 m; default puts the '
                        'throw plane 28 m above the tunnel.')
    p.add_argument('--theta-max', type=float, default=np.deg2rad(80),
                   help='Maximum zenith angle (rad) for sampling. Default 80°.')
    p.add_argument('--seed', type=int, default=42)
    p.add_argument('--interactive', action='store_true')
    p.add_argument('--out', default='cosmic_transit_observables.png')
    return p.parse_args()


def sample_cosmics(n_rays, y_throw, theta_max, seed):
    """Cosmic direction and origin sampling.

    Direction PDF: f(θ) = 3 cos²θ sinθ on [0, π/2]  (per unit solid angle ~ cos²θ).
    CDF: F(θ) = 1 − cos³θ  →  cosθ = (1 − U)^{1/3} for U ∈ [0, 1].
    """
    rng = np.random.default_rng(seed)
    u = rng.uniform(0, 1, n_rays)
    # Clamp to theta_max
    cos_min = np.cos(theta_max)
    # Truncated inverse CDF on [0, theta_max]: cos θ ∈ [cos_min, 1]
    cos_z = ((1 - u) * (1 - cos_min**3) + cos_min**3) ** (1.0/3.0)
    sin_z = np.sqrt(1 - cos_z**2)
    phi   = rng.uniform(0, 2*np.pi, n_rays)
    # CMS convention: Y up. Downgoing: dy = -cos_z
    dx = sin_z * np.cos(phi)
    dy = -cos_z
    dz = sin_z * np.sin(phi)
    dirs = np.column_stack([dx, dy, dz])

    # Origin sampling: a horizontal plane at Y = y_throw covering the tunnel
    # plus padding so steep rays don't miss it.
    bbox = mesh_fiducial.bounds  # (2, 3): [min, max] in (x, y, z)
    # Pad in X and Z by y_throw * tan(theta_max)
    pad = (y_throw - bbox[1, 1]) * np.tan(theta_max) + 5.0  # extra margin
    x_min, x_max = bbox[0, 0] - pad, bbox[1, 0] + pad
    z_min, z_max = bbox[0, 2] - pad, bbox[1, 2] + pad
    origins = np.column_stack([
        rng.uniform(x_min, x_max, n_rays),
        np.full(n_rays, y_throw),
        rng.uniform(z_min, z_max, n_rays),
    ])
    return origins, dirs


def run_cosmic_mc(n_rays, y_throw, theta_max, seed, sigma_hit, L):
    origins, dirs = sample_cosmics(n_rays, y_throw, theta_max, seed)
    rng = np.random.default_rng(seed + 1)

    intersector = trimesh.ray.ray_triangle.RayMeshIntersector(mesh_fiducial)
    print(f"Ray-casting {n_rays} cosmic rays against fiducial mesh...")
    locs, ray_idx, _ = intersector.intersects_location(
        origins, dirs, multiple_hits=True)
    ray_hits = defaultdict(list)
    for loc, r in zip(locs, ray_idx):
        # distance from origin along the (unit) direction
        d = np.dot(loc - origins[r], dirs[r])
        if d > 0:
            ray_hits[int(r)].append((d, loc))
    valid = [(r, dirs[r], origins[r],
              min(h[0] for h in v), max(h[0] for h in v))
             for r, v in ray_hits.items() if len(v) >= 2]
    nv = len(valid)
    print(f"  {nv} / {n_rays} cosmics traverse the fiducial "
          f"({nv/n_rays*100:.4f}% acceptance × throw-box).")
    if nv == 0:
        return None

    mu_dir = np.array([v[1] for v in valid])
    mu_ori = np.array([v[2] for v in valid])
    d_in   = np.array([v[3] for v in valid])
    d_out  = np.array([v[4] for v in valid])
    chord  = d_out - d_in
    # Sanity: also keep zenith for diagnostic
    cos_z  = -mu_dir[:, 1]   # downgoing → dy negative
    zenith = np.arccos(np.clip(cos_z, -1, 1))

    # 4 true hit positions along each cosmic's straight line
    P_eo = mu_ori + (d_in  - L)[:, None] * mu_dir
    P_ei = mu_ori + (d_in     )[:, None] * mu_dir
    P_xi = mu_ori + (d_out    )[:, None] * mu_dir
    P_xo = mu_ori + (d_out + L)[:, None] * mu_dir

    # Per-cosmic transverse basis for hit smearing
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

    d_t1 = H_eo - H_ei
    d_t2 = H_xo - H_xi
    mag1 = np.linalg.norm(d_t1, axis=1)
    mag2 = np.linalg.norm(d_t2, axis=1)
    sep_in  = np.linalg.norm(H_ei - H_xi, axis=1)
    sep_out = np.linalg.norm(H_eo - H_xo, axis=1)
    open_a  = np.arccos(np.clip((d_t1*d_t2).sum(axis=1)/(mag1*mag2), -1, 1))

    nvec = np.cross(d_t1, d_t2)
    nmag = np.linalg.norm(nvec, axis=1)
    w    = H_ei - H_xi
    dca  = np.where(nmag > 1e-20,
                    np.abs((w*nvec).sum(axis=1)) / nmag,
                    np.linalg.norm(np.cross(w, d_t1/mag1[:, None]), axis=1))

    pts = np.stack([H_eo, H_ei, H_xi, H_xo], axis=1)
    ctr = pts.mean(axis=1, keepdims=True)
    _, S, _ = np.linalg.svd(pts - ctr, full_matrices=False)
    collin = np.sqrt((S[:, 1]**2 + S[:, 2]**2) / 4)

    return dict(nv=nv, chord=chord, sep_in=sep_in, sep_out=sep_out,
                open_a=open_a, dca=dca, collin=collin, zenith=zenith)


def make_plot(res, out_path, interactive):
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # collinearity (linear)
    ax = axes[0, 0]
    co_mm = res['collin'] * 1000
    bins = np.linspace(0, max(15, COLLIN_MIN*1000*1.5), 80)
    ax.hist(co_mm, bins=bins, color='darkgreen', edgecolor='black',
            linewidth=0.3, alpha=0.85)
    ax.axvline(COLLIN_MIN*1000, color='red', linestyle='--', linewidth=2,
               label=f'cut = {COLLIN_MIN*1000:.0f} mm')
    ax.axvline(np.percentile(co_mm, 99), color='gray', linestyle=':',
               linewidth=1.5,
               label=f'99% = {np.percentile(co_mm, 99):.2f} mm')
    ax.set_xlabel('collinearity (mm)')
    ax.set_ylabel('Counts')
    ax.set_title(f'Cosmic transit collinearity '
                 f'(N = {res["nv"]}, σ_hit = {HIT_RESOLUTION*1000:.1f} mm)')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # collinearity log-log
    ax = axes[0, 1]
    positive = co_mm[co_mm > 0]
    lo = np.log10(max(positive.min(), 0.05)) if len(positive) else -1
    hi = np.log10(max(co_mm.max(), 100))
    bins_log = np.logspace(lo, hi, 80)
    ax.hist(co_mm, bins=bins_log, color='darkgreen', edgecolor='black',
            linewidth=0.3, alpha=0.85)
    ax.axvline(COLLIN_MIN*1000, color='red', linestyle='--', linewidth=2,
               label=f'cut = {COLLIN_MIN*1000:.0f} mm')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel('collinearity (mm)')
    ax.set_ylabel('Counts')
    ax.set_title(f'Log-log (max observed = {co_mm.max():.2f} mm)')
    ax.legend(fontsize=9)
    ax.grid(True, which='both', alpha=0.3)

    # sep_outer vs collinearity
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
    ax.set_xlabel('sep_outer (cm)'); ax.set_ylabel('collinearity (mm)')
    ax.set_title('Conditional cut region')
    ax.legend(fontsize=9, loc='lower right')
    plt.colorbar(h[3], ax=ax, label='Counts')

    # zenith angle and chord length (cosmic-specific diagnostic)
    ax = axes[1, 1]
    h = ax.hist2d(np.rad2deg(res['zenith']), res['chord'],
                  bins=[np.linspace(0, 90, 50),
                        np.linspace(0, max(res['chord'].max(), 5), 50)],
                  cmap='magma', cmin=1)
    ax.set_xlabel('zenith angle (deg)')
    ax.set_ylabel('fiducial chord (m)')
    ax.set_title('Cosmic geometry: zenith vs chord length')
    plt.colorbar(h[3], ax=ax, label='Counts')

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    print(f"Saved {out_path}")
    if interactive:
        plt.show()
    else:
        plt.close('all')


def cutflow_summary(res):
    sep_in, sep_out = res['sep_in'], res['sep_out']
    open_a, dca, collin = res['open_a'], res['dca'], res['collin']
    nv = res['nv']
    print(f"\nCutflow on cosmic transits (N = {nv}):")
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
    print("Cosmic observable distributions:")
    def pct(x): return tuple(round(v, 4) for v in np.percentile(x, [50, 90, 99]))
    print(f"  zenith       [deg] med/90/99 = {pct(np.rad2deg(res['zenith']))}")
    print(f"  chord         [m] med/90/99 = {pct(res['chord'])}")
    print(f"  sep_inner     [m] med/90/99 = {pct(sep_in)}")
    print(f"  sep_outer     [m] med/90/99 = {pct(sep_out)}")
    print(f"  open_angle  [rad] med/90/99 = {pct(open_a)}")
    print(f"  DCA           [m] med/90/99 = {pct(dca)}")
    print(f"  collinearity  [m] med/90/99 = {pct(collin)}, "
          f"max = {res['collin'].max():.4f}")


def main():
    args = parse_args()
    if not args.interactive and not os.environ.get('MPLBACKEND'):
        matplotlib.use('Agg')

    res = run_cosmic_mc(args.n_rays, args.y_throw, args.theta_max,
                        args.seed, HIT_RESOLUTION, DETECTOR_THICKNESS)
    if res is None:
        print("No cosmics traversed the fiducial.")
        return

    cutflow_summary(res)
    make_plot(res, args.out, args.interactive)


if __name__ == '__main__':
    main()
