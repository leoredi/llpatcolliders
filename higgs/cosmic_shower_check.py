"""
Cosmic-shower (2 parallel muons) transit MC.

Cosmic showers can deliver multiple muons into the detector within one
trigger window. With N=2 the reconstruction sees 8 hits and can form
several "vertex" candidates by pairing tracks across the two muons.
The single-cosmic collinearity cut catches the correct pairings (each
muon's own 4 hits are collinear) but the cross-pairings have 4 hits on
*two* parallel lines, separated by the muon-to-muon distance d.

This script:
  - Samples cosmic direction from cos²θ_zenith.
  - Places 2 parallel muons separated by Δr in the plane perpendicular
    to the shower direction; both required to traverse the fiducial.
  - Reconstructs all four pairings (2 correct + 1 anti-parallel cross
    + 1 same-side parallel) and applies the analysis cutflow to each.
  - Outputs:
       * Cutflow per pairing as a function of separation d.
       * Figure cosmic_shower_observables.png showing survival vs d.

The separation is swept on a log axis from 1 mm to 10 m so the cut
boundaries are visible. Two reference shower-LDF marks are overlaid so
the survival curves can be convolved with a rate model later.

Usage:
    python cosmic_shower_check.py [--n-per-bin 5000] [--seed 42]
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
    p.add_argument('--n-per-bin', type=int, default=5000,
                   help='Number of cosmic showers per separation bin.')
    p.add_argument('--n-bins', type=int, default=40,
                   help='Number of separation bins (log-spaced 1 mm – 10 m).')
    p.add_argument('--y-throw', type=float, default=50.0)
    p.add_argument('--theta-max', type=float, default=np.deg2rad(80))
    p.add_argument('--seed', type=int, default=42)
    p.add_argument('--interactive', action='store_true')
    p.add_argument('--out', default='cosmic_shower_observables.png')
    return p.parse_args()


def sample_cosmic_direction(n, theta_max, rng):
    """Sample cos²θ_zenith direction vectors (Y up, downgoing)."""
    cos_min = np.cos(theta_max)
    u = rng.uniform(0, 1, n)
    cos_z = ((1 - u) * (1 - cos_min**3) + cos_min**3) ** (1.0/3.0)
    sin_z = np.sqrt(1 - cos_z**2)
    phi = rng.uniform(0, 2*np.pi, n)
    return np.column_stack([sin_z*np.cos(phi), -cos_z, sin_z*np.sin(phi)])


def random_perp_basis(direction, rng):
    """Two orthonormal vectors perpendicular to `direction`."""
    n = direction.shape[0]
    helper = np.tile(np.array([0., 1., 0.]), (n, 1))
    parallel = np.abs((direction * helper).sum(axis=1)) > 0.95
    helper[parallel] = np.array([1., 0., 0.])
    e1 = np.cross(direction, helper); e1 /= np.linalg.norm(e1, axis=1, keepdims=True)
    e2 = np.cross(direction, e1)
    # Pick a random azimuth so the offset is isotropic in the perp plane
    psi = rng.uniform(0, 2*np.pi, n)
    return (np.cos(psi)[:, None]*e1 + np.sin(psi)[:, None]*e2)


def find_traversal(origins, dirs, intersector):
    """Ray-cast against the fiducial mesh; return entry/exit distances
    for rays that traverse (≥ 2 intersections)."""
    locs, ray_idx, _ = intersector.intersects_location(
        origins, dirs, multiple_hits=True)
    rh = defaultdict(list)
    for loc, r in zip(locs, ray_idx):
        d = np.dot(loc - origins[r], dirs[r])
        if d > 0:
            rh[int(r)].append(d)
    return {r: (min(v), max(v)) for r, v in rh.items() if len(v) >= 2}


def reco_observables(P_eo, P_ei, P_xi, P_xo, L):
    """Compute reco observables from 4 hit position arrays of shape (N, 3)."""
    d_t1 = P_eo - P_ei
    d_t2 = P_xo - P_xi
    mag1 = np.linalg.norm(d_t1, axis=1)
    mag2 = np.linalg.norm(d_t2, axis=1)
    sep_in  = np.linalg.norm(P_ei - P_xi, axis=1)
    sep_out = np.linalg.norm(P_eo - P_xo, axis=1)
    open_a  = np.arccos(np.clip((d_t1*d_t2).sum(1)/(mag1*mag2), -1, 1))
    nvec = np.cross(d_t1, d_t2); nmag = np.linalg.norm(nvec, axis=1)
    w = P_ei - P_xi
    dca = np.where(nmag > 1e-20,
                   np.abs((w*nvec).sum(1))/nmag,
                   np.linalg.norm(np.cross(w, d_t1/mag1[:, None]), axis=1))
    pts = np.stack([P_eo, P_ei, P_xi, P_xo], axis=1)
    ctr = pts.mean(axis=1, keepdims=True)
    _, S, _ = np.linalg.svd(pts - ctr, full_matrices=False)
    collin = np.sqrt((S[:, 1]**2 + S[:, 2]**2) / 4)
    return dict(sep_in=sep_in, sep_out=sep_out, open_angle=open_a,
                dca=dca, collinearity=collin)


def apply_cutflow(o):
    """Return boolean mask of events passing the full analysis cutflow."""
    m = (o['sep_in']  >= SEP_MIN) & (o['sep_out'] >= SEP_MIN)
    m &= (o['sep_in'] <= SEP_MAX)
    m &= (o['dca']    <= DCA_CUT)
    par = o['open_angle'] < THETA_PARALLEL
    m &= (~par) | (o['sep_out'] < SEP_OUT_MAX_PARALLEL)
    gated = o['sep_out'] > SEP_OUT_COLLIN_GATE
    m &= (~gated) | (o['collinearity'] > COLLIN_MIN)
    return m


def sample_muons_hitting_tunnel(n_target, theta_max, rng, intersector):
    """Generate `n_target` cosmic muons that are guaranteed to traverse the
    fiducial volume. Anchors a point inside the volume per muon, then back-
    projects to the throw plane along a cos²θ-sampled direction."""
    # Random points inside the fiducial volume (a guaranteed-hit anchor)
    pts_inside = trimesh.sample.volume_mesh(mesh_fiducial, n_target)
    if len(pts_inside) < n_target:
        # Retry until we have enough
        more = []
        while sum(len(x) for x in more) + len(pts_inside) < n_target:
            more.append(trimesh.sample.volume_mesh(mesh_fiducial, n_target))
        pts_inside = np.concatenate([pts_inside] + more)[:n_target]
    dirs = sample_cosmic_direction(n_target, theta_max, rng)
    # Choose origins way above the tunnel along each ray (doesn't matter where,
    # we just need a back-projected starting point). Use 100 m back along -dir.
    origins = pts_inside - 100.0 * dirs
    travA = find_traversal(origins, dirs, intersector)
    # Keep only those that actually traverse (should be ~all by construction)
    valid_idx = sorted(travA.keys())
    return (origins[valid_idx], dirs[valid_idx],
            np.array([travA[i] for i in valid_idx]))


def simulate_pair_at_separation(d, n, theta_max, sigma_hit, L, rng,
                                 intersector,
                                 muon_A_cache=None):
    """Simulate cosmic showers with 2 parallel muons separated by d.
    Muon A is anchored inside the fiducial; muon B is offset perpendicular
    by d and may or may not hit."""

    if muon_A_cache is None:
        ori_A_all, dirs_all, dA_all = sample_muons_hitting_tunnel(
            n, theta_max, rng, intersector)
    else:
        ori_A_all, dirs_all, dA_all = muon_A_cache

    n_A = len(ori_A_all)
    # Offset muon B: distance d in a random direction perpendicular to muon A
    perp = random_perp_basis(dirs_all, rng)
    origins_B = ori_A_all + d * perp

    travB = find_traversal(origins_B, dirs_all, intersector)
    common = sorted(travB.keys())
    if not common:
        return None, 0, n_A

    common = np.array(common, dtype=int)
    n_valid = len(common)

    mu_dir = dirs_all[common]
    ori_A  = ori_A_all[common]
    ori_B  = origins_B[common]
    dA = dA_all[common]
    dB = np.array([[travB[r][0], travB[r][1]] for r in common])

    # Per-muon transverse basis (independent for hit smearing)
    e1 = random_perp_basis(mu_dir, rng);
    e2 = np.cross(mu_dir, e1)
    e2 /= np.linalg.norm(e2, axis=1, keepdims=True)

    def smear(P):
        a = rng.normal(0, sigma_hit, (n_valid, 1))
        b = rng.normal(0, sigma_hit, (n_valid, 1))
        return P + a*e1 + b*e2

    # 4 true hit positions per muon
    A_eo = ori_A + (dA[:, 0]-L)[:, None]*mu_dir; A_eo = smear(A_eo)
    A_ei = ori_A + (dA[:, 0]   )[:, None]*mu_dir; A_ei = smear(A_ei)
    A_xi = ori_A + (dA[:, 1]   )[:, None]*mu_dir; A_xi = smear(A_xi)
    A_xo = ori_A + (dA[:, 1]+L)[:, None]*mu_dir; A_xo = smear(A_xo)
    B_eo = ori_B + (dB[:, 0]-L)[:, None]*mu_dir; B_eo = smear(B_eo)
    B_ei = ori_B + (dB[:, 0]   )[:, None]*mu_dir; B_ei = smear(B_ei)
    B_xi = ori_B + (dB[:, 1]   )[:, None]*mu_dir; B_xi = smear(B_xi)
    B_xo = ori_B + (dB[:, 1]+L)[:, None]*mu_dir; B_xo = smear(B_xo)

    # reco_observables takes (P_eo, P_ei, P_xi, P_xo) where eo,ei are the
    # outer/inner of track 1 and xi,xo are the inner/outer of track 2.
    # For each pairing we have to pick which physical hit plays each role.
    pairings = {
        'A self (correct)':  (A_eo, A_ei, A_xi, A_xo),
        'B self (correct)':  (B_eo, B_ei, B_xi, B_xo),
        # cross: track 1 = entry of one muon, track 2 = exit of the other →
        # the two reconstructed tracks point in opposite directions (open=π).
        'cross A_entry+B_exit': (A_eo, A_ei, B_xi, B_xo),
        'cross B_entry+A_exit': (B_eo, B_ei, A_xi, A_xo),
        # same-side entries: BOTH tracks are entry-side, both point inward
        # (outer→inner = anti-muon-dir) → parallel tracks, open ≈ 0.
        'same-side entries':    (A_eo, A_ei, B_ei, B_eo),
        'same-side exits':      (A_xo, A_xi, B_xi, B_xo),
    }
    survival = {}
    for name, (eo, ei, xi, xo) in pairings.items():
        o = reco_observables(eo, ei, xi, xo, L)
        passed = apply_cutflow(o)
        survival[name] = int(passed.sum())
    return survival, n_valid, n_A


def plot(separations, valid, surv, out_path, interactive):
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    d_cm = np.array(separations) * 100  # convert to cm for x-axis
    valid = np.array(valid)
    eff = {k: np.array(v) / np.maximum(valid, 1) for k, v in surv.items()}

    # Survival per pairing
    ax = axes[0]
    pairing_colors = {
        'A self (correct)':         'gray',
        'B self (correct)':         'lightgray',
        'cross A_entry+B_exit':     'crimson',
        'cross B_entry+A_exit':     'orange',
        'same-side entries':        'royalblue',
        'same-side exits':          'steelblue',
    }
    for name, e in eff.items():
        ax.plot(d_cm, e, lw=2, color=pairing_colors[name], label=name)
    ax.axvline(DCA_CUT*100, color='black', ls=':', alpha=0.6,
               label=f'DCA cut = {DCA_CUT*100:.0f} cm')
    ax.axvline(COLLIN_MIN*100*2, color='gray', ls=':', alpha=0.6,
               label=f'collin/2 = {COLLIN_MIN*100*2:.0f} cm')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel('muon-to-muon separation Δr (cm)')
    ax.set_ylabel('per-pairing survival probability')
    ax.set_title(f'Cosmic-shower 2-muon survival vs separation\n'
                 f'(per-bin N traversing both = {np.median(valid):.0f} median)')
    ax.set_ylim(max(1.0/valid.max() / 10, 1e-4), 2.0)
    ax.legend(fontsize=8, loc='lower right')
    ax.grid(True, which='both', alpha=0.3)

    # Total fake-vertex yield per 2-muon event (sum over wrong pairings)
    ax = axes[1]
    wrong_keys = ['cross A_entry+B_exit', 'cross B_entry+A_exit',
                  'same-side entries', 'same-side exits']
    total = sum(eff[k] for k in wrong_keys)
    ax.plot(d_cm, total, lw=2, color='crimson',
            label='wrong-pairing fakes summed')
    correct = eff['A self (correct)'] + eff['B self (correct)']
    ax.plot(d_cm, correct, lw=1.5, color='gray', alpha=0.7,
            label='correct pairings (should be 0)')
    ax.axvline(DCA_CUT*100, color='black', ls=':', alpha=0.6,
               label=f'DCA = {DCA_CUT*100:.0f} cm')
    ax.set_xscale('log')
    ax.set_xlabel('muon-to-muon separation Δr (cm)')
    ax.set_ylabel('expected fake-vertex yield per 2-muon event')
    ax.set_title('Fake-vertex yield per 2-muon event')
    ax.legend(fontsize=9, loc='upper right')
    ax.grid(True, which='both', alpha=0.3)
    # Mark the "danger zone" identified analytically
    ax.axvspan(0.1, 10, color='red', alpha=0.08,
               label='danger zone (1 mm – 10 cm)')

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    print(f"Saved {out_path}")
    if interactive:
        plt.show()
    else:
        plt.close('all')


def main():
    args = parse_args()
    if not args.interactive and not os.environ.get('MPLBACKEND'):
        matplotlib.use('Agg')

    rng = np.random.default_rng(args.seed)
    intersector = trimesh.ray.ray_triangle.RayMeshIntersector(mesh_fiducial)
    separations = np.logspace(-3, 1, args.n_bins)  # 1 mm – 10 m
    surv = defaultdict(list)
    valid_per_bin = []
    nA_per_bin   = []

    # Pre-sample muon A's that hit the tunnel (reuse across all d bins so the
    # per-pairing fraction = #pass / #(both hit) is a clean conditional rate)
    print(f"Pre-sampling {args.n_per_bin} muon-A's anchored inside fiducial...")
    muon_A_cache = sample_muons_hitting_tunnel(
        args.n_per_bin, args.theta_max, rng, intersector)

    print(f"Cosmic-shower 2-muon MC: {args.n_bins} bins × "
          f"{len(muon_A_cache[0])} muon-As, sweeping Δr ∈ [1 mm, 10 m]")
    for i, d in enumerate(separations):
        s, nv, n_A = simulate_pair_at_separation(
            d, args.n_per_bin, args.theta_max,
            HIT_RESOLUTION, DETECTOR_THICKNESS, rng, intersector,
            muon_A_cache=muon_A_cache)
        valid_per_bin.append(nv)
        nA_per_bin.append(n_A)
        if s is None:
            for k in ['A self (correct)', 'B self (correct)',
                      'cross A_entry+B_exit', 'cross B_entry+A_exit',
                      'same-side entries', 'same-side exits']:
                surv[k].append(0)
        else:
            for k, c in s.items():
                surv[k].append(c)
        if i % 5 == 0 or i == len(separations)-1:
            wrong = sum(surv[k][-1] for k in [
                'cross A_entry+B_exit', 'cross B_entry+A_exit',
                'same-side entries', 'same-side exits'])
            print(f"  Δr = {d*100:>8.3f} cm | both traverse: "
                  f"{nv:>5d}/{n_A} ({nv/n_A*100:>5.1f}%) | "
                  f"wrong-pair pass: {wrong:>5d}")

    plot(separations, valid_per_bin, surv, args.out, args.interactive)

    # Summary at a few sentinel separations
    print()
    print(f"Summary: per-2-muon-event probability of any wrong-pairing "
          f"vertex surviving the full cutflow:")
    print(f"  {'Δr':>10} {'cross':>10} {'same-side':>12}")
    for d, nv in zip(separations, valid_per_bin):
        if nv == 0: continue
        idx = list(separations).index(d)
        cross = (surv['cross A_entry+B_exit'][idx] +
                 surv['cross B_entry+A_exit'][idx]) / nv
        sames = (surv['same-side entries'][idx] +
                 surv['same-side exits'][idx]) / nv
        if d*100 in [0.1, 1, 3, 10, 30, 100, 300] or \
           any(abs(d*100 - t) / max(d*100, 1) < 0.1 for t in [0.1, 1, 3, 10, 30, 100, 300]):
            print(f"  {d*100:>8.2f}cm  {cross*100:>8.3f}%  {sames*100:>10.3f}%")


if __name__ == '__main__':
    main()
