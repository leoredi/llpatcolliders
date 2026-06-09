"""
Compare the signal selection computed with the OLD idealized local-frame
reconstruction (decayProbPerEvent_2body.sample_separations) vs the unified
real-3D reconstruction (reco_common.reconstruct_3d), on the same events.

Keeps both definitions side by side so the effect of the inconsistency we found
(idealized hits at z in {0,L} vs real 3D wall hits) can be quantified per cut.

    python compare_reco.py [--csv LLPSmall.csv] [--n 250] [--tau 1e-7]
"""
import argparse
import numpy as np
import decayProbPerEvent_2body as sig
from grendel_geometry import mesh_fiducial, DETECTOR_THICKNESS
import reco_common as rc

L = DETECTOR_THICKNESS


def cumulative_cutflow(obs, w, tag, collin_frac=None):
    """Weighted cumulative efficiency through the signal selection.
    obs: dict with sep, sep_outer, open_angle, dca, collin, pointing, vtx_in,
    on_tracker, p_soft.  collin_frac: if set, use targeted collin>frac*L,
    gate sep_out>L; else nominal COLLIN_MIN / SEP_OUT_COLLIN_GATE."""
    m = np.ones(len(w), bool)
    tot = w.sum()
    rows = []

    def step(name, cond):
        nonlocal m
        m = m & cond
        rows.append((name, w[m].sum() / tot))

    step('p_soft>P_CUT', obs['p_soft'] >= sig.P_CUT)
    step('on_tracker', obs['on_tracker'])
    step('sep>SEP_MIN', (obs['sep'] >= sig.SEP_MIN) & (obs['sep_outer'] >= sig.SEP_MIN))
    step('sep<SEP_MAX', obs['sep'] <= sig.SEP_MAX)
    step('dca<DCA_CUT', obs['dca'] <= sig.DCA_CUT)
    par = obs['open_angle'] < sig.THETA_PARALLEL
    step('parallel', (~par) | (obs['sep_outer'] < sig.SEP_OUT_MAX_PARALLEL))
    if collin_frac is None:
        gated = obs['sep_outer'] > sig.SEP_OUT_COLLIN_GATE
        step('collin(nominal 30mm)', (~gated) | (obs['collin'] > sig.COLLIN_MIN))
    else:
        gated = obs['sep_outer'] > L
        step(f'collin>{collin_frac:.2f}L', (~gated) | (obs['collin'] > collin_frac * L))
    step('vtx_in', obs['vtx_in'])
    gpt = obs['sep'] < sig.SEP_IN_POINT_GATE
    step('pointing', (~gpt) | (obs['pointing'] < sig.POINT_TIGHT_SEP_IN))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--csv', default='LLPSmall.csv')
    ap.add_argument('--n', type=int, default=250)
    ap.add_argument('--tau', type=float, default=1e-7)
    ap.add_argument('--collin-frac', type=float, default=0.48)
    args = ap.parse_args()

    geo = sig.cache_geometry(args.csv, mesh_fiducial, [0, 0, 0])
    # idealized observables: explicitly request the OLD local-frame reco
    mc = sig.sample_separations(geo, lifetime_seconds=args.tau,
                                n_samples_per_particle=args.n, use_3d_reco=False)
    N = len(mc['collin'])
    rng = np.random.default_rng(7)
    print(f"\nsignal {args.csv}: {N} samples, L={L*100:.0f} cm, tau={args.tau:.0e}s")

    # ---- real-3D hits; keep only events where all 4 hits are well-defined ----
    p1, p2, d1, d2 = mc['exit_pt_1'], mc['exit_pt_2'], mc['dir1'], mc['dir2']
    in1, out1 = rc.wall_inner_outer(p1, d1, L)
    in2, out2 = rc.wall_inner_outer(p2, d2, L)
    fin = (np.all(np.isfinite(in1), 1) & np.all(np.isfinite(out1), 1)
           & np.all(np.isfinite(in2), 1) & np.all(np.isfinite(out2), 1))
    print(f"  {fin.sum()}/{N} events have well-defined 3D hits "
          f"(compared on this common set)")

    # idealized observables (old reconstruction), restricted to the common set
    ideal = dict(sep=mc['sep'][fin], sep_outer=mc['sep_outer'][fin],
                 open_angle=mc['open_angle'][fin], dca=mc['dca'][fin],
                 collin=mc['collin'][fin], pointing=mc['pointing'][fin],
                 vtx_in=mc['vtx_in'][fin], on_tracker=mc['on_tracker'][fin],
                 p_soft=mc['p_soft'][fin])

    # unified real-3D observables (shared reconstruction) on the same events
    r3 = rc.reconstruct_3d(out1[fin], in1[fin], in2[fin], out2[fin],
                           sig.HIT_RESOLUTION, rng)
    three = dict(sep=r3['sep'], sep_outer=r3['sep_outer'],
                 open_angle=r3['open_angle'], dca=r3['dca'], collin=r3['collin'],
                 pointing=r3['pointing'], vtx_in=r3['vtx_in'],
                 on_tracker=mc['on_tracker'][fin], p_soft=mc['p_soft'][fin])

    w = mc['weights'][fin]
    for frac, label in ((None, 'NOMINAL collin (30mm)'),
                        (args.collin_frac, f'TARGETED collin>{args.collin_frac}L')):
        ri = cumulative_cutflow(ideal, w, 'ideal', frac)
        r3c = cumulative_cutflow(three, w, '3d', frac)
        print(f"\n=== {label} : cumulative signal efficiency ===")
        print(f"  {'cut':<22} {'idealized':>11} {'real-3D':>11}")
        for (nm, ei), (_, e3) in zip(ri, r3c):
            print(f"  {nm:<22} {ei*100:>10.2f}% {e3*100:>10.2f}%")


if __name__ == '__main__':
    main()
