"""
2D map of inner separation (sep_in) vs pointing angle for the two-body signal,
to judge whether a separation-dependent (sloped) pointing cut makes sense
instead of the current hard gate (pointing < 50 mrad applied only for
sep_in < 10 cm).

Reuses the signal reconstruction (decayProbPerEvent_2body.sample_separations,
which calls reco_common) so the observables are defined identically to the
analysis. Overlays:
  - the weighted signal pointing envelope per sep_in bin (high percentile),
  - the current gate (sep_in < 10 cm -> pointing < 50 mrad),
  - the cosmic decay-in-flight survivors that leaked through the gate crack.

Usage:  python plot_sep_vs_pointing.py
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import decayProbPerEvent_2body as sig
from grendel_geometry import mesh_fiducial

LIFETIME_S = 1e-7          # 100 ns; pointing/sep geometry is ~lifetime-insensitive
N_PER = 200
COSMIC_CSV = 'cosmic_decay_collin_passers_empirical.csv'
ALL_MASSES = {'0.5': ('0.5 GeV', 'LLP0p5GeVSmall.csv'),
              '1':   ('1 GeV',   'LLP1GeV.csv'),
              '15':  ('15 GeV',  'LLPSmall.csv'),
              '40':  ('40 GeV',  'LLP40GeVSmall.csv')}
import sys
_keys = sys.argv[1:] if len(sys.argv) > 1 else ['0.5', '15']
MASSES = [ALL_MASSES[k] for k in _keys]


def weighted_quantile(x, w, q):
    if len(x) == 0 or w.sum() <= 0:
        return np.nan
    i = np.argsort(x)
    x, w = x[i], w[i]
    c = np.cumsum(w)
    return np.interp(q * c[-1], c, x)


def signal_observables(csv):
    """Return sep_in [cm], pointing [mrad], weights for decays passing the full
    baseline selection EXCEPT the pointing cut (so the natural sep-vs-pointing
    distribution is exposed)."""
    geo = sig.cache_geometry(csv, mesh_fiducial, [0, 0, 0])
    mc = sig.sample_separations(geo, LIFETIME_S, n_samples_per_particle=N_PER)
    sep, sepo = mc['sep'], mc['sep_outer']
    w = mc['weights']
    L = sig.DETECTOR_THICKNESS
    m = mc['on_tracker'].copy()
    m &= mc['p_soft'] >= sig.P_CUT
    m &= (sep >= sig.SEP_MIN) & (sepo >= sig.SEP_MIN) & (sep <= sig.SEP_MAX)
    m &= mc['dca'] <= sig.DCA_CUT
    par = mc['open_angle'] < sig.THETA_PARALLEL
    m &= (~par) | (sepo < sig.SEP_OUT_MAX_PARALLEL)
    gated = sepo > L
    m &= (~gated) | (mc['collin'] > sig.COLLIN_FRAC * L)
    m &= mc['vtx_in']
    return sep[m] * 100, mc['pointing'][m] * 1000, w[m]


def cosmic_survivors():
    """sep_in [cm], pointing [mrad] for cosmic collin-passers; split survivors."""
    if not os.path.exists(COSMIC_CSV):
        return None
    import pandas as pd
    d = pd.read_csv(COSMIC_CSV)
    return d


def main():
    cos = cosmic_survivors()
    fig, axes = plt.subplots(1, len(MASSES), figsize=(7.5 * len(MASSES), 6),
                             sharey=True, squeeze=False)
    axes = axes[0]
    ybins = np.logspace(-1, np.log10(3000), 70)   # 0.1 .. 3000 mrad
    xbins = np.linspace(0, 60, 61)                  # 0 .. 60 cm

    for ax, (label, csv) in zip(axes, MASSES):
        sx, py, w = signal_observables(csv)
        h = ax.hist2d(sx, py, bins=[xbins, ybins], weights=w,
                      norm=LogNorm(), cmap='viridis')
        fig.colorbar(h[3], ax=ax, label='signal yield (a.u.)')

        # signal pointing envelope per sep_in bin (weighted 99th percentile)
        xc = 0.5 * (xbins[:-1] + xbins[1:])
        p99, p50 = [], []
        for lo, hi in zip(xbins[:-1], xbins[1:]):
            sel = (sx >= lo) & (sx < hi)
            p99.append(weighted_quantile(py[sel], w[sel], 0.99))
            p50.append(weighted_quantile(py[sel], w[sel], 0.50))
        ax.plot(xc, p99, 'w-', lw=2, label='signal 99% envelope')
        ax.plot(xc, p50, 'w--', lw=1.2, alpha=0.8, label='signal median')

        # current gate: pointing < 50 mrad, only for sep_in < 10 cm
        ax.plot([0, 10], [50, 50], 'r-', lw=2.5, label='current cut (50 mrad)')
        ax.plot([10, 10], [50, 3000], 'r-', lw=2.5)
        ax.axvline(10, color='red', ls=':', lw=1, alpha=0.5)

        # cosmic survivors / collin-passers
        if cos is not None:
            surv = cos[cos['pass_all']]
            other = cos[~cos['pass_all']]
            ax.scatter(other['sep_in_cm'], other['pointing_mrad'], s=12,
                       c='lightgray', edgecolor='k', linewidth=0.2, alpha=0.6,
                       label='cosmic collin-passers')
            ax.scatter(surv['sep_in_cm'], surv['pointing_mrad'], s=120,
                       marker='*', c='orange', edgecolor='k',
                       label='cosmic survivors (final)')

        ax.set_yscale('log')
        ax.set_xlabel('inner separation sep_in (cm)')
        ax.set_title(f'{label} signal: sep_in vs pointing  '
                     f'(τ = {LIFETIME_S*1e9:.0f} ns)')
        ax.set_xlim(0, 60)
        ax.set_ylim(0.1, 3000)
        ax.legend(fontsize=8, loc='lower right')
    axes[0].set_ylabel('pointing angle (mrad)')
    plt.tight_layout()
    out = 'sep_vs_pointing_' + '_'.join(_keys).replace('.', 'p') + 'GeV.png'
    plt.savefig(out, dpi=150)
    print(f'wrote {out}')


if __name__ == '__main__':
    main()
