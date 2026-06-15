"""
Pointing distribution of the cosmic decay-in-flight survivors (events passing
the FULL selection), to inform a separation-dependent pointing cut.

Left  : rate-weighted pointing histogram of the survivors.
Right : survivor pointing vs sep_in (rate-coloured), with the signal 99%
        pointing envelopes (0.5 / 15 / 40 GeV) and the current gate overlaid.

Usage:  python plot_survivor_pointing.py [survivors.csv]
"""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

import plot_sep_vs_pointing as P

CSV = sys.argv[1] if len(sys.argv) > 1 else \
    'cosmic_decay_collin_passers_empirical_p1_c45.csv'
SIG = [('0.5 GeV', 'LLP0p5GeVSmall.csv', 'tab:green'),
       ('15 GeV',  'LLPSmall.csv',       'tab:blue'),
       ('40 GeV',  'LLP40GeVSmall.csv',  'tab:purple')]


def main():
    d = pd.read_csv(CSV)
    s = d[d.pass_all]
    pt = s.pointing_mrad.values
    w = s.rate_hz.values
    print(f'{len(s)} survivors; rate {w.sum():.3e} Hz; '
          f'pointing median {np.median(pt):.0f} mrad, '
          f'min {pt.min():.0f}, max {pt.max():.0f}')

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(15, 6))

    # --- left: rate-weighted pointing histogram of survivors ---
    bins = np.logspace(1, np.log10(3000), 30)
    a1.hist(pt, bins=bins, weights=w * 2e8, color='orange',
            edgecolor='k', alpha=0.8)
    a1.axvline(np.median(pt), color='red', ls='--',
               label=f'median {np.median(pt):.0f} mrad')
    a1.set_xscale('log')
    a1.set_xlabel('pointing angle (mrad)')
    a1.set_ylabel('survivor events / 2e8 s')
    a1.set_title(f'Cosmic survivor pointing  ({len(s)} MC, '
                 f'{w.sum()*2e8:.0f} events/2e8 s)')
    a1.legend()

    # --- right: survivor pointing vs sep_in + signal envelopes ---
    sc = a2.scatter(s.sep_in_cm, pt, c=w, cmap='Oranges',
                    s=40, edgecolor='k', linewidth=0.3, zorder=5,
                    label='cosmic survivors')
    fig.colorbar(sc, ax=a2, label='rate weight [Hz]')

    xbins = np.linspace(0, 70, 36)
    xc = 0.5 * (xbins[:-1] + xbins[1:])
    for label, csv, col in SIG:
        sx, py, ws = P.signal_observables(csv)
        p99 = [P.weighted_quantile(py[(sx >= lo) & (sx < hi)],
                                   ws[(sx >= lo) & (sx < hi)], 0.99)
               for lo, hi in zip(xbins[:-1], xbins[1:])]
        a2.plot(xc, p99, color=col, lw=2, label=f'{label} signal 99%')

    # current gate
    a2.plot([0, 10], [50, 50], 'r-', lw=2.5, label='current cut (50 mrad)')
    a2.plot([10, 10], [50, 3000], 'r-', lw=2.5)

    a2.set_yscale('log')
    a2.set_xlabel('inner separation sep_in (cm)')
    a2.set_ylabel('pointing angle (mrad)')
    a2.set_title('Survivor pointing vs sep_in, with signal envelopes')
    a2.set_xlim(0, 70)
    a2.set_ylim(10, 3000)
    a2.legend(fontsize=8, loc='lower right')

    plt.tight_layout()
    out = 'survivor_pointing.png'
    plt.savefig(out, dpi=150)
    print(f'wrote {out}')


if __name__ == '__main__':
    main()
