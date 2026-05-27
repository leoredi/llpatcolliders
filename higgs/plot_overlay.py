"""
Overlay MC distributions and cutflow tables from multiple signal points.

Usage:
    python plot_overlay.py mc_distributions_1GeV.npz mc_distributions_15GeV.npz ...

Reads .npz files produced by decayProbPerEvent_2body.py and generates:
  - overlay_separation.png    : separation distributions
  - overlay_pointing.png      : pointing angle distributions
  - overlay_cutflow.csv       : combined cutflow table
"""

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from decayProbPerEvent_2body import build_cutflow

COLORS = ['steelblue', 'darkorange', 'forestgreen', 'crimson',
          'mediumpurple', 'goldenrod', 'teal', 'hotpink']


def load_signals(file_list):
    """Load .npz files into a list of dicts."""
    signals = []
    for f in file_list:
        d = np.load(f, allow_pickle=True)
        signals.append({k: d[k] for k in d.files})
    return signals


def plot_separation_overlay(signals, output='overlay_separation.png'):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    for i, sig in enumerate(signals):
        label = str(sig['label'])
        color = COLORS[i % len(COLORS)]
        seps = sig['seps']
        weights = sig['weights']
        # Normalise so distributions are comparable
        w_norm = weights / weights.sum() if weights.sum() > 0 else weights

        bins_lin = np.linspace(0, 0.5, 80)
        ax1.hist(seps, bins=bins_lin, weights=w_norm, histtype='step',
                 color=color, linewidth=2, label=label)

        bins_log = np.logspace(-4, 1, 80)
        ax2.hist(seps, bins=bins_log, weights=w_norm, histtype='step',
                 color=color, linewidth=2, label=label)

    sep_min = float(signals[0].get('sep_min', 0.001))
    sep_max = float(signals[0].get('sep_max', 10.0))
    for ax in (ax1, ax2):
        ax.axvline(sep_min, color='red', linestyle='--', linewidth=1.5,
                   label=f'min sep = {sep_min*1000:.0f} mm')
        ax.axvline(sep_max, color='red', linestyle='-', linewidth=1.5,
                   label=f'max sep = {sep_max*100:.0f} cm')
        ax.set_ylabel('Normalised weighted counts')
        ax.legend(fontsize=8)

    ax1.set_xlabel('Separation at detector (m)')
    ax1.set_xlim(0, 0.5)
    ax1.set_title('Separation (linear)')

    ax2.set_xlabel('Separation at detector (m)')
    ax2.set_xscale('log')
    ax2.set_title('Separation (log)')

    plt.tight_layout()
    plt.savefig(output, dpi=150)
    plt.show()
    print(f"Saved {output}")


def plot_pointing_overlay(signals, output='overlay_pointing.png'):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    for i, sig in enumerate(signals):
        label = str(sig['label'])
        color = COLORS[i % len(COLORS)]
        seps = sig['seps']
        pointing = sig['pointing']
        weights = sig['weights']

        # Apply separation acceptance
        sep_min = float(sig.get('sep_min', 0.001))
        sep_max = float(sig.get('sep_max', 10.0))
        mask = (seps >= sep_min) & (seps <= sep_max)
        pt_mrad = pointing[mask] * 1000
        w = weights[mask]
        w_norm = w / w.sum() if w.sum() > 0 else w

        pct99 = np.percentile(pt_mrad, 99.5) if len(pt_mrad) > 0 else 10
        bins_lin = np.linspace(0, pct99, 80)
        ax1.hist(pt_mrad, bins=bins_lin, weights=w_norm, histtype='step',
                 color=color, linewidth=2, label=label)

        if len(pt_mrad) > 0 and pt_mrad.min() > 0:
            bins_log = np.logspace(np.log10(max(pt_mrad.min(), 1e-3)),
                                   np.log10(pt_mrad.max()), 80)
        else:
            bins_log = np.logspace(-3, 2, 80)
        ax2.hist(pt_mrad, bins=bins_log, weights=w_norm, histtype='step',
                 color=color, linewidth=2, label=label)

    ax1.set_xlabel('Pointing angle (mrad)')
    ax1.set_ylabel('Normalised weighted counts')
    ax1.set_title('Pointing angle (linear)')
    ax1.legend(fontsize=8)

    ax2.set_xlabel('Pointing angle (mrad)')
    ax2.set_ylabel('Normalised weighted counts')
    ax2.set_xscale('log')
    ax2.set_title('Pointing angle (log)')
    ax2.legend(fontsize=8)

    plt.tight_layout()
    plt.savefig(output, dpi=150)
    plt.show()
    print(f"Saved {output}")


def plot_dca_overlay(signals, output='overlay_dca.png'):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    for i, sig in enumerate(signals):
        label = str(sig['label'])
        color = COLORS[i % len(COLORS)]
        seps = sig['seps']
        dca = sig['dca']
        weights = sig['weights']

        sep_min = float(sig.get('sep_min', 0.001))
        sep_max = float(sig.get('sep_max', 10.0))
        mask = (seps >= sep_min) & (seps <= sep_max)
        dca_cm = dca[mask] * 100
        w = weights[mask]
        w_norm = w / w.sum() if w.sum() > 0 else w

        pct99 = np.percentile(dca_cm, 99.5) if len(dca_cm) > 0 else 10
        bins_lin = np.linspace(0, pct99, 80)
        ax1.hist(dca_cm, bins=bins_lin, weights=w_norm, histtype='step',
                 color=color, linewidth=2, label=label)

        pos = dca_cm[dca_cm > 0]
        if len(pos) > 0:
            bins_log = np.logspace(np.log10(max(pos.min(), 1e-4)),
                                   np.log10(pos.max()), 80)
        else:
            bins_log = np.logspace(-4, 2, 80)
        ax2.hist(dca_cm, bins=bins_log, weights=w_norm, histtype='step',
                 color=color, linewidth=2, label=label)

    dca_cut = float(signals[0].get('dca_cut', 0.01))
    for ax in (ax1, ax2):
        ax.axvline(dca_cut * 100, color='red', linestyle='--', linewidth=1.5,
                   label=f'DCA cut = {dca_cut*100:.1f} cm')
        ax.set_ylabel('Normalised weighted counts')
        ax.legend(fontsize=8)

    ax1.set_xlabel('DCA (cm)')
    ax1.set_title('DCA (linear)')
    ax2.set_xlabel('DCA (cm)')
    ax2.set_xscale('log')
    ax2.set_title('DCA (log)')

    plt.tight_layout()
    plt.savefig(output, dpi=150)
    plt.show()
    print(f"Saved {output}")


def build_combined_cutflow(signals, output='overlay_cutflow.csv'):
    rows = []
    for sig in signals:
        label = str(sig['label'])
        p_cut = float(sig.get('p_cut', 0.6))
        sep_min = float(sig.get('sep_min', 0.001))
        sep_max = float(sig.get('sep_max', 10.0))
        dca_cut = float(sig.get('dca_cut', 0.01))
        theta_parallel = float(sig.get('theta_parallel', 0.050))
        sep_out_max_parallel = float(sig.get('sep_out_max_parallel', 1.5))

        if 'open_angle' in sig.files:
            open_angle = sig['open_angle']
        else:
            open_angle = np.full(len(sig['seps']), np.inf)

        # Older npz files may lack sep_outer — fall back to "always passes"
        # both the lower-sep cut and the conditional max
        if 'sep_outer' in sig.files:
            sep_outer = sig['sep_outer']
        else:
            sep_outer = np.full(len(sig['seps']),
                                max(sep_min, 1e-3) * 10)

        cutflow = build_cutflow(sig['seps'], sig['pointing'], sig['weights'],
                                sig['momenta'], sig['p_soft'],
                                sig['dca'], sig['vtx_in'],
                                open_angle, sep_outer,
                                p_cut=p_cut, sep_min=sep_min, sep_max=sep_max,
                                dca_cut=dca_cut,
                                theta_parallel=theta_parallel,
                                sep_out_max_parallel=sep_out_max_parallel)
        for row in cutflow:
            row['signal'] = label
        rows.extend(cutflow)

    df = pd.DataFrame(rows)
    pivot = df.pivot(index='cut', columns='signal',
                     values='efficiency')
    cut_order = [r['cut'] for r in rows if r['signal'] == str(signals[0]['label'])]
    pivot = pivot.reindex(cut_order)

    print("\n" + "=" * 60)
    print("COMBINED CUTFLOW (cumulative efficiency)")
    print("=" * 60)
    print(pivot.to_string(float_format='%.4f'))
    print()

    df.to_csv(output, index=False)
    print(f"Saved {output}")
    return pivot


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    npz_files = sys.argv[1:]
    signals = load_signals(npz_files)

    print(f"Loaded {len(signals)} signal(s): "
          + ", ".join(str(s['label']) for s in signals))

    for s in signals:
        res = float(s.get('hit_resolution', 0))
        nl = int(s.get('n_layers', 0))
        if res > 0:
            print(f"  {s['label']}: hit_resolution={res*1000:.1f} mm, "
                  f"n_layers={nl}")

    plot_separation_overlay(signals)
    plot_pointing_overlay(signals)
    has_dca = all('dca' in s for s in signals)
    if has_dca:
        plot_dca_overlay(signals)
    build_combined_cutflow(signals)
