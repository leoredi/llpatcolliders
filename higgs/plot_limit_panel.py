"""
Standalone exclusion-limit panel (the bottom-right panel of the full exclusion
figure) for given masses, with the current default selection. One figure per mass.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import decayProbPerEvent_2body as sig
from grendel_geometry import mesh_fiducial

LIFETIMES = np.logspace(-10.5, -3.5, 20)


def panel(csv, mass, out):
    geo = sig.cache_geometry(csv, mesh_fiducial, [0, 0, 0])
    scan = sig.analyze_decay_vs_lifetime(csv, geo, LIFETIMES)
    mc = sig.sample_separations(geo, 1e-6, n_samples_per_particle=200)
    mc_scan = sig.mc_exclusion_vs_lifetime(mc, LIFETIMES, scan['total_events'])

    fig, ax = plt.subplots(figsize=(7, 6))
    ctm = LIFETIMES * sig.SPEED_OF_LIGHT
    ax.loglog(ctm, mc_scan['exclusion'], color='blue', linewidth=2,
              label='GRENDEL (full selection)')
    ax.loglog(ctm, scan['exclusion'], color='blue', linewidth=2,
              linestyle='--', alpha=0.5, label='GRENDEL (acceptance only)')

    def overlay(path, color, label, ls='-'):
        if os.path.exists(path):
            d = np.loadtxt(path, delimiter=',')
            ax.loglog(d[:, 0], d[:, 1], color=color, linewidth=2, ls=ls, label=label)
    if mass == 15:
        overlay('external/MATHUSLA.csv', 'green', 'MATHUSLA')
        overlay('external/CODEX.csv', 'cyan', 'CODEX-b')
        overlay('external/CMS.csv', 'purple', 'CMS')
        overlay('external/ANUBISUpdateCons.csv', 'magenta', 'ANUBIS Cons', '--')
        ax.set_ylim([1e-5, 1])
    elif mass == 0.5:
        overlay('external/CODEX0p5.csv', 'cyan', 'CODEX-b')

    ax.set_xlabel(r'$c\tau$ (m)')
    ax.set_ylabel('BR')
    ax.set_title(f'$m = {mass}$ GeV')
    ax.grid(True, which='both', ls='-', alpha=0.2)
    ax.legend(fontsize=9, loc='lower right')
    plt.tight_layout()
    plt.savefig(out, dpi=150)
    print('wrote', out, ' best excl BR =', f"{np.nanmin(mc_scan['exclusion']):.2e}")


panel('LLP0p5GeVSmall.csv', 0.5, 'exclusion_panel_0p5GeV.png')
panel('LLPSmall.csv', 15, 'exclusion_panel_15GeV.png')
