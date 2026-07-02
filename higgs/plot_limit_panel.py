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


def panel(csv, mass, out,mini=False):
    geo = sig.cache_geometry(csv, mesh_fiducial, [0, 0, 0])
    scan = sig.analyze_decay_vs_lifetime(csv, geo, LIFETIMES)
    mc = sig.sample_separations(geo, 1e-6, n_samples_per_particle=200)
    mc_scan = sig.mc_exclusion_vs_lifetime(mc, LIFETIMES, scan['total_events'])

    fig, ax = plt.subplots(figsize=(7, 6))
    ctm = LIFETIMES * sig.SPEED_OF_LIGHT
    if not mini:
        ax.loglog(ctm, mc_scan['exclusion'], color='blue', linewidth=2,
                  label='GRENDEL (full selection)')
        ax.loglog(ctm, scan['exclusion'], color='blue', linewidth=2,
                  linestyle='-', alpha=0.5, label='GRENDEL (acceptance only)')
    else:
        for i in range(len( mc_scan['exclusion'])):
            mc_scan['exclusion'][i] *= 100
        for i in range(len( mc_scan['exclusion'])):
            scan['exclusion'][i] *= 100
        ax.loglog(ctm, mc_scan['exclusion'], color='blue', linewidth=2,
                  label='GRENDEL (full selection)')
        ax.loglog(ctm, scan['exclusion'], color='blue', linewidth=2,
                  linestyle='-', alpha=0.5, label='GRENDEL (acceptance only)')

    def overlay(path, color, label, ls='-',scale=1):
        if os.path.exists(path):
            d = np.loadtxt(path, delimiter=',')
            ax.loglog(d[:, 0], d[:, 1]*scale, color=color, linewidth=2, ls=ls, label=label)
    if not mini:
        if mass == 15:
            overlay('external/MATHUSLA.csv', 'green', 'MATHUSLA','--')
            overlay('external/CODEX.csv', 'cyan', 'CODEX-b','--')
            # overlay('external/CMS.csv', 'purple', 'CMS')
            overlay('external/ANUBISPBC.csv', 'magenta', 'ANUBIS', '--')
            ax.set_ylim([1e-5, 1])
        elif mass == 0.5:
            overlay('external/CODEX0p5.csv', 'cyan', 'CODEX-b',"--")
    if mini:
        if mass == 15:
            overlay('external/CMS_current.csv', 'cyan', 'CMS current','-')
            overlay('external/ATLAS_current.csv', 'cyan', 'ATLAS current','-')

    ax.set_xlabel(r'$c\tau$ (m)')
    ax.set_ylabel('BR')
    ax.set_title(f'$m = {mass}$ GeV')
    ax.grid(True, which='both', ls='-', alpha=0.2)
    ax.legend(fontsize=9, loc='lower right')
    plt.tight_layout()
    plt.savefig(out, dpi=150)
    print('wrote', out, ' best excl BR =', f"{np.nanmin(mc_scan['exclusion']):.2e}")


# panel('LLP0p5GeVSmall.csv', 0.5, 'exclusion_panel_0p5GeV.png')
# panel('LLPSmall.csv', 15, 'exclusion_panel_15GeV.png')
panel('LLPSmall.csv', 15, 'exclusion_panel_15GeV_mini.png',True)
