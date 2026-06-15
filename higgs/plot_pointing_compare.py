"""
Pointing comparison (cosmic vs signals) for events surviving everything EXCEPT
the pointing cut, under the corrected radial-spacing geometry. Shows whether
pointing still separates the cosmic decay-in-flight background from signal.

Cosmic: 500k volume events (p<1 GeV window) for a quick test.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import cosmic_decay_check as cdc
import decayProbPerEvent_2body as sig
from grendel_geometry import mesh_fiducial

L = sig.DETECTOR_THICKNESS

# ---------------- cosmic: mask through vertex stage (pre-pointing) ----------
cdc.GEN_PMIN, cdc.GEN_PMAX = 0.0, 1.0
flo = cdc.empirical_window_fraction(0.0, 1.0)
r = cdc.run_volume(500_000, np.deg2rad(80), 1.0, 'empirical', 70.0,
                   cdc.HIT_RESOLUTION, 42, sigma_t=1.5e-9,
                   target_muon_rate_hz=600.0, chunk=200_000, spectrum_weight=flo)
stages = cdc.cutflow_stages(r)
m_cos = next(mask for name, mask in stages if 'vertex' in name)   # pre-pointing
cos_pt = r['pointing'][m_cos] * 1000
cos_w = r['w'][m_cos]
print(f"cosmic pre-pointing survivors: {int(m_cos.sum())}")


# ---------------- signals: baseline+collin+vtx, NO pointing -----------------
def sig_pointing(csv):
    geo = sig.cache_geometry(csv, mesh_fiducial, [0, 0, 0])
    mc = sig.sample_separations(geo, 1e-7, n_samples_per_particle=200)
    sep, sepo = mc['sep'], mc['sep_outer']
    m = mc['on_tracker'] & (mc['p_soft'] >= sig.P_CUT)
    m &= (sep >= sig.SEP_MIN) & (sepo >= sig.SEP_MIN) & (sep <= sig.SEP_MAX)
    m &= mc['dca'] <= sig.DCA_CUT
    par = mc['open_angle'] < sig.THETA_PARALLEL
    m &= (~par) | (sepo < sig.SEP_OUT_MAX_PARALLEL)
    gated = sepo > sig.SEP_OUT_GATE
    m &= (~gated) | (mc['collin'] > sig.COLLIN_FRAC * L)
    m &= mc['vtx_in']
    return mc['pointing'][m] * 1000, mc['weights'][m]


SIG = [('0.5 GeV', 'LLP0p5GeVSmall.csv', 'tab:green'),
       ('1 GeV',   'LLP1GeVSmall.csv',   'tab:olive'),
       ('15 GeV',  'LLPSmall.csv',       'tab:blue'),
       ('40 GeV',  'LLP40GeVSmall.csv',  'tab:purple')]

bins = np.logspace(0, np.log10(3000), 55)
plt.figure(figsize=(9, 6))
plt.hist(cos_pt, bins=bins, weights=cos_w, density=True, histtype='step',
         lw=2.5, color='k', label='cosmic')
print(f"cosmic: median {np.median(cos_pt):.0f} mrad")
for lbl, csv, c in SIG:
    pt, w = sig_pointing(csv)
    plt.hist(pt, bins=bins, weights=w, density=True, histtype='step',
             lw=2, color=c, label=lbl)
    print(f"{lbl}: median {np.median(pt):.0f} mrad, 99% {np.percentile(pt,99):.0f}")

plt.axvline(sig.POINT_GLOBAL * 1000, color='red', ls='-', lw=1.5,
            label=f'global cut = {sig.POINT_GLOBAL*1000:.0f} mrad')
plt.axvline(sig.POINT_TIGHT_SEP_IN * 1000, color='blue', ls=':', lw=1.5,
            label=f'tight cut (sep_in<10cm) = {sig.POINT_TIGHT_SEP_IN*1000:.0f} mrad')
plt.xscale('log')
plt.xlabel('pointing angle (mrad)')
plt.ylabel('density (rate-weighted)')
plt.title('Pointing for surviving events (pre-pointing selection)')
plt.legend()
plt.tight_layout()
plt.savefig('pointing_compare.png', dpi=150)
print('wrote pointing_compare.png')
