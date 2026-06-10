"""
Per-track lever length (inner-hit -> outer-hit distance = L/|d.n_hat|) for events
surviving the full selection, comparing the cosmic decay-in-flight background to
several signal masses. Tests whether the long-lever (oblique/skimming) tracks
that now dominate are a cosmic-only feature or also present in signal.

Cosmic: 500k volume events (p<1 GeV window) for a quick test.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import cosmic_decay_check as cdc
import decayProbPerEvent_2body as sig
import reco_common
from grendel_geometry import mesh_fiducial

L = sig.DETECTOR_THICKNESS

# ---------------- cosmic ----------------
cdc.GEN_PMIN, cdc.GEN_PMAX = 0.0, 1.0
flo = cdc.empirical_window_fraction(0.0, 1.0)
r = cdc.run_volume(500_000, np.deg2rad(80), 1.0, 'empirical', 70.0,
                   cdc.HIT_RESOLUTION, 42, sigma_t=1.5e-9,
                   target_muon_rate_hz=600.0, chunk=200_000, spectrum_weight=flo)
m_cos = cdc.cutflow_stages(r)[-1][1]          # final (timing) survivors
stub_mu = np.linalg.norm(r['P_eo'] - r['P_ei'], axis=1)
stub_e = np.linalg.norm(r['P_xo'] - r['P_xi'], axis=1)
cos_len = np.concatenate([stub_mu[m_cos], stub_e[m_cos]])
cos_w = np.concatenate([r['w'][m_cos], r['w'][m_cos]])
print(f"cosmic full-selection survivors: {int(m_cos.sum())} events")


# ---------------- signals ----------------
def sig_lengths(csv):
    geo = sig.cache_geometry(csv, mesh_fiducial, [0, 0, 0])
    mc = sig.sample_separations(geo, 1e-7, n_samples_per_particle=200)
    m = sig.selection_mask(mc)
    in1, out1 = reco_common.wall_inner_outer(mc['exit_pt_1'], mc['dir1'], L)
    in2, out2 = reco_common.wall_inner_outer(mc['exit_pt_2'], mc['dir2'], L)
    s1 = np.linalg.norm(out1 - in1, axis=1)
    s2 = np.linalg.norm(out2 - in2, axis=1)
    lens = np.concatenate([s1[m], s2[m]])
    ws = np.concatenate([mc['weights'][m], mc['weights'][m]])
    fin = np.isfinite(lens)
    return lens[fin], ws[fin]


SIG = [('0.5 GeV', 'LLP0p5GeVSmall.csv', 'tab:green'),
       ('1 GeV',   'LLP1GeVSmall.csv',   'tab:olive'),
       ('15 GeV',  'LLPSmall.csv',       'tab:blue'),
       ('40 GeV',  'LLP40GeVSmall.csv',  'tab:purple')]

bins = np.logspace(np.log10(L), np.log10(20.0), 55)   # L .. 20 m
plt.figure(figsize=(9, 6))
fin = np.isfinite(cos_len)
plt.hist(cos_len[fin], bins=bins, weights=cos_w[fin], density=True,
         histtype='step', lw=2.5, color='k', label='cosmic')
for lbl, csv, c in SIG:
    l, w = sig_lengths(csv)
    plt.hist(l, bins=bins, weights=w, density=True, histtype='step',
             lw=2, color=c, label=lbl)
    print(f"{lbl}: median track length {np.median(l):.2f} m, "
          f"90% {np.percentile(l, 90):.2f} m, max {l.max():.2f} m")
print(f"cosmic: median {np.median(cos_len[fin]):.2f} m, "
      f"90% {np.percentile(cos_len[fin], 90):.2f} m, max {cos_len[fin].max():.2f} m")

plt.axvline(L, color='gray', ls=':', label=f'layer spacing L = {L*100:.0f} cm')
plt.xscale('log')
plt.xlabel('per-track length: inner -> outer hit (m)')
plt.ylabel('density (rate-weighted)')
plt.title('Track length (lever arm) for surviving events')
plt.legend()
plt.tight_layout()
plt.savefig('track_lengths.png', dpi=150)
print('wrote track_lengths.png')
