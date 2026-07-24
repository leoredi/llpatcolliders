"""Independent cross-check of the h->A'A' GRENDEL reach (no heavy MC).

Recomputes the geometric acceptance and an analytic best-case excluded BR from
first principles, and compares to the saved curve, to confirm the headline
number is self-consistent.
"""
import sys
from pathlib import Path
import numpy as np

_HERE = Path(__file__).resolve().parent
_HIGGS = _HERE.parents[1] / "higgs"
sys.path.insert(0, str(_HIGGS))

import decayProbPerEvent_2body as sig
from grendel_geometry import mesh_fiducial

CSV = _HERE / "tmp" / "llp_hAA_2GeV.csv"
CURVE = _HERE / "tmp" / "grendel_hAA_2GeV.csv"
C = sig.SPEED_OF_LIGHT
N_H = 3000 * 52e3   # HL-LHC ggF Higgs count used by the machinery

geo = sig.cache_geometry(str(CSV), mesh_fiducial, [0, 0, 0])
hits = geo["hits"]
n_tot = len(hits); n_hit = int(hits.sum())
entry = geo["entry_d"][hits]; exit_ = geo["exit_d"][hits]
bg = geo["gamma"][hits] * geo["beta"][hits]
path = exit_ - entry

print(f"A' total={n_tot}  hitting tunnel={n_hit}  ({100*n_hit/n_tot:.2f}%)")
print(f"  entry distance [m]:  median={np.median(entry):.1f}  mean={entry.mean():.1f}"
      f"  range=[{entry.min():.1f},{entry.max():.1f}]")
print(f"  fiducial path  [m]:  median={np.median(path):.2f}  mean={path.mean():.2f}")
print(f"  betagamma of hitters: median={np.median(bg):.1f}  mean={bg.mean():.1f}")

# Per-event geometric prob that >=1 A' points into the tunnel (2 A'/event).
n_events = 20000
frac_evt_ge1 = n_hit / n_events  # ~ P(>=1 in tunnel) for small hit prob
print(f"\nP(>=1 A' in tunnel) per Higgs ~ {frac_evt_ge1:.4f}")

# Best-case decay-in-fiducial probability, scanning ctau analytically per hitter:
#   P_decay_in = exp(-entry/lam) - exp(-exit/lam),  lam = bg*ctau
ctaus = np.logspace(-2, 2, 200)
best_mean_pdecay = 0.0; best_ct = None
for ct in ctaus:
    lam = bg * ct
    pdec = np.exp(-entry / lam) - np.exp(-exit_ / lam)
    m = pdec.mean() * n_hit / n_events      # per-Higgs, summed over hitters
    if m > best_mean_pdecay:
        best_mean_pdecay, best_ct = m, ct
# no selection eff here -> upper bound on sensitivity (lower bound on BR)
br_geom_only = 3.0 / (best_mean_pdecay * N_H)
print(f"Best per-Higgs decay-in-tunnel prob (geom, no cuts) = {best_mean_pdecay:.3e}"
      f"  at ctau~{best_ct:.2g} m")
print(f"=> geometry-only best excluded BR (no selection) = {br_geom_only:.2e}")

# Compare to the saved full-selection and acceptance-only curve.
d = np.loadtxt(CURVE, delimiter=",", skiprows=1)
ct, brf, bra = d[:, 0], d[:, 1], d[:, 2]
bi = np.nanargmin(brf); ai = np.nanargmin(bra)
print(f"\nSaved curve:")
print(f"  full selection : best BR={brf[bi]:.2e} at ctau={ct[bi]:.2g} m")
print(f"  acceptance only: best BR={bra[ai]:.2e} at ctau={ct[ai]:.2g} m")
print(f"  implied selection eff at min (acc/full) = {bra[ai]/brf[bi]:.2f}")
print(f"  ratio full/geom-only = {brf[bi]/br_geom_only:.2f}  "
      f"(should be >1: selection + reco cost)")
