# Higgs-portal dark photons — `h -> A'A'` long-lived reach (GRENDEL)

Purely-additive sub-package that evaluates GRENDEL's sensitivity to a
**long-lived dark photon `A'` produced in exotic Higgs decays**, `h -> A'A'`,
in the PBC "Higgs decay to long-lived dark photons" plane: **BR(`h -> A'A'`)
vs. proper decay length `c*tau_{A'}`** at fixed `m_{A'} = 2 GeV`.

Everything geometric, reconstruction-related and the exclusion normalisation is
**imported, not copied**, from the shared single source:

* `higgs/decayProbPerEvent_2body.py` — `cache_geometry` (ray-cast into the PX56
  fiducial volume), `sample_separations` (one uniform MC pass over decay points),
  `analyze_decay_vs_lifetime`, and `mc_exclusion_vs_lifetime` (reweight the single
  pass to every `c*tau`). The LLP is modelled as a 2-body `A' -> e+e-` decay with
  the same `P_CUT / SEP_MIN / SEP_MAX` selection used for the dark-scalar curves.
* `higgs/grendel_geometry.py` — GRENDEL PX56 mesh + ray-casting.

## Why this plane, and not the BC1 island

BC1 (the minimal kinetic-mixing dark photon, `dark_photon/` core) is
**coupling-controlled**: a single `eps^2`, at fixed `m_{A'}`, fixes the
production rate *and* the lifetime simultaneously, so the reach collapses to a
closed `(m_{A'}, eps^2)` island — a poor probe of a decay-length detector.

The **Higgs portal decouples them**: the production BR(`h -> A'A'`) is a free
parameter and the `A'` lifetime is set independently by its own (small) coupling.
That gives the genuine 2-D BR-vs-`c*tau` plane, which is exactly what maps out a
displaced detector's fiducial acceptance. The exclusion is set
model-agnostically, identically to BC5 (`higgs/`, `h -> SS`):

```
excluded BR = 3 / (N_h * P(>=1 A' decays visibly & passes selection))
N_h = sigma_h(ggF) * L = 52 pb * 3000 fb^-1 = 1.56e8
```

i.e. 3 signal events (background-free 95% CL) at HL-LHC. This is exactly the PBC
y-axis.

## How to run

```bash
# 1. generate the production sample (Pythia8 ggF Higgs, h -> A'A' forced 100%)
#    run in an env with pythia8 (e.g. llpatcolliders_FONLL)
python gen_higgs_darkphoton.py            # -> tmp/llp_hAA_2GeV.csv (40k A')

# 2. cast, sample, and reweight to every c*tau -> the reach curve
python run_higgs_reach.py                 # -> tmp/grendel_hAA_2GeV.csv (+plots)

# 3. overlay on the PBC reference curves (CODEX-b / MATHUSLA / ANUBIS + h->inv)
python plot_overlay_reach.py              # -> tmp/grendel_hAA_2GeV_overlay.{png,pdf}
```

Raw MC and working plots land in `tmp/` (gitignored); the published snapshot is
committed under `results/`.

## Result (`m_{A'} = 2 GeV`, HL-LHC, 3 ab^-1)

Best excluded **BR(`h -> A'A'`) ~ 7.5e-5 at `c*tau ~ 0.95 m`** (full selection),
with BR < 1e-3 over `c*tau` in ~[0.14, 41] m. GRENDEL turns over at the
**shortest `c*tau` of the HL-LHC set** (its ~22 m stand-off + metre-scale fiducial
depth): it beats CODEX-b (~5x) and MATHUSLA (~2x) at their optima; only ANUBIS
(a much larger PX14-shaft volume) reaches deeper, at longer `c*tau`. All curves
are HL-LHC 3 ab^-1 projections at their respective IPs, so the comparison is
apples-to-apples on luminosity — it is a geometry/acceptance statement.

## Layout

| file | role |
|------|------|
| `gen_higgs_darkphoton.py` | Pythia8 ggF Higgs generator; `h` forced 100% to a `6000113` (`A'`) pair, `A'` made stable so its production 4-vector is recorded (the machinery supplies the lifetime). Writes `event,id,pt,eta,phi,momentum,mass`. |
| `run_higgs_reach.py` | driver: `cache_geometry -> sample_separations -> mc_exclusion_vs_lifetime` over a `c*tau` grid covering [1e-3, 1e3] m. Writes the curve CSV + plots. |
| `plot_reach.py` | fast re-plot of the GRENDEL-only curve from the saved CSV. |
| `plot_overlay_reach.py` | overlay vs the digitized PBC curves in `higgs/external/*.csv` (CODEX-b, MATHUSLA, ANUBIS) + the `h->inv` HL-LHC cap + an existing-constraint band from `ATLAS_current`/`CMS_current`. |
| `verify_reach.py` | independent geometric cross-check of the best-point number. |
| `results/` | published snapshot: the reach curve CSV and the GRENDEL/overlay plots. |
