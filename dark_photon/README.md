# BC1 — kinetic-mixing dark photon (GRENDEL)

Purely-additive package implementing the PBC benchmark **BC1**: a dark photon
`A'` that couples to the SM through kinetic mixing `epsilon` with the photon,
parameters `(m_{A'}, epsilon^2)`. Like BC4/BC10 the minimal model is
**coupling-controlled** — the single `epsilon^2`, at fixed `m_{A'}`, fixes the
production rate *and* the lifetime *and* the visible BRs simultaneously — so the
reach is a closed `(m_{A'}, epsilon^2)` island.

Everything geometric and reconstruction-related is **imported, not copied**, from
the shared single sources (`higgs/grendel_geometry.py`, `higgs/reco_common.py`,
`hnl/analysis/…`, `hnl/production/…`) exactly as the `scalar/` (BC4) package does.

> **Import gotcha:** BC1 modules must be imported package-qualified
> (`from dark_photon import …`). Never put `dark_photon/` on `sys.path` — its
> `production.py` would shadow the `hnl` `production/` package the shared engine
> imports.

## Two production channels

| channel | mass range | how | status |
|---------|-----------|-----|--------|
| **Meson** | low (`~0.02–0.14 GeV`) | `pi0/eta/omega` from Pythia8 SoftQCD → `P → A' gamma`, `omega → A' pi0`. Committed spectrum `data/spectra/meson_softqcd.npz` (500k events, Pythia 8.315). | done + validated |
| **Drell-Yan** | higher (`~1.65–50 GeV`) | MG5 `p p > zp` hard process + Pythia ISR shower (Berlin/Kling style — the 2→1 parton process alone gives `pT=0` and misses `|eta|<0.5`). ATLAS **HAHM_darkphoton_LJmod_UFO_GFfix** UFO (dark photon = `zp`, **pdg 3000001**; mixing = HIDDEN block `epsilon`). Spectra `data/spectra/dy/dy_*.npz` (272 masses). | spectra produced |

## Higgs portal — `h → A'A'` (sub-package `higgs_portal/`)

The minimal-mixing island is a poor probe of a decay-length detector (production
and lifetime are locked together). The **Higgs portal decouples them**: BR(`h →
A'A'`) is free and the `A'` lifetime is independent, giving the genuine
BR-vs-`c*tau` plane. See `higgs_portal/README.md`. This is the strongest GRENDEL
story for BC1-type dark photons; result: best excluded BR ~ `7.5e-5` at
`c*tau ~ 0.95 m` (`m_{A'} = 2 GeV`), published under `higgs_portal/results/`.

## Layout

| file | role |
|------|------|
| `model.py` | model layer: `c*tau(m_{A'}, epsilon^2)`, production and visible BRs, from the DeLiVeR VMD table `data/br_tables/dp_brs_deliver.csv`. The BC1 analogue of HNLCalc / the BC4 `model.py`. |
| `production.py` | meson + Drell-Yan sampling → weighted `A'` four-vector CSVs. |
| `acceptance.py` | decay-in-flight × two-body acceptance on the shared GRENDEL geometry/reco (mirrors `scalar/acceptance.py`). |
| `run_sensitivity.py` | mass × `epsilon^2` scan → the closed exclusion island. |
| `plot_exclusion.py` | island plot, incl. the existing-exclusion overlay `data/constraints/existing_exclusion_senscalc.npz`. |
| `generator/` | `dp_meson_softqcd.cc`, `make_meson_spectrum.py`, `dy_madgraph.py` + cards. |
| `higgs_portal/` | the `h → A'A'` long-lived reach sub-study (see its README). |
| `tests/` | package tests. |

See `STATUS.md` for the publication state and remaining work.
