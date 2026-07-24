# BC1 (kinetic-mixing dark photon) — status

_Last reviewed 2026-07-24. Active worktree:
`signal-models/bc1-darkphoton/llpatcolliders_BC1`, branch `bc1-darkphoton`
(3 commits on top of `bc4-scalar`, pushed to `myfork`; no PR yet — intended PR
base `bc4-scalar`). Package `dark_photon/`._

## Physics definition

BC1 is the minimal kinetic-mixing dark photon: one parameter, `epsilon^2`, at
fixed `m_{A'}` controls production, lifetime, and all visible branching
fractions. The signal criterion is `N_signal >= 3` at `3000 fb^-1` under the
paper's explicit zero-background working assumption. Because production and
lifetime are locked to the same coupling, the reach is a **closed island** in
`(m_{A'}, epsilon^2)` (a lower edge = too little production, an upper edge =
decays before reaching the detector).

## What is done

- **Meson production** (low mass, `~0.02–0.14 GeV`): validated; spectrum
  `data/spectra/meson_softqcd.npz` (500k Pythia 8.315 SoftQCD events) committed.
  Deepest reach `epsilon^2 ~ 4.4e-14` at `m_{A'} ~ 0.07 GeV`.
- **Drell-Yan production** (higher mass): spectra `data/spectra/dy/dy_*.npz`
  committed for 272 masses spanning `1.65–50 GeV`, generated with the ATLAS
  HAHM `GFfix` UFO via MG5 `p p > zp` + Pythia ISR shower (provenance in `README.md`).
- **Model layer** (`model.py`): `c*tau`, production and visible BRs from the
  DeLiVeR VMD table `data/br_tables/dp_brs_deliver.csv`.
- **Existing-constraint overlay**: `data/constraints/existing_exclusion_senscalc.npz`
  (SensCalc v1.3.3 DP contour, now vendored under `shared/vendored/`).
- **Higgs portal `h → A'A'`** (`higgs_portal/`): **complete and published.**
  Best excluded BR ~ `7.5e-5` at `c*tau ~ 0.95 m` (`m_{A'} = 2 GeV`, HL-LHC
  3 ab⁻¹). Curve committed as `higgs_portal/results/grendel_hAA_2GeV.csv`;
  method, run steps and comparison in `higgs_portal/README.md`.

## Remaining work

- **P0 — publish the BC1 island.** `run_sensitivity.py` computes the closed
  `(m_{A'}, epsilon^2)` island, but the result currently lives only in the
  gitignored `tmp/`. To reach BC4 fidelity it must be promoted to a tracked
  canonical CSV under `data/published/` with an adjacent `MANIFEST.json`
  (source, config, checksum, grid, convergence), mirroring
  `scalar/data/published/bc4_island.csv`.
- **P0 — join the meson and DY channels** into a single continuous exclusion
  across the `~0.14 GeV` handover, with the boundary/overlap policy documented.
- **P1 — DY low-mass integration.** The 2→1 MG5 integration was unstable at low
  mass ("Bjorken x>1"); confirm DY is only relied on above the meson endpoint,
  or stabilise it.
- **P1 — uncertainty band** analogous to BC4's variation bundle (not yet built
  for BC1).
- **P2 — `h → A'A'` cross-checks:** run `higgs_portal/verify_reach.py` (the
  independent geometric cross-check of the `7.5e-5` point, still un-run); extend
  the portal scan to other `m_{A'}` (only 2 GeV done). The `A' → e+e-` topology
  is conservative on the min-separation cut — a real 2 GeV `A'` decays mostly to
  hadrons/`mumu` (wider angle, `BR_visible ~ 1`).
- **P2 — detector response and backgrounds** are not yet validated at
  publication fidelity; zero background is a working assumption (shared with all
  channels).

For rerun commands see `README.md` (core) and `higgs_portal/README.md` (portal);
production provenance is in `README.md`.
