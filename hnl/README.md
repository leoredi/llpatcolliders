# `hnl/` — Heavy Neutral Lepton production from rare meson decays

HNL 4-vector production pipeline for GRENDEL sensitivity studies at the
HL-LHC (pp, sqrt(s) = 14 TeV).

## Scope of this PR

**In scope** — production from rare meson decays:

- Direct: `B -> N + X`, `D -> N + X`, `Bc -> N + X` (2-body and 3-body channels)
- Induced tau: `Ds -> tau nu -> N + X`, `B+ -> tau nu -> N + X`

**Out of scope** — deferred to follow-up PRs:

- W/Z -> ell N direct production
- Drell-Yan, `pp -> tau tau`
- Analysis, sensitivity scan, exclusion plots
- Geometry, decay-probability, daughter ray-casting

## Directory layout

    hnl/
    |-- README.md                this file
    |-- run_all.py               parallel orchestrator (entry point)
    |-- config_mass_grid.py      MASS_GRID definition (116 points)
    |-- .gitignore               ignores output/
    |
    |-- production/
    |   |-- constants.py         meson masses, FRAG_B / FRAG_C, sigma_Bc, BRs
    |   |-- combine_channels.py  vstack {Bmeson,Dmeson,Bc,tau} -> combined/
    |   |
    |   |-- fonll/
    |   |   |-- fonll_parser.py   parse 2D dsigma/dpT/dy table, integrate
    |   |   |-- meson_sampler.py  inverse-CDF sample meson 4-vectors + species assign
    |   |
    |   |-- decay_engine/
    |       |-- kinematics.py            2-body, 3-body flat phase space, boost
    |       |-- generate_meson_csvs.py   driver: B, D, Bc -> N
    |       |-- generate_induced_tau.py  driver: Ds, B+ -> tau -> N
    |
    |-- tests/                   parser, sampler, kinematics, end-to-end smoke
    |
    |-- vendored/
        |-- PROVENANCE.md                FONLL form params, HNLCalc upstream
        |-- HNLCalc/                     pure-Python HNL BR computation
        |-- fonll_pp14tev_nnpdf40_nlo_as_01180_..._charm.dat
        |-- fonll_pp14tev_nnpdf40_nlo_as_01180_..._bottom.dat
        |-- fonll_pp14tev_cteq66_..._charm.dat
        |-- fonll_pp14tev_cteq66_..._bottom.dat

## Pipeline (per (flavor, channel, mass) point)

    1. Sample N meson 4-vectors from FONLL dsigma/dpT/dy table
                                                  (production/fonll/meson_sampler.py)
    2. Assign species (D0/D+/Ds  or  B+/B0/Bs) via FRAG_C / FRAG_B
                                                  (production/constants.py)
    3. For each species, compute BR(meson -> N + X) at U^2 = 1 via HNLCalc
                                                  (production/decay_engine/generate_meson_csvs.py
                                                   compute_production_br_components)
    4. Sample 2-body vs 3-body channel weighted by BR
                                                  (_sample_hnl_from_mesons)
    5. Decay meson -> ell N (2-body) or meson -> ell N H' (3-body) in lab frame
                                                  (production/decay_engine/kinematics.py)
    6. Emit HNL 4-vector with weight (formula below)
                                                  (CSV writer)

For induced tau (`generate_induced_tau.py`): steps 1-2 are the same but step 3
is a 2-body decay `Ds/B+ -> tau nu_tau`, then steps 3-6 are repeated with the
tau as the parent and `tau -> N + X` channels via HNLCalc's
`get_2body_br_tau` / `get_3body_dbr_tau`.

## Weight formula

Per HNL 4-vector i, in pb at `U^2 = 1`:

**Direct meson channel (B, D, Bc -> N):**

    w_i = 2 * sigma_FONLL(quark) * f_species * BR(meson -> N+X | U^2=1) / N_species_sampled

- `sigma_FONLL(quark)` from `production/fonll/fonll_parser.py::get_sigma_total`
  (numerical 2D integration of the vendored grid)
- Factor 2: FONLL gives `(quark + antiquark)/2`; we want both -> x 2
- `f_species` from `FRAG_C` / `FRAG_B` in `constants.py`
- For Bc: `sigma_FONLL(bottom) * f_Bc` is replaced by `SIGMA_BC_PB` (separate
  measured input, ~0.9 ub)

**Induced tau channel (Ds, B+ -> tau -> N):**

    w_i = 2 * sigma_FONLL(quark) * f_parent * BR(parent -> tau nu)
            * BR(tau -> N+X | U^2=1) / N_tau_sampled

with `BR(Ds -> tau nu) = 5.35e-2`, `BR(B+ -> tau nu) = 1.09e-4` (PDG 2024).

**Final yield at chosen `U^2` and luminosity:**

    N_HNL_per_event = w_i * U_alpha^2 * (epsilon_decay_at_ctau(U^2))
    N_signal = L_int * sum_i  w_i * U_alpha^2 * P_decay_i(ctau(U^2))

The `U^2`-dependent decay probability is the consumer's responsibility
(deferred to the analysis PR). Production weights stored here factor out
`U_alpha^2` cleanly.

## Channel inventory

**2-body meson production** (`compute_production_br_components`,
`_eval_2body_br`): all parents in `MESON_MASSES` with
`m_N < m_parent - m_lepton` get a direct `meson+ -> ell+ N` channel.

**3-body meson production** (`THREEBODY_CHANNELS` in
`generate_meson_csvs.py`): semi-leptonic channels with a pseudoscalar or
vector hadronic daughter. Inventory:

| Parent | 3-body daughters considered (X in `parent -> X ell N`) |
|---|---|
| D0   (421) | K-, K*-, pi-, rho-                       |
| D+   (411) | Kbar0, Kbar*0, pi0, rho0, eta, eta'      |
| Ds+  (431) | eta, eta', K0, K*0, phi                  |
| B+   (521) | Dbar0, Dbar*0, pi0, rho0, eta, omega, eta' |
| B0   (511) | D-, D*-, pi-, rho-                       |
| Bs   (531) | Ds-, Ds*-, K-, K*-                       |
| Bc+  (541) | B0, Bs, B*0, Bs*, D0, etac, D*0, J/psi   |

**Tau decay channels into N** (`generate_induced_tau.py`):

- 2-body hadronic: `tau- -> {pi-, K-, rho-, K*-} + N`
- 3-body leptonic: `tau- -> ell- nu_tau N` and `tau- -> ell- nubar_ell N`
  for `ell in {e, mu}`

All BRs are evaluated at `U^2 = 1` via HNLCalc and stored that way; the
explicit `U^2` scan is downstream.

## Output format

All drivers write headerless CSVs with five columns:

    weight, E, px, py, pz

- `weight` is in pb at `U^2 = 1`.
- `(E, px, py, pz)` is the HNL 4-momentum in the lab frame (CMS origin),
  GeV units, beam along z.

Files end up at:

    hnl/output/llp_4vectors/{Ue,Umu,Utau}/{Bmeson,Dmeson,Bc,tau}/mN_{mass}.csv
    hnl/output/llp_4vectors/{Ue,Umu,Utau}/combined/mN_{mass}.csv

`mN_{mass}` uses the encoding from `config_mass_grid.format_mass_for_filename`
(e.g. `1.025 GeV -> mN_1p03.csv`). Empty file = mass point below all
production thresholds for that channel.

## Setup

Requires Python 3.10+ with: `numpy`, `pandas`, `scipy`, `sympy`, `mpmath`,
`particle`, `numba`, `pytest` (tests only). Recommended via conda
(matching the existing repo convention):

    conda env list  # check if you already have it
    # or set up a new env with the project packages

For development the maintainer uses `conda env llpatcolliders_FONLL`;
substitute your own env name in the run commands below.

The vendored `HNLCalc` package is bundled in `hnl/vendored/HNLCalc/` and
imported via a `sys.path.insert` from each driver. No pip install needed
for the vendored deps.

## Usage

### Full grid in parallel (recommended)

    cd hnl
    python run_all.py                       # all flavors, all 116 masses, all cores
    python run_all.py --workers 8           # cap concurrency
    python run_all.py --masses 0.5 1.0 2.0  # subset of masses for a quick check
    python run_all.py --n-pool 50000        # smaller per-worker pool (faster, noisier)
    python run_all.py --skip-combine        # skip final combine step

`run_all.py` submits one process per (flavor, channel) pair (12 jobs total:
3 flavors x 3 meson channels + 3 flavors x 1 tau channel), then runs the
combine step at the end. The bottleneck is HNLCalc's 3-body BR integration,
which dominates over meson sampling; pool generation is duplicated per worker
but cheap.

### Individual drivers (serial, one channel at a time)

    python -m production.decay_engine.generate_meson_csvs \
        --flavor Ue --channel Dmeson --masses 1.0 --n-pool 100000
    python -m production.decay_engine.generate_induced_tau \
        --flavor Ue --masses 1.0 --n-pool 100000
    python -m production.combine_channels --flavor Ue --masses 1.0

Each driver takes `--flavor`, `--masses`, `--n-pool`, `--seed`.

## Tests

    pytest hnl/tests/

The test suite covers parser shape and total xsec, sampled-pool moments,
2-body and 3-body energy/momentum conservation, and a single end-to-end
smoke for `D -> N` at `m_N = 1 GeV`.

## Methodology references

- FONLL: Cacciari, Greco, Nason (NLO+NLL heavy-quark production)
- Central production tables: local FONLL+LHAPDF with
  `NNPDF40_nlo_as_01180`
- MATHUSLA reference files: davidrcurtin/MATHUSLA_LLPfiles_RHN_U{e,mu,tau}
  — legacy FONLL + CTEQ6.6 comparison baseline
- HNL branching ratios: Bondarenko, Boyarsky, Mikulenko, Naumov et al.
  1805.08567 "Phenomenology of GeV-scale HNLs" (implemented in HNLCalc)
- Tau-decay HNL channels: standard charged-current + neutral-current,
  evaluated through HNLCalc's `get_2body_br_tau` and `get_3body_dbr_tau`

See `hnl/vendored/PROVENANCE.md` for exact FONLL form parameters and
HNLCalc upstream commit.

## Known limitations

These are deliberate trade-offs in this PR and should be tracked explicitly
if revisited:

- **Default FONLL backend: NNPDF4.0 NLO** — the production parser now uses
  local FONLL+LHAPDF `NNPDF40_nlo_as_01180` tables by default. The old
  CTEQ6.6 tables remain available as the `cteq66_legacy` backend for
  reference comparisons by setting `HNL_FONLL_SET=cteq66_legacy` before
  starting Python.
- **FONLL grid pT_max = 50 GeV** — high-pT contribution above this is
  dropped. Negligible for forward / off-axis acceptance but matters for
  high-pT analyses; extend by regenerating tables on a wider grid.
- **Central scale, no PDF or scale uncertainty** — single central
  prediction in the vendored tables. Systematics deferred.
- **Species shapes** — charm uses the public-FONLL `D0` convention with
  calibrated `D*` feeddown as the central pT-y shape for D0, D+, and Ds.
  D+ and Ds species fractions are applied in weights, but dedicated D+/Ds
  pT-y shape variations are still a systematic, not a separate central grid.
  The FONLL output convention used here is `ifrframe = 1` (y=0-frame
  fragmentation, the public-FONLL default for `dsigma/dpT/dy`). On that
  branch `xmh` enters only the `szmin` kinematic guard in `fragmfonll.f`,
  which stays inactive across the (`pT <= 50 GeV`, `|y| <= 3`) grid here,
  so the table is empirically insensitive to `xmh` in this phase space.
  The sampler reconstructs each meson four-vector with the physical
  on-shell PDG mass, which is the kinematically consistent thing to do
  for this table convention. Extending the grid toward the kinematic
  edge of `sh/4` would re-activate the guard and require re-verifying
  this insensitivity.
- **Static bottom fragmentation fractions** — the omitted bottom-baryon
  remainder is derived from LHCb's pT-averaged Lambda_b ratio over
  `4 < pT < 25 GeV`, `2 < eta < 5` and applied over the full table. Its
  measured pT dependence and extrapolation outside that acceptance are not
  modeled. A detector-level refinement must update sampling and event weights
  consistently.
- **No NLO matching for the production decay** — `decay_3body_flat`
  is flat phase space (no matrix-element weighting). Acceptable for
  total BR-weighted yields; biased for differential distributions in
  the decay daughters.
- **Tau polarization neglected** — Ds -> tau nu and B+ -> tau nu produce
  longitudinally polarized taus, but `decay_2body` / `decay_3body_flat`
  decay them isotropically in the tau rest frame. This washes out the
  parent-N angular correlation and biases the HNL energy spectrum at
  fixed parent boost by O(20-30 %). Acceptable for total yields; biased
  for differential angular distributions.
- **B0 -> tau nu omitted** — helicity-suppressed in SM; only B+ -> tau nu
  retained in the induced-tau chain.
- **W -> tau nu omitted** — moves to the W/Z PR. This is the dominant
  prompt-tau source at LHC and the induced-tau chain alone underestimates
  the tau parent yield.
- **No kaon (K+ -> lN) production** — MASS_GRID starts at 0.20 GeV, but
  the two-body kaon decay K+ -> lN is the dominant HNL source below
  ~0.5 GeV and is not included. Sub-0.5 GeV yields are therefore
  underestimated.
- **No Lambda_b / Lambda_c / Xi_c baryon channels** — pp fragmentation
  fractions now track the omitted baryon component explicitly, but only meson
  parents are simulated in this production layer.
