# `hnl/` — Heavy Neutral Lepton production from meson decays and electroweak (W/Z) processes

HNL 4-vector production pipeline for GRENDEL sensitivity studies at the
HL-LHC (pp, sqrt(s) = 14 TeV).

## Scope of this PR

**In scope** — production from rare meson decays plus electroweak production:

- Direct meson: `B -> N + X`, `D -> N + X`, `Bc -> N + X` (2-body and 3-body)
- Kaon: `K+ -> ell+ N` (2-body, dominant source below ~0.5 GeV)
- Induced tau: `Ds -> tau nu -> N + X`, `B+ -> tau nu -> N + X`
- Electroweak: `W/Z -> ell N` via MadGraph (`production/madgraph/`, default-on; opt out with `--no-wz`)
- Prompt tau: `pp -> W -> tau nu` + Drell-Yan `pp -> tau tau` -> `tau -> N + X`
  via MadGraph (`production/madgraph/run_tau_production.py`, default-on; opt out
  with `--no-prompt-tau`)

**Out of scope** — deferred to follow-up PRs:

- Direct Drell-Yan production `pp -> N ell` (without the tau intermediary)
- Analysis, sensitivity scan, exclusion plots
- Geometry, decay-probability, daughter ray-casting

## Directory layout

    hnl/
    |-- README.md                this file
    |-- run_all.py               parallel orchestrator (entry point)
    |-- config_mass_grid.py      MASS_GRID + filename encoder/decoder
    |-- .gitignore               ignores output/
    |
    |-- production/
    |   |-- constants.py         meson masses, FRAG_B / FRAG_C, sigma_Bc/K, BRs
    |   |-- combine_channels.py  vstack channels -> combined/
    |   |
    |   |-- fonll/
    |   |   |-- fonll_parser.py   parse 2D dsigma/dpT/dy table, integrate
    |   |   |-- meson_sampler.py  inverse-CDF sample meson 4-vectors + species assign
    |   |
    |   |-- decay_engine/
    |   |   |-- kinematics.py            2-body (incl. polarized), 3-body flat, boost
    |   |   |-- tau_decay.py             shared tau -> N + X sampler (HNLCalc-driven)
    |   |   |-- generate_meson_csvs.py   driver: B, D, Bc -> N
    |   |   |-- generate_induced_tau.py  driver: Ds, B+ -> tau -> N   (induced_tau/)
    |   |   |-- generate_kaon_csvs.py    driver: K+ -> N + X
    |   |
    |   |-- madgraph/                    MG5-driven pipelines (default-on)
    |       |-- _mg5_common.py           shared LHAPDF/dyld plumbing
    |       |-- run_wz_production.py     W/Z -> l N                    (WZ/)
    |       |-- run_tau_production.py    prompt tau (W,Z -> tau)       (tau/)
    |       |-- lhe_to_csv.py
    |       |-- cards/                   proc/run/param templates
    |
    |-- tests/                   parser, sampler, kinematics, smoke, kaon, mass-grid
    |
    |-- vendored/
        |-- PROVENANCE.md                FONLL form params, HNLCalc, HeavyN UFO, tau pool
        |-- HNLCalc/                     pure-Python HNL BR computation
        |-- SM_HeavyN_CKM_AllMasses_LO/  MadGraph UFO model for W/Z driver
        |-- tau_pool.csv                 prompt-tau LHE pool (Stage 1 product)
        |-- fonll_pp14tev_nnpdf40_nlo_as_01180_..._{charm,bottom}.dat
        |-- fonll_pp14tev_cteq66_..._{charm,bottom}.dat   # legacy comparison

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

**Induced tau channel (Ds, B+ -> tau -> N, output: `induced_tau/`):**

    w_i = 2 * sigma_FONLL(quark) * f_parent * BR(parent -> tau nu)
            * BR(tau -> N+X | U^2=1) / N_tau_sampled

with `BR(Ds -> tau nu) = 5.35e-2`, `BR(B+ -> tau nu) = 1.09e-4` (PDG 2024).

**Prompt tau channel (W -> tau nu, gamma*/Z -> tau tau -> tau -> N, output: `tau/`):**

    w_i = (sigma_LO_MG5_taupool / N_pool) * K_FACTOR_EW
            * BR(tau -> N+X | U^2=1)

The first factor is the per-event LHE weight from a single shared MG5 tau
pool (Stage 1, vendored at `vendored/tau_pool.csv`); Stage 2 multiplies by
`BR(tau -> N+X)` per (flavor, m_N) point. The tau pool stores the mother
PDG so per-event polarisation can be applied: `W+/- mothers -> fully
polarised` (asymmetry = TAU_2BODY_ASYMMETRY = -1, same as Ds -> tau nu by
the V-A + CP argument in `tau_decay.py`); any other mother (Z, gamma*,
direct-DY initial-state quarks) -> `unpolarised` (asymmetry = 0).

**Kaon channel (K+ -> ell+ N, K+ -> pi0 ell+ N):**

    w_i = SIGMA_KAON_PB * BR(K -> N+X | U^2=1) / N_pool

`SIGMA_KAON_PB` already bundles both charges (K+ and K-) — no extra factor 2.
At `U^2 = 1` the HNLCalc "BR" exceeds 1 for K (the partial width is much
larger than the SM total K width), which is by construction: the per-event
weight is a rate-equivalent that becomes a physical rate once multiplied by
the realistic `U^2 << 1` downstream. See `production/constants.py` for the
approximation status of `SIGMA_KAON_PB`.

**Electroweak channel (W/Z -> ell N, default-on):**

    w_i = (sigma_LO_MG5 / N) * K_FACTOR_EW

The per-row weight comes from the MG5 unweighted-event weight (`XWGTUP`) read
out of the LHE, then scaled by `K_FACTOR_EW` so the summed weight is the
NLO-corrected cross-section.

The W/Z driver follows three explicit conventions worth flagging because they
affect how the per-row weight should be interpreted:

- **PDF set** — `pdlabel = lhapdf`, `lhaid = 331700` (NNPDF40_nlo_as_01180,
  same set the FONLL meson tables use). The LO matrix element paired with an
  NLO PDF is the standard practice in the HNL literature; the residual
  ME-side correction is absorbed into `K_FACTOR_EW`.
- **W width left at the SM value** in `param_card_template.dat`
  (`DECAY 24 = 2.085 GeV`), even though the driver writes `U_alpha = 1`
  for the active flavor. Pascoli-Ruiz convention: keep Gamma_W fixed at
  SM so the MG5 cross-section factorises cleanly as
  `sigma = sigma_SM_prod * BR(W -> ell N | U^2 = 1)`, and downstream
  scaling `sigma -> sigma * U_alpha^2` recovers the physical rate at the
  consumer's chosen mixing.
- **K-factor is a flat 1.3 multiplier** (`K_FACTOR_EW` in
  `production/constants.py`), applied per event after LHE parsing. The
  true NLO/LO QCD K-factor for W/Z -> ell N drifts mildly with m_N
  (~1.30 in the resonant region, ~1.15 at the high-m_N tail); the flat
  value is the same approximation the upstream llpatcolliders_FONLL
  pipeline uses.

**Final yield at chosen `U^2` and luminosity:**

    N_HNL_per_event = w_i * U_alpha^2 * (epsilon_decay_at_ctau(U^2))
    N_signal = L_int * sum_i  w_i * U_alpha^2 * P_decay_i(ctau(U^2))

The `U^2`-dependent decay probability is the consumer's responsibility
(deferred to the analysis PR). Production weights stored here factor out
`U_alpha^2` cleanly.

## Channel inventory

**2-body meson production** (`compute_production_br_components`,
`_eval_2body_br`): all charged-pseudoscalar parents (K+, D+, Ds+, B+, Bc+)
with `m_N < m_parent - m_lepton` get a direct `meson+ -> ell+ N` channel.

**3-body meson production** (`THREEBODY_CHANNELS` in
`generate_meson_csvs.py`): semi-leptonic channels with a pseudoscalar or
vector hadronic daughter. Inventory:

| Parent | 3-body daughters considered (X in `parent -> X ell N`) |
|---|---|
| K+   (321) | pi0                                       |
| D0   (421) | K-, K*-, pi-, rho-                       |
| D+   (411) | Kbar0, Kbar*0, pi0, rho0, eta, eta'      |
| Ds+  (431) | eta, eta', K0, K*0, phi                  |
| B+   (521) | Dbar0, Dbar*0, pi0, rho0, eta, omega, eta' |
| B0   (511) | D-, D*-, pi-, rho-                       |
| Bs   (531) | Ds-, Ds*-, K-, K*-                       |
| Bc+  (541) | B0, Bs, B*0, Bs*, D0, etac, D*0, J/psi   |

For Bc+ the (B0, Bs, B*0, Bs*) entries are listed for completeness — at every
`m_N` in `MASS_GRID` they are kinematically closed (`m_Bc - m_B - m_lep < m_N`)
and contribute nothing in practice.

**Tau decay channels into N** (shared sampler: `production/decay_engine/tau_decay.py`):

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

    hnl/output/llp_4vectors/{Ue,Umu,Utau}/{Bmeson,Dmeson,Bc,Kmeson,induced_tau,tau}/mN_{mass}.csv
    hnl/output/llp_4vectors/{Ue,Umu,Utau}/WZ/mN_{mass}.csv         # MadGraph (default-on)
    hnl/output/llp_4vectors/{Ue,Umu,Utau}/combined/mN_{mass}.csv

`mN_{mass}` uses the encoding from `config_mass_grid.format_mass_for_filename`
(e.g. `1.025 GeV -> mN_1p025.csv`; three decimals so the 15/25-MeV grid
spacing is encoded without drift). `config_mass_grid.parse_mass_from_filename`
is the exact inverse — consumers that pair a CSV with a `ctau(m_N)` table
should use it instead of re-parsing the label by hand.

**Empty (zero-byte) files are intentional sentinels, not job failures.** A
zero-byte CSV means the mass point is below every production threshold for
that channel (e.g. `*/tau/mN_*.csv` above `m_tau`, or `Utau/Dmeson` where
`m_N > m_D - m_tau` everywhere). A glob-based consumer should treat
`stat().st_size == 0` as "channel closed here" and skip it; this is exactly
what `combine_channels.py` does.

Near the high-mass edge of a channel a **non-empty** file may contain fewer
than `n_pool` rows: only meson species with `m_parent - m_lepton > m_N`
survive the per-species threshold, so as `m_N` rises the surviving subset of
the pool shrinks (e.g. `Ue,Umu/Dmeson/mN_1p730.csv` keeps only the D0/D+/Ds
tail above threshold and lands around 40k of 100k rows). This is physics, not
truncation.

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

`run_all.py` submits one process per (flavor, channel) pair. Default channel
set per flavor: 3 meson channels (Bmeson, Dmeson, Bc) + Kmeson + induced_tau
+ prompt-tau (Stage 2) + W/Z = 7 jobs/flavor, so 21 jobs in the parallel
pool by default. Before the pool, prompt-tau Stage 1 runs once as a
serialised MG5 step that produces `vendored/tau_pool.csv` (gated on file
presence and row-count threshold), and after the pool the combine step runs
serially. Heavy paths opt-out with `--no-wz` and `--no-prompt-tau`. The
bottleneck is HNLCalc's 3-body BR integration, which dominates over meson
sampling; pool generation is duplicated per meson/kaon/induced-tau worker
but cheap.

### Individual drivers (serial, one channel at a time)

    python -m production.decay_engine.generate_meson_csvs \
        --flavor Ue --channel Dmeson --masses 1.0 --n-pool 100000
    python -m production.decay_engine.generate_induced_tau \
        --flavor Ue --masses 1.0 --n-pool 100000
    python -m production.decay_engine.generate_kaon_csvs \
        --flavor Umu --masses 0.30 --n-pool 100000
    python -m production.combine_channels --flavor Ue --masses 1.0

Each driver takes `--flavor`, `--masses`, `--n-pool`, `--seed`.

### Electroweak W/Z -> ell N (default-on, needs MadGraph)

    # one quick point (Umu, 1.0 GeV, 1000 events)
    python production/madgraph/run_wz_production.py --test
    # full grid for one flavor
    python production/madgraph/run_wz_production.py --flavor Umu
    # default pipeline; opt out with --no-wz if you don't have MG5
    python run_all.py

### Prompt tau (default-on, needs MadGraph)

    # smoke: builds vendored/tau_pool.csv (1k events) + 1 mass × 1 flavor
    python production/madgraph/run_tau_production.py --test
    # reuse cached pool, decay across all (flavor, mass)
    python production/madgraph/run_tau_production.py --skip-mg5
    # opt out of the prompt-tau Stage 1 + Stage 2 entirely
    python run_all.py --no-prompt-tau

### Migration note (induced tau)

Old runs wrote the induced-τ output to `output/.../tau/`. That folder is now
reserved for *prompt* τ. Before re-running `combine_channels` on a pre-
existing checkout: either delete the old `tau/` CSVs or move them to
`induced_tau/`, otherwise you will double-count between the two paths once
the new prompt-τ job repopulates `tau/`. A fresh `run_all.py` regenerates
everything in the correct folders and is the safest path.

The MadGraph executable is resolved at runtime (see
`production/madgraph/_mg5_common.py::_resolve_mg5_exe`): `$HNL_MG5_EXE`,
then a drop-in `vendored/MG5_aMC_v3_6_6/`, then the sibling
`llpatcolliders_FONLL` install. The 148 MB MG5 tree is not committed; only the
460 KB HeavyN UFO model (`vendored/SM_HeavyN_CKM_AllMasses_LO/`) is vendored.
MG5 is invoked as a plain subprocess — no Docker or other container runtime is
needed. Output lands at `output/llp_4vectors/{flavor}/WZ/mN_*.csv` and is
picked up automatically by the combine step.

## Tests

    pytest hnl/tests/

The test suite covers parser shape and total xsec, sampled-pool moments,
2-body and 3-body energy/momentum conservation, a kaon-channel smoke,
mass-grid encoder round-trip, and an end-to-end smoke for `D -> N` at
`m_N = 1 GeV`.

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
- **Tau polarization (2-body modeled, 3-body flat)** — Ds -> tau nu and
  B+ -> tau nu produce ~100% longitudinally polarized taus. The 2-body
  hadronic modes (`tau -> {pi,K,rho,K*} N`) now use
  `kinematics.decay_2body_polarized`, sampling the N polar angle from
  `1 + (alpha·P_tau) cosθ` about the tau momentum axis. P_tau = -1 is exact
  for `P+ -> tau+ nu`; the analyzing power uses the chiral (maximal) limit
  `|alpha| = 1` (exact as m_N -> 0, an upper bound at finite m_N — see
  `TAU_2BODY_ANALYZING_POWER` in `production/decay_engine/tau_decay.py`). The 3-body
  leptonic modes remain flat phase space (their fully-correct treatment needs
  the decay matrix element), so a residual differential bias persists there.
- **B0 -> tau nu omitted** — helicity-suppressed in SM; only B+ -> tau nu
  retained in the induced-tau chain.
- **Prompt tau channel** (W -> tau nu, gamma*/Z -> tau tau) is now produced
  by `production/madgraph/run_tau_production.py` and lands in `tau/`; the
  earlier "dominant prompt-tau missing" caveat is retired. Its absolute
  rate is sub-dominant to induced tau at LHC -- charm production
  (sigma_ccbar ~ 2.7 mb) makes Ds -> tau nu the leading tau source by
  ~2 orders of magnitude over W -> tau nu (sigma_W ~ 0.1 mb x BR). Prompt
  tau still matters for completeness and for forward-rapidity / high-pT
  signatures.
- **Approximate kaon (K+ -> lN) flux** — the kaon channel is now included
  (`generate_kaon_csvs.py`), but FONLL supplies no light-meson spectrum, so
  the K± production uses a parametrized soft-QCD flux (Tsallis pT + Gaussian
  rapidity) and an order-of-magnitude inclusive cross-section
  (`SIGMA_KAON_PB`). The kaon *shape* and especially the *absolute
  normalization* below ~0.5 GeV are a systematic, not a precision input;
  regenerate from a measured K± spectrum to remove this caveat.
- **No Lambda_b / Lambda_c / Xi_c baryon channels** — pp fragmentation
  fractions now track the omitted baryon component explicitly, but only meson
  parents are simulated in this production layer.
