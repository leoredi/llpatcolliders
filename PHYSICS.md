# Physics Guide

This file describes the physics model used in this repository.
For implementation details and run order, use `CODING.md`.

## 1. Goal

We project HL-LHC sensitivity to Heavy Neutral Leptons (HNLs) at a transverse LLP detector near CMS.

- Collision energy: `sqrt(s) = 14 TeV`.
- Integrated luminosity: `L = 3000 fb^-1`.
- Exclusion criterion: `N_sig >= 2.996` (95% CL, zero-background Poisson).

The output is an exclusion island in the `(m_N, |U_l|^2)` plane for each flavour.

## 2. Benchmarks and parameters

Single-flavour benchmarks:

- `electron` (`benchmark = "100"`): `|U_e|^2 = eps2`, others zero.
- `muon` (`benchmark = "010"`).
- `tau` (`benchmark = "001"`).

Per mass point, the scan runs over `eps2 = |U_l|^2` and extracts lower and upper crossings of `N_sig = 2.996`.

## 3. Production channels

### 3.1 Meson and baryon production (Pythia)

Production regimes used by the code:

- `kaon`.
- `charm`.
- `beauty`.
- `Bc`.

Tau has two production modes:

- `direct`: parent meson/boson produces `tau + N` at the primary decay.
- `fromTau`: parent produces a tau first, then `tau -> N X`.

The batch launcher emits `fromTau` only for `m_N < 1.78 GeV`.
In nominal `auto` runs this is an inclusive sample (no implicit hard-slice split).

### 3.2 Electroweak production (MadGraph)

For higher masses, the pipeline can add `ew` files from `W/Z -> N` production.
These are merged with meson channels at the same `(mass, flavour)`.

## 4. QCD production strategy

Standard PBC projections (MATHUSLA, ANUBIS, CODEX-b) use FONLL differential
cross-sections for meson kinematics. Our pipeline uses Pythia inclusive
generation with FONLL cross-sections for normalisation — a simpler approach
that gives equivalent results for meson-decay channels.

Production uses two Pythia passes:

- `auto`: inclusive SoftQCD (kaon regime) or hard-QCD (charm/beauty) depending
  on mass. Provides sufficient statistics across the full mass grid except
  for the Bc-only regime at high mass.
- `hardBc`: enriched bb̄ production (`HardQCD:gg2bbbar` + `qqbar2bbbar`,
  `pTHatMin = 15 GeV`) with parent filter `|parent_pdg| = 541`. Required
  for `m_N > 5 GeV` where Bc is the only meson parent and `auto` yields
  O(100) HNLs vs O(100k) at lower masses.

Legacy modes `hardccbar` and `hardbbbar` are still supported by the generator
for diagnostics, but nominal combine/limits runs ignore them.

The overlap stage resolves parent ownership by normalization keys so duplicated
parents are counted once. Tau channels are kept as explicit components
(`direct`, `fromTau`, `ew`) in nominal analysis; legacy tau `_all/_combined`
inputs are disabled by default.

## 5. Cross-section normalisation

The code normalises with FONLL/LHCb 13/14 TeV reference values in `analysis_pbc/config/production_xsecs.py`:

- `sigma(ccbar) = 23.6 mb`.
- `sigma(bbbar) = 495 microbarn`.
- `sigma(Bc) = 0.9 microbarn`.

Per-parent cross-sections are derived from these base values plus fragmentation fractions.

## 6. Decay and detector acceptance

For each produced HNL:

1. Geometry determines whether the trajectory intersects the detector and gives path length.
2. Decay probability in fiducial volume uses `beta*gamma`, path geometry, and `ctau`.
3. Decay visibility requires charged-track separation above the configured threshold (default `1 mm`).
   Baseline policy is `all-pairs-min` with no upper-separation cut.

The runtime pipeline uses precomputed decay libraries and an HNLCalc interface for lifetime/BR inputs.

Decay-library source priority:

- overlay (`output/decay/generated/`) first,
- external reference libraries (`analysis_pbc/decay/external/`) second.

Hybrid source routing is enforced:

- low-mass analytical regime (`mass <= low-mass threshold`): external analytical files.
- hadronized `mass < 5 GeV`: external MATHUSLA hadronized files.
- hadronized `mass >= 5 GeV`: generated overlay files (all-inclusive).

Parsing convention:

- External decay libraries may store PDG IDs as integral float text
  (for example `16.0`). The loader normalizes integral float PID tokens to
  integer IDs before applying visibility logic.
- Agreement fixes should preserve physics-complete generated content; neutrino
  stripping is not used as a compatibility patch.

Large decay-mass extrapolation is blocked by default (`|m_requested - m_file| > 0.5 GeV` fails).
The `4-5 GeV` overlap region is used to validate generated libraries against external references.

The `4-5 GeV` overlap region is used to validate generated libraries against external references.

## 7. Signal model

The expected signal is computed as a per-parent sum:

```text
N_sig = L * sum_parents [ sigma_parent * BR(parent -> N + X) * P_decay * epsilon_sep ]
```

Physical interpretation of exclusion island boundaries:

- Lower `eps2` boundary: HNL is too long-lived and decays after the detector.
- Upper `eps2` boundary: HNL is too short-lived and decays before reaching the detector.

## 8. Mass grid and event budget

`config_mass_grid.py` defines the shared configuration used by production and analysis:

- `MASS_GRID`: `116` points, `0.20` to `10.00 GeV`.
- `N_EVENTS_DEFAULT`: pp collisions to simulate per production job (default `100k`).
- `MAX_SIGNAL_EVENTS`: max HNL signal events per production job (`0` = unlimited, no early stop).

Tau `direct` runs use the full grid. Tau `fromTau` is only generated below `1.78 GeV`.

## 9. Dominant approximations

Main projection-level approximations to keep in mind:

- Soft-QCD kaon normalisation is approximate.
- Decay kinematics rely on precomputed libraries; coverage gaps now fail fast when
  the nearest available decay sample differs by more than `0.5 GeV` (unless explicitly overridden for diagnostics).
- 3-body decay kinematics in production are simplified relative to full matrix-element treatments.

These are acceptable for fast sensitivity projections, but they are not a full experimental systematics model.
