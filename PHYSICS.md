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

The batch launcher emits `fromTau` only for `m_N < 1.77 GeV`.

### 3.2 Electroweak production (MadGraph)

For higher masses, the pipeline can add `ew` files from `W/Z -> N` production.
These are merged with meson channels at the same `(mass, flavour)`.

## 4. QCD slicing for transverse detectors

At O(100 m) transverse distance, acceptance is dominated by boosted HNLs.
To populate the high-`pT` tail efficiently, production uses four Pythia passes:

- `auto` (baseline SoftQCD pass).
- `hardccbar` with default `pTHatMin = 10 GeV`.
- `hardbbbar` with default `pTHatMin = 10 GeV`.
- `hardBc` with default `pTHatMin = 15 GeV` and parent filter `|parent_pdg| = 541`.

The combination stage keeps the best variant per regime/mode and concatenates regimes.

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

The runtime pipeline uses precomputed decay libraries and an HNLCalc interface for lifetime/BR inputs.

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

- `MASS_GRID`: `130` points, `0.20` to `17.00 GeV`.
- `N_EVENTS_DEFAULT`: pp collisions to simulate per production job (default `100k`).
- `MAX_SIGNAL_EVENTS`: max HNL signal events per channel (default `100k`). Production stops early when this is reached; analysis downsamples if exceeded.

Tau `direct` runs use the full grid. Tau `fromTau` is only generated below `1.77 GeV`.

## 9. Dominant approximations

Main projection-level approximations to keep in mind:

- Soft-QCD kaon normalisation is approximate.
- Some decay-kinematics inputs use nearest available precomputed decay samples.
- 3-body decay kinematics in production are simplified relative to full matrix-element treatments.

These are acceptable for fast sensitivity projections, but they are not a full experimental systematics model.
