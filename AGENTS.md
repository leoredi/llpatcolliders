# Agent Quick-Start

This document gets AI agents productive fast. For physics context read PHYSICS.md. For code architecture read CODING.md.

## What this repo does

Computes sensitivity projections for Heavy Neutral Leptons (HNLs) at a drainage-gallery detector near CMS at the HL-LHC. The pipeline: simulate HNL production from meson decays (Pythia) and EW production (MadGraph) → compute geometric acceptance in a tube-shaped detector → scan over coupling |U|² to find the exclusion region → plot.

## Do not touch

These directories contain vendored or external code. Never modify files inside them:

- `production/pythia_production/pythia8315/` — vendored Pythia 8.315
- `production/madgraph_production/mg5/` — vendored MadGraph5
- `analysis_pbc/decay/external/` — MATHUSLA RHN decay event files
- `analysis_pbc/HNLCalc/` — vendored HNLCalc decay width calculator

## Key entry points

| What you want to do | File |
|---|---|
| Change the mass grid | `config_mass_grid.py` |
| Change production cross-sections or fragmentation | `analysis_pbc/config/production_xsecs.py` |
| Change detector geometry | `analysis_pbc/geometry/per_parent_efficiency.py` |
| Change how signal events are counted | `analysis_pbc/limits/expected_signal.py` |
| Change the analysis orchestration | `analysis_pbc/limits/run.py` |
| Change the HNLCalc interface (not HNLCalc itself) | `analysis_pbc/models/hnl_model_hnlcalc.py` |
| Change decay acceptance logic | `analysis_pbc/decay/decay_detector.py` |
| Change production event generation | `production/pythia_production/main_hnl_production.cc` |
| Change the money plot | `money_plot/plot_money_island.py` |

## Data flow

```
config_mass_grid.py
       │
       ▼
Pythia (main_hnl_production.cc)    MadGraph (run_hnl_scan.py)
       │                                    │
       ▼                                    ▼
output/csv/simulation/HNL_*GeV_{electron,muon}_{kaon,charm,beauty}.csv
output/csv/simulation/HNL_*GeV_tau_{charm,beauty}_{direct,fromTau}.csv
output/csv/simulation/HNL_*GeV_*_ew.csv
       │
       ▼
combine_production_channels.py  →  HNL_*GeV_*_combined.csv
       │
       ▼
per_parent_efficiency.py  →  output/csv/geometry/  (cached ray-mesh intersections)
       │
       ▼
expected_signal.py  ←  hnl_model_hnlcalc.py (ctau, BRs)
                    ←  decay_detector.py (decay acceptance)
                    ←  production_xsecs.py (σ, fragmentation)
       │
       ▼
run.py  →  output/csv/analysis/HNL_U2_limits_summary.csv
       │
       ▼
plot_money_island.py  →  output/images/HNL_moneyplot_island.png
```

## CSV column conventions

**Production CSVs** (`output/csv/simulation/`):
`event, weight, hnl_id, parent_pdg, tau_parent_id, pt, eta, phi, p, E, mass, prod_x_mm, prod_y_mm, prod_z_mm, beta_gamma`

**Geometry cache** (`output/csv/geometry/`):
Adds: `hits_tube, entry_distance, path_length` (all in metres)

**Results CSV** (`output/csv/analysis/`):
`mass_GeV, flavour, benchmark, eps2_min, eps2_max, peak_events, separation_mm`

## Benchmark codes

- `"100"` = pure electron coupling (Ue² = eps2, Umu² = 0, Utau² = 0)
- `"010"` = pure muon coupling
- `"001"` = pure tau coupling

## Key physics constants in code

| Constant | Value | Location |
|---|---|---|
| HL-LHC luminosity | 3000 fb⁻¹ | `run.py` |
| Exclusion threshold | N = 2.996 events (95% CL Poisson) | `expected_signal.py` |
| Track separation cut | 1 mm default | `run.py` CLI arg |
| σ(bb̄) at 14 TeV | 5×10⁸ pb | `production_xsecs.py` |
| σ(cc̄) at 14 TeV | 2.4×10¹⁰ pb | `production_xsecs.py` |
| fromTau mass threshold | 1.77 GeV | `run_parallel_production.sh`, `main_hnl_production.cc` |

## Scaling optimisation

HNLCalc is slow (~0.5s per call). The key optimisation: ctau ∝ 1/|U|² and BR ∝ |U|². So we call HNLCalc once at a reference eps2, then scale analytically for all other eps2 values. This is validated by `test_scaling_vs_per_eps2.py`. The flag `--hnlcalc-per-eps2` disables this and recomputes from scratch (40x slower).

## Environment

```bash
conda activate llpatcolliders   # Python 3.11, numpy, pandas, matplotlib, trimesh, scipy
```

C++ compilation:
```bash
cd production/pythia_production && PYTHIA8=$(pwd)/pythia8315 make
```

## Running tests

```bash
cd analysis_pbc
python tests/test_scaling_vs_per_eps2.py
python scripts/check_hnlcalc_scaling.py
```

## Common pitfalls

1. Geometry cache staleness: if you change the detector geometry or CSV format, delete `output/csv/geometry/` to force recomputation.
2. The `parent_pdg` values can be signed in production CSVs; analysis uses `abs(parent_pdg)` for BR and cross-section lookups.
3. Tau production has two modes: "direct" (meson/W → τ N, mixing at the production vertex) and "fromTau" (meson/W → τ ν, then τ → N X; only for m_HNL < 1.77 GeV).
4. MadGraph EW production runs in Docker, not natively.
5. The `weight` column in production CSVs is the Pythia event weight (usually 1.0) or MadGraph cross-section weight.
