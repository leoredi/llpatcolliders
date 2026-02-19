# Coding Guide

This file describes how the pipeline is implemented.
For physics assumptions and equations, use `PHYSICS.md`.

## 1. Terminology (canonical)

Use these terms consistently across code and docs:

- `mass point`: one HNL mass on the unified `MASS_GRID`.
- `flavour`: one benchmark coupling (`electron`, `muon`, `tau`).
- `production mode` (tau only): `direct` or `fromTau`.
- `production regime`: `kaon`, `charm`, `beauty`, `Bc`, `ew`, `combined`.
- `QCD mode` (Pythia): `auto`, `hardccbar`, `hardbbbar`, `hardBc`.
- `meta sidecar`: `<simulation_csv>.meta.json` with run-level metadata.
- `analysis_pbc/hnl_models/`: analysis Python modules (code).
- `analysis_pbc/model/`: HNLCalc runtime cache (generated data; not analysis code).

## 2. Pipeline map

```text
config_mass_grid.py  (MASS_GRID, N_EVENTS_DEFAULT, MAX_SIGNAL_EVENTS)
  -> production/pythia_production/main_hnl_production.cc
  -> production/madgraph_production/scripts/run_hnl_scan.py (EW)
  -> analysis_pbc/limits/combine_production_channels.py
  -> analysis_pbc/limits/run.py
  -> money_plot/plot_money_island.py
```

Main outputs:

- `output/csv/simulation/` production CSVs.
- `output/csv/geometry/` geometry cache CSVs.
- `output/csv/analysis/HNL_U2_limits_summary.csv` final limits table.
- `output/images/HNL_moneyplot_island.png` final plot.

## 3. Full production recipe (transverse detector)

The full production chain consists of two main parts: Pythia production for meson-decay channels, and MadGraph production for high-mass Electroweak (EW) channels. Both `N_EVENTS_DEFAULT` and `MAX_SIGNAL_EVENTS` are read from `config_mass_grid.py` (`MAX_SIGNAL_EVENTS=0` means unlimited).

### Bash

**1. Pythia Production (Meson Decays)**

```bash
cd production/pythia_production
PYTHIA8=$(pwd)/pythia8315 make main_hnl_production

./run_parallel_production.sh all both auto
./run_parallel_production.sh all direct hardBc 15
cd ../..
```

**2. MadGraph Production (EW)**

```bash
docker run --rm -v "$(pwd):/work" mg5-hnl bash -c \
  "cd /work/production/madgraph_production && \
   python3 scripts/run_hnl_scan.py --flavour electron --min-mass 3 && \
   python3 scripts/run_hnl_scan.py --flavour muon --min-mass 3 && \
   python3 scripts/run_hnl_scan.py --flavour tau --min-mass 3"
```

**3. Analysis and Plotting**

```bash
cd analysis_pbc
python limits/combine_production_channels.py --allow-variant-drop
python limits/run.py --parallel --workers 12
cd ..
python money_plot/plot_money_island.py
```

### Fish

**1. Pythia Production (Meson Decays)**

```fish
set -gx REPO /path/to/llpatcolliders
cd $REPO

source (conda info --base)/etc/fish/conf.d/conda.fish
conda activate llpatcolliders

set -gx MPLCONFIGDIR /tmp/mpl-cache
mkdir -p $MPLCONFIGDIR

cd $REPO/production/pythia_production; or exit 1
set -lx PYTHIA8 (pwd)/pythia8315
make

cd $REPO/production/pythia_production; and begin
    ./run_parallel_production.sh all both auto
    and ./run_parallel_production.sh all direct hardBc 15
end
```

**2. MadGraph Production (EW)**

```fish
docker run --rm -v "$REPO:/work" mg5-hnl:latest /bin/bash -lc "
set -euo pipefail
cd /work/production/madgraph_production
python3 scripts/run_hnl_scan.py --flavour electron --min-mass 3
python3 scripts/run_hnl_scan.py --flavour muon --min-mass 3
python3 scripts/run_hnl_scan.py --flavour tau --min-mass 3
"
```

**3. Analysis and Plotting**

```fish
cd $REPO/analysis_pbc
python limits/combine_production_channels.py --allow-variant-drop
python limits/run.py --parallel --workers 12

cd $REPO
python money_plot/plot_money_island.py
```

### Notes

Only `auto` + `hardBc` Pythia passes are needed. The `auto` mode gives O(100k) HNL events per mass point throughout the kaon, charm, and beauty regimes. `hardBc` enriches Bc production for `m_N > 5 GeV` where auto yields only O(100) events.

`run_hnl_scan.py` flags: `--flavour electron|muon|tau`, `--min-mass <GeV>` (skip low masses where mesons dominate), `--masses <list>`, `--nevents <N>`, `--test` (single point: 15 GeV muon, 1k events). EW production is essential for masses above ~2 GeV.

Tau-only rerun does not need a dedicated script:

```bash
python limits/run.py --parallel --workers 12 --flavour tau
```

## 4. Production layer details

`production/pythia_production/main_hnl_production.cc`:

- CLI: `./main_hnl_production <mass_GeV> <flavor> [nEvents] [mode] [qcdMode] [pTHatMin] [maxSignalEvents]`.
- `qcdMode=auto` (default): inclusive SoftQCD or hard-QCD depending on mass. Sufficient for kaon, charm, and beauty regimes.
- `qcdMode=hardBc` uses the dedicated Bc card and filters to `|parent_pdg| = 541` in the event loop. Required for `m_N > 5 GeV` (Bc-only regime).
- Legacy modes `hardccbar` and `hardbbbar` are still supported but not needed (auto gives O(100k) HNLs throughout charm/beauty).
- `maxSignalEvents` stops the event loop early once enough HNLs are found (0 = unlimited).
- Every run writes a metadata sidecar containing `qcd_mode`, `sigma_gen_pb`, and `pthat_min_gev`.

`production/pythia_production/run_parallel_production.sh`:

- Batch launcher for full mass scans.
- Reads `N_EVENTS_DEFAULT` and `MAX_SIGNAL_EVENTS` from `config_mass_grid.py` via `load_mass_grid.sh`.
- Tau `fromTau` jobs are emitted only for `mass < 1.78 GeV`.
- In `auto` mode, tau `fromTau` is inclusive (one job per mass; no implicit `hardccbar/hardbbbar` split).

## 5. Combination and file selection

`analysis_pbc/limits/combine_production_channels.py`:

- Parses simulation filenames and groups by `(mass, flavour)`.
- Skips tiny files (`<1000` bytes).
- Nominal combine ignores `hardccbar` and `hardbbbar` samples.
- Chooses the best variant per `(regime, mode)` using priority:
  - `hardBc > auto` for `Bc`.
  - `_ff` over non-`_ff` when both exist.
  - higher `pTHat` wins when modes otherwise tie.
- Resolves overlapping channels by normalization keys (per-parent ownership), then filters each dataframe to owned keys.
- By default, keeps tau channels as explicit components (`direct`, `fromTau`, `ew`) and does not emit tau `_all.csv` unless `--allow-tau-all` is set.
- Writes `HNL_<mass>GeV_<flavour>_all.csv` for merged non-tau overlaps.

By default, the script raises on ambiguous same-regime variants (for example multiple `pTHat` slices) to prevent accidental drops. Pass `--allow-variant-drop` to warn and keep the top-priority variant.

`analysis_pbc/limits/run.py` applies the same nominal selection guards:

- skips `hardccbar` / `hardbbbar` files,
- rejects legacy tau `_all/_combined` inputs,
- resolves overlap ownership by parent normalization keys before computing acceptance.

## 6. Normalisation path

This is the critical implementation chain:

1. `main_hnl_production.cc` stores generator slice info in `.meta.json`.
2. `combine_production_channels.py` and `run.py` attach `qcd_mode`, `sigma_gen_pb`, `pthat_min_gev` to dataframes.
3. `expected_signal.py` computes signal using parent cross-sections from `analysis_pbc/config/production_xsecs.py` and decay/geometry acceptance.
4. `run.py` scans `|U|^2` and records exclusions at `N_sig >= 2.996` with `L = 3000 fb^-1`.

## 7. Cross-section source (code)

`analysis_pbc/config/production_xsecs.py` currently defines:

- `SIGMA_CCBAR_PB = 23.6 mb` (FONLL/LHCb 13/14 TeV reference).
- `SIGMA_BBBAR_PB = 495 microbarn` (FONLL/LHCb 13/14 TeV reference).
- `SIGMA_BC_PB = 0.9 microbarn` (FONLL/LHCb 13/14 TeV reference).

Treat this file as the normalisation source in code.

## 8. Output file contracts

Production CSV columns:

```text
event,weight,hnl_id,parent_pdg,tau_parent_id,pt,eta,phi,p,E,mass,prod_x_mm,prod_y_mm,prod_z_mm,beta_gamma
```

Geometry cache adds:

```text
hits_tube,entry_distance,path_length
```

Final limits CSV contains at minimum:

```text
mass_GeV,flavour,benchmark,eps2_min,eps2_max,peak_events,separation_mm
```

## 9. Cache and rerun rules

- If geometry definition or geometry cache schema changes, clear `output/csv/geometry/`.
- If production settings change (mass grid, QCD mode, pTHat, card-level content), regenerate affected simulation CSVs.
- Keep sidecar metadata files with their CSVs; deleting them weakens hard-QCD provenance.

Clean restart (recommended after interrupted or mixed production runs):

```bash
# Run from the repository root directory
cd /path/to/your/llpatcolliders

# stop in-flight jobs first (fish-safe)
pkill -f main_hnl_production; or true

# remove staging leftovers in production dir
find production/pythia_production -maxdepth 1 -name 'HNL_*.csv' -delete
find production/pythia_production -maxdepth 1 -name 'HNL_*.csv.meta.json' -delete

# remove canonical simulation outputs
find output/csv/simulation -maxdepth 1 -name 'HNL_*.csv' -delete
find output/csv/simulation -maxdepth 1 -name 'HNL_*.csv.meta.json' -delete

# remove downstream caches/results
find output/csv/geometry -maxdepth 1 -name 'HNL_*_geom*.csv' -delete
rm -f output/csv/analysis/HNL_U2_limits_summary.csv output/csv/analysis/HNL_U2_timing.csv
```

Important: delete CSV and `.meta.json` sidecars together.

### Decay overlay workflow

The runtime decay selector supports a generated overlay root while keeping
external repositories read-only:

- `output/decay/generated/` (overlay, highest priority)
- `analysis_pbc/decay/external/` (fallback)

Hybrid source-routing policy:

- `mass <= low_mass_threshold`: prefer external `analytical2and3bodydecays`.
- `low_mass_threshold < mass < 5.0 GeV`: use external hadronized files with legacy category priorities.
- `mass >= 5.0 GeV`: require generated overlay files (all-inclusive, category bypass).

Global mass grid (`config_mass_grid.py`) is `116` points from `0.20` to `10.00 GeV`.

Strict mismatch policy:

- Selection fails if `|m_requested - m_file| > 0.5 GeV`.
- Diagnostics override: `HNL_ALLOW_DECAY_MASS_MISMATCH=1` (warning instead of fail).


Limits/combination mass discovery is file-driven. When mass-grid boundaries
change, manually archive/remove stale `>10 GeV` simulation and generated decay
files before runs.

Parser note:

- External MATHUSLA decay files sometimes encode PID tokens as integral floats
  (for example `16.0`). `analysis_pbc/decay/rhn_decay_library.py::load_decay_events`
  normalizes integral float tokens to integer PDG IDs and warns when malformed
  rows are dropped.
- Do not strip neutrino species from generated overlay files as a compatibility
  workaround; fix parser/format issues at source instead.


For projection realism, set `--reco-efficiency` explicitly in limits runs
(for example `0.5`) instead of relying on the default idealized `1.0`.

## 10. HNLCalc scaling assumption

The fast limit scan uses the scaling relation `ctau ∝ 1/eps2` and `BR ∝ eps2`.
This is valid in the single-coupling benchmark regime used by this pipeline.
