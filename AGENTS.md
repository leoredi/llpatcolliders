# AGENT CONTRACT (MACHINE-ORIENTED)

Scope: repository-local execution contract for AI agents.
Conflict rule: if this file conflicts with code, code wins.

## DOC_ROLE_MAP

- `README.md`: minimal run entrypoint only.
- `CODING.md`: implementation source of truth.
- `PHYSICS.md`: physics source of truth.
- `AGENTS.md`: agent execution constraints and fast lookup.

## CANONICAL_TERMS

- `mass point`: one value from `config_mass_grid.py::MASS_GRID`.
- `flavour`: one of `electron|muon|tau`.
- `production mode` (tau-only): `direct|fromTau`.
- `production regime`: `kaon|charm|beauty|Bc|ew|combined`.
- `QCD mode`: `auto|hardccbar|hardbbbar|hardBc`.
- `meta sidecar`: `<simulation_csv>.meta.json`.
- `analysis_pbc/hnl_models/`: analysis Python code modules.
- `analysis_pbc/model/`: HNLCalc runtime cache directory.

## DO_NOT_EDIT

- `production/pythia_production/pythia8315/`
- `production/madgraph_production/mg5/`
- `analysis_pbc/decay/external/`
- `analysis_pbc/HNLCalc/`

## ENTRYPOINT_INDEX

- mass grid: `config_mass_grid.py`
- pythia generator: `production/pythia_production/main_hnl_production.cc`
- pythia batch launcher: `production/pythia_production/run_parallel_production.sh`
- production combine: `analysis_pbc/limits/combine_production_channels.py`
- limit runner: `analysis_pbc/limits/run.py`
- signal kernel: `analysis_pbc/limits/expected_signal.py`
- hnlcalc scaling validator: `tools/analysis/check_hnlcalc_scaling.py`
- xsecs/frag config: `analysis_pbc/config/production_xsecs.py`
- money plot: `money_plot/plot_money_island.py`

## NON_MAIN_UTILITIES

- EW xsec validator: `tools/madgraph/validate_xsec.py`
- pythia monitor: `tools/pythia/monitor_production.sh`
- custom decay event generation: `tools/decay/generate_hnl_decay_events.py`

## RUNTIME_CONSTANTS (CODE-ANCHORED)

- luminosity: `L_HL_LHC_FB = 3000` (`analysis_pbc/limits/run.py`)
- exclusion threshold: `N_limit = 2.996` (`analysis_pbc/limits/expected_signal.py`)
- default separation: `--separation-mm 1.0` (`analysis_pbc/limits/run.py`)
- tau fromTau threshold: `1.77 GeV` (`production/pythia_production/run_parallel_production.sh`)
- FONLL/LHCb reference cross-sections (`analysis_pbc/config/production_xsecs.py`):
- `sigma(ccbar) = 23.6 mb`
- `sigma(bbbar) = 495 microbarn`
- `sigma(Bc) = 0.9 microbarn`

## PRODUCTION_PROTOCOL (TRANSVERSE DETECTOR)

Required Pythia passes:

```bash
cd production/pythia_production
PYTHIA8=$(pwd)/pythia8315 make
./run_parallel_production.sh all both
./run_parallel_production.sh all direct hardccbar 10
./run_parallel_production.sh all direct hardbbbar 10
./run_parallel_production.sh all direct hardBc 15
```

Optional EW pass (Docker image name: `mg5-hnl`):

```bash
docker run --rm -it -v "$(pwd):/work" mg5-hnl bash
# inside container: run production/madgraph_production/scripts/run_hnl_scan.py
```

Post-production:

```bash
cd analysis_pbc
python limits/combine_production_channels.py
python limits/run.py --parallel --workers 12
cd ../
python money_plot/plot_money_island.py
```

## FILE_DISCOVERY_RULES

- ignore simulation CSVs with size `<1000 bytes`.
- simulation filename grammar supports qcd suffixes:
- `..._<hardBc|hardccbar|hardbbbar>_pTHatX.csv`
- combined output: `HNL_<mass>GeV_<flavour>_combined.csv`

## NORMALIZATION_PATH

1. Pythia emits CSV + meta sidecar with `qcd_mode`, `sigma_gen_pb`, `pthat_min_gev`.
2. combine/run attach sidecar metadata to dataframes.
3. expected-signal uses per-parent cross-sections from `production_xsecs.py`.
4. limits are extracted from `N_sig(eps2)` crossings at `N=2.996`.

## SAFETY_RULES

- never apply absolute cross-section twice; CSV `weight` is relative MC weight.
- preserve meta sidecars when moving/copying simulation CSVs.
- on restart cleanup, delete `HNL_*.csv` and `HNL_*.csv.meta.json` together.
- canonical simulation location is `output/csv/simulation/`; files left in `production/pythia_production/` are staging leftovers.
- if geometry definition changes, invalidate `output/csv/geometry/` cache.
- use `abs(parent_pdg)` for parent BR/cross-section lookup.
- if `hardBc` stats are sparse, increase event count before interpreting gaps.

## VALIDATION_COMMANDS

```bash
python tools/docs/check_docs_sync.py
python tools/analysis/check_hnlcalc_scaling.py
python tools/madgraph/validate_xsec.py production/madgraph_production/summary_HNL_EW_production.csv
```
