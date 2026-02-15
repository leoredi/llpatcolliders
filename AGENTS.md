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
- madgraph EW driver: `production/madgraph_production/scripts/run_hnl_scan.py`
- production combine: `analysis_pbc/limits/combine_production_channels.py`
- limit runner: `analysis_pbc/limits/run.py`
- signal kernel: `analysis_pbc/limits/expected_signal.py`
- hnlcalc scaling validator: `tools/analysis/check_hnlcalc_scaling.py`
- xsecs/frag config: `analysis_pbc/config/production_xsecs.py`
- money plot: `money_plot/plot_money_island.py`

## NON_MAIN_UTILITIES

- production input tests: `tools/tests/production/test_production_inputs.py`
- EW xsec validator: `tools/madgraph/validate_xsec.py`
- LHE→CSV converter: `production/madgraph_production/scripts/lhe_to_csv.py`
- pythia monitor: `tools/pythia/monitor_production.sh`
- custom decay event generation: `tools/decay/generate_hnl_decay_events.py`
- decay overlay precompute: `tools/decay/precompute_decay_library_overlay.py`
- decay coverage audit: `tools/decay/audit_decay_coverage.py`
- decay overlap validator: `tools/decay/validate_decay_overlap.py`
- brvis-kappa calibrator: `tools/decay/calibrate_brvis_kappa.py`
- brvis-kappa validator: `tools/decay/validate_brvis_kappa.py`

## RUNTIME_CONSTANTS (CODE-ANCHORED)

Central production config (`config_mass_grid.py`):
- `MASS_GRID = 116 points`: `0.20` to `10.00 GeV`.
- `N_EVENTS_DEFAULT = 100_000`: pp collisions to simulate per production job.
- `MAX_SIGNAL_EVENTS = 1_000`: max HNL signal events per channel (caps production early + analysis downsampling).

Other constants:
- luminosity: `L_HL_LHC_FB = 3000` (`analysis_pbc/limits/run.py`)
- exclusion threshold: `N_limit = 2.996` (`analysis_pbc/limits/expected_signal.py`)
- default separation: `--separation-mm 1.0` (`analysis_pbc/limits/run.py`)
- default p-min: `--p-min-gev 0.6` (`analysis_pbc/limits/run.py`)
- decay modes: `library|brvis-kappa` (`analysis_pbc/limits/run.py`)
- tau fromTau threshold: `1.77 GeV` (`production/pythia_production/run_parallel_production.sh`)
- FONLL/LHCb reference cross-sections (`analysis_pbc/config/production_xsecs.py`):
- `sigma(ccbar) = 23.6 mb`
- `sigma(bbbar) = 495 microbarn`
- `sigma(Bc) = 0.9 microbarn`

## PRODUCTION_PROTOCOL (TRANSVERSE DETECTOR)

The production sequence has two main phases, followed by analysis. Both `N_EVENTS_DEFAULT` and `MAX_SIGNAL_EVENTS` are read from `config_mass_grid.py`.

**1. Pythia Production (auto + hardBc only)**
```bash
cd production/pythia_production
PYTHIA8=$(pwd)/pythia8315 make main_hnl_production
./run_parallel_production.sh all both auto
./run_parallel_production.sh all direct hardBc 15
cd ../..
```

**2. Electroweak (EW) Production**

EW production is run from the project root using the `mg5-hnl` Docker image.
```bash
docker run --rm -v "$(pwd):/work" mg5-hnl bash -c \
  "cd /work/production/madgraph_production && \
   python3 scripts/run_hnl_scan.py --flavour electron --min-mass 3 && \
   python3 scripts/run_hnl_scan.py --flavour muon --min-mass 3 && \
   python3 scripts/run_hnl_scan.py --flavour tau --min-mass 3"
# flags: --flavour <e|mu|tau>, --min-mass <GeV>, --masses <list>, --nevents <N>, --test (single point: 15 GeV muon, 1k events)
```

**3. Post-production**
```bash
cd analysis_pbc
python limits/combine_production_channels.py --allow-variant-drop
python limits/run.py --parallel --workers 12
cd ..
python money_plot/plot_money_island.py
```

## FILE_DISCOVERY_RULES

- ignore simulation CSVs with size `<1000 bytes`.
- simulation filename grammar supports qcd suffixes:
- `..._<hardBc|hardccbar|hardbbbar>_pTHatX.csv`
- combined output: `HNL_<mass>GeV_<flavour>_all.csv`

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
python tools/decay/validate_decay_overlap.py --flavours electron,muon,tau --from-mass-grid --min-mass 4.0 --max-mass 5.0 --out output/decay/overlap_check_now.csv
python tools/decay/audit_decay_coverage.py --flavours electron,muon,tau --from-mass-grid --overlay-switch-mass 5.0 --out output/decay/coverage_check_now.csv
python tools/decay/validate_brvis_kappa.py --flavours electron,muon,tau --from-mass-grid --p-min-gev 0.6 --separation-mm 1.0 --kappa-table output/csv/analysis/decay_kappa_table.csv --out output/csv/analysis/decay_kappa_validation_check.csv
pytest tools/tests/production/test_production_inputs.py -v
# after EW production completes:
python tools/madgraph/validate_xsec.py production/madgraph_production/summary_HNL_EW_production.csv
```
