# Pipeline Runbook

Everything a human needs to run the full HNL sensitivity pipeline from scratch,
in copy-pasteable commands. Fish shell throughout; bash equivalents noted where
they differ.

---

## 0. Prerequisites

**Conda environment** (create once):

```fish
conda env create -f environment.yml
conda activate llpatcolliders
```

**Pythia 8** is vendored at `production/pythia_production/pythia8315/` — no
separate install needed.

**MadGraph / EW production** requires the `mg5-hnl` Docker image (built from
`production/madgraph_production/Dockerfile`). If you only need meson-regime
results (mass ≲ 5 GeV), skip the MadGraph step.

---

## 1. Build the Pythia executable

Run once per checkout (or after editing `main_hnl_production.cc`):

```fish
cd production/pythia_production
set -lx PYTHIA8 (pwd)/pythia8315
make
cd ../..
```

Bash equivalent:

```bash
cd production/pythia_production
PYTHIA8=$(pwd)/pythia8315 make
cd ../..
```

Quick sanity check (1 000 events, muon, 2 GeV):

```bash
cd production/pythia_production
./main_hnl_production 2.0 muon 1000
```

---

## 2. Pythia production (meson decays)

Two passes are needed: `auto` mode (covers kaon/charm/beauty regimes) and
`hardBc` mode (enriches Bc statistics for m_N > 5 GeV). Run both from
`production/pythia_production/`.

**Pass 1 — auto, all flavours, both modes:**

```fish
cd production/pythia_production
./run_parallel_production.sh all both auto
```

**Pass 2 — hardBc, all flavours, direct only, pTHatMin = 15 GeV:**

```fish
./run_parallel_production.sh all direct hardBc 15
cd ../..
```

### What this runs

- 116 mass points × 3 flavours × (direct + fromTau where applicable) = O(700)
  jobs in pass 1, O(350) jobs in pass 2.
- `fromTau` jobs are automatically skipped for `mass >= 1.776 GeV` (tau mass).
- Up to 12 jobs run in parallel (hardcoded in the script as `MAX_PARALLEL=12`).
- Each job simulates `N_EVENTS_DEFAULT = 100 000` pp collisions (set in
  `config_mass_grid.py`).

### Outputs

```
output/csv/simulation/HNL_<mass>GeV_<flavour>_<regime>[_<mode>][_<qcdMode>_pTHat<N>].csv
output/csv/simulation/HNL_*.csv.meta.json   ← keep alongside every CSV
output/logs/simulation/production_run_<timestamp>.log
```

### Monitoring a running pass

```fish
# live count of completed CSVs
watch "ls output/csv/simulation/HNL_*.csv 2>/dev/null | wc -l"

# tail the latest log
ls -t output/logs/simulation/*.log | head -1 | xargs tail -f

# check for failed jobs
grep -l "FAILED\|ERROR" output/logs/simulation/*.log
```

---

## 3. MadGraph production (EW, W/Z channels)

Required for m_N ≳ 5 GeV. Runs inside the `mg5-hnl` Docker image.

```bash
docker run --rm -v "$(pwd):/work" mg5-hnl bash -c \
  "cd /work/production/madgraph_production && \
   python3 scripts/run_hnl_scan.py --flavour electron --min-mass 3 && \
   python3 scripts/run_hnl_scan.py --flavour muon --min-mass 3 && \
   python3 scripts/run_hnl_scan.py --flavour tau --min-mass 3"
```

Fish variant (set `REPO` first):

```fish
set -gx REPO (pwd)
docker run --rm -v "$REPO:/work" mg5-hnl:latest /bin/bash -lc "
set -euo pipefail
cd /work/production/madgraph_production
python3 scripts/run_hnl_scan.py --flavour electron --min-mass 3
python3 scripts/run_hnl_scan.py --flavour muon    --min-mass 3
python3 scripts/run_hnl_scan.py --flavour tau     --min-mass 3
"
```

Available flags for `run_hnl_scan.py`:

| flag | default | meaning |
|---|---|---|
| `--flavour electron\|muon\|tau` | required | which coupling benchmark |
| `--min-mass <GeV>` | 0 | skip masses below this |
| `--masses <list>` | from mass grid | override specific mass points |
| `--nevents <N>` | from config | events per mass point |
| `--test` | off | single test point (15 GeV muon, 1 000 events) |

Outputs land in `output/csv/simulation/` as `HNL_*_ew.csv`.

---

## 4. Combine production channels

Merges per-regime files for each (mass, flavour), resolves parent ownership
across overlapping channels, and writes consolidated `_all.csv` files.

```fish
cd analysis_pbc
python limits/combine_production_channels.py --allow-variant-drop
cd ..
```

`--allow-variant-drop` is needed whenever both a `hardBc` and `auto` variant
exist for the same (mass, flavour, regime): the script keeps the higher-priority
one and discards the other.

### Useful flags

| flag | effect |
|---|---|
| `--dry-run` | print what would happen, write nothing |
| `--flavour electron\|muon\|tau` | process only one flavour |
| `--keep-originals` | don't delete per-regime source files after merging |
| `--allow-tau-all` | also emit a `_tau_all.csv` (off by default; tau components stay separate) |

### Outputs

```
output/csv/simulation/HNL_<mass>GeV_<flavour>_all.csv   ← electron & muon
output/csv/simulation/HNL_<mass>GeV_tau_direct.csv      ← tau components kept separate
output/csv/simulation/HNL_<mass>GeV_tau_fromTau.csv
output/csv/simulation/HNL_<mass>GeV_tau_ew.csv
```

---

## 5. Run the limits scan

Computes geometric acceptance (ray-tracing) and scans |U|² per mass point to
find exclusion intervals.

```fish
cd analysis_pbc
python limits/run.py --parallel --workers 12
cd ..
```

### Key flags

| flag | default | meaning |
|---|---|---|
| `--parallel` | off | enable multiprocessing |
| `--workers N` | all cores | number of parallel workers |
| `--mass <GeV>` | all | run a single mass point |
| `--max-mass <GeV>` | all | run only up to this mass |
| `--flavour electron\|muon\|tau` | all | run one flavour only |
| `--separation-mm <float>` | 1.0 | min charged-track separation at detector (mm) |
| `--p-min-gev <float>` | 0.6 | min charged-track momentum (GeV/c) |
| `--reco-efficiency <float>` | 1.0 | flat reconstruction efficiency; use 0.5 for realistic projections |
| `--dirac` | off | Dirac HNL interpretation (×2 yield vs Majorana) |
| `--decay-seed <int>` | 12345 | RNG seed for decay sampling |
| `--timing` | off | record per-mass timing; saved to `--timing-out` |
| `--allow-variant-drop` | off | suppress error on ambiguous QCD variants |

### Rerun for a single flavour

```bash
cd analysis_pbc
python limits/run.py --parallel --workers 12 --flavour tau
```

### Geometry cache

Geometry results are cached under `output/csv/geometry/`. If the detector
geometry changes, clear the cache before rerunning:

```bash
find output/csv/geometry -name 'HNL_*_geom*.csv' -delete
```

### Output

```
output/csv/analysis/HNL_U2_limits_summary.csv
```

Columns: `mass_GeV, flavour, benchmark, eps2_min, eps2_max, peak_events, separation_mm`

---

## 6. Make the money plot

```fish
python money_plot/plot_money_island.py
```

Reads `output/csv/analysis/HNL_U2_limits_summary.csv`, writes two files:

```
output/images/HNL_moneyplot_island.png                  ← fixed filename (overwritten each run)
output/images/<timestamp>_HNL_moneyplot_island.png      ← timestamped copy
```

---

## 7. Clean restart

Use this before re-running production from scratch, or after an interrupted run
that left mixed files:

```bash
# Kill any in-flight production jobs
pkill -f main_hnl_production; true

# Clear staging area inside production directory
find production/pythia_production -maxdepth 1 -name 'HNL_*.csv' -delete
find production/pythia_production -maxdepth 1 -name 'HNL_*.csv.meta.json' -delete

# Clear canonical simulation outputs
find output/csv/simulation -maxdepth 1 -name 'HNL_*.csv' -delete
find output/csv/simulation -maxdepth 1 -name 'HNL_*.csv.meta.json' -delete

# Clear downstream caches and results
find output/csv/geometry -maxdepth 1 -name 'HNL_*_geom*.csv' -delete
rm -f output/csv/analysis/HNL_U2_limits_summary.csv
rm -f output/csv/analysis/HNL_U2_timing.csv
```

> Always delete a CSV and its `.meta.json` sidecar together.

---

## 8. Full pipeline, one block

Assumes the Pythia binary is already built and the conda env is active.

```bash
# --- Pythia (meson) ---
cd production/pythia_production
./run_parallel_production.sh all both auto
./run_parallel_production.sh all direct hardBc 15
cd ../..

# --- MadGraph (EW) ---
docker run --rm -v "$(pwd):/work" mg5-hnl bash -c \
  "cd /work/production/madgraph_production && \
   python3 scripts/run_hnl_scan.py --flavour electron --min-mass 3 && \
   python3 scripts/run_hnl_scan.py --flavour muon --min-mass 3 && \
   python3 scripts/run_hnl_scan.py --flavour tau --min-mass 3"

# --- Analysis ---
cd analysis_pbc
python limits/combine_production_channels.py --allow-variant-drop
python limits/run.py --parallel --workers 12
cd ..

# --- Plot ---
python money_plot/plot_money_island.py
```

---

## 9. Output directory map

```
output/
├── csv/
│   ├── simulation/          # per-regime CSVs + .meta.json sidecars
│   │                        # + combined _all.csv / tau component files
│   ├── geometry/            # geometry cache (ray-trace hits); safe to delete
│   └── analysis/
│       ├── HNL_U2_limits_summary.csv   ← main result
│       └── HNL_U2_timing.csv           ← optional, from --timing
├── logs/
│   └── simulation/          # one log per production run
└── images/
    └── HNL_moneyplot_island.png
```

---

## 10. Physics quick-reference

| quantity | value |
|---|---|
| Luminosity | 3 000 fb⁻¹ (HL-LHC) |
| Exclusion threshold | N_sig ≥ 2.996 (95 % CL Poisson, zero background) |
| HNL type (default) | Majorana; use `--dirac` for ×2 |
| σ(cc̄) | 23.6 mb |
| σ(bb̄) | 495 μb |
| σ(Bc) | 0.9 μb |
| Default p_min | 0.6 GeV/c |
| Default separation | 1.0 mm |
| Tau fromTau ceiling | 1.776 GeV (τ mass) |
| Mass grid | 116 points, 0.20–10.00 GeV |

Cross-sections are defined in `analysis_pbc/config/production_xsecs.py`.
Pythia weights are kinematic only; all physical normalisation is applied in
the analysis step.
