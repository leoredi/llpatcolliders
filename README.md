# HNL Sensitivity at the LHC — MATHUSLA Drainage Gallery

Sensitivity projections for Heavy Neutral Leptons (HNLs) at a proposed drainage-gallery detector near CMS at the HL-LHC. Covers electron, muon, and tau couplings over 0.2–17 GeV, combining meson/baryon production (Pythia), electroweak production (MadGraph), and a full decay + geometry acceptance pipeline.

## Repository layout

```
.
├── config_mass_grid.py                  # Central mass grid (130 points, 0.2–17 GeV)
├── environment.yml                      # Conda environment (llpatcolliders)
│
├── production/
│   ├── pythia_production/               # Meson/baryon → HNL (Pythia 8.315)
│   │   ├── main_hnl_production.cc       # C++ event generator
│   │   ├── Makefile                     # Build with Pythia8 linking
│   │   ├── run_parallel_production.sh   # Parallel job launcher
│   │   ├── load_mass_grid.sh            # Loads mass grid into bash
│   │   ├── monitor_production.sh        # Live dashboard
│   │   └── pythia8315/                  # Vendored Pythia (do not edit)
│   │
│   └── madgraph_production/             # EW production W/Z → ℓN (MadGraph5)
│       ├── MG5_aMC_v3.6.6.tar.gz        # Vendored MG5 source used by Docker build (keeps MG offline)
│       ├── scripts/
│       │   ├── run_hnl_scan.py          # MadGraph driver (Docker)
│       │   ├── lhe_to_csv.py            # LHE → CSV converter
│       │   └── validate_xsec.py         # Cross-section sanity checks
│       └── mg5/                         # Vendored MadGraph (do not edit)
│
├── analysis_pbc/
│   ├── HNLCalc/                         # Vendored HNLCalc (do not edit)
│   │   └── HNLCalc.py                   # Decay widths, BRs, ctau
│   ├── config/
│   │   └── production_xsecs.py          # σ(cc̄), σ(bb̄), fragmentation, τ BRs
│   ├── decay/
│   │   ├── external/                    # MATHUSLA RHN decay files (do not edit)
│   │   ├── decay_detector.py            # Decay acceptance + track separation
│   │   ├── rhn_decay_library.py         # Decay file selection/loading
│   │   └── generate_hnl_decay_events.py # MG5+Pythia decay sample generation
│   ├── geometry/
│   │   └── per_parent_efficiency.py     # Drainage gallery mesh + ray intersection
│   ├── limits/
│   │   ├── run.py                       # Main analysis orchestrator
│   │   ├── expected_signal.py           # N_sig computation + eps2 scanning
│   │   ├── combine_production_channels.py # Merge kaon/charm/beauty/ew CSVs
│   │   └── timing_utils.py              # Performance timing
│   ├── models/
│   │   └── hnl_model_hnlcalc.py         # Safe HNLCalc interface (ctau, BRs)
│   ├── scripts/
│   │   └── check_hnlcalc_scaling.py     # Validates BR ∝ |U|² scaling
│   ├── tests/
│   │   └── test_scaling_vs_per_eps2.py  # Unit test for scaling optimisation
│   └── run_tau_only.py                  # Re-run tau channel only
│
├── money_plot/
│   └── plot_money_island.py             # Exclusion island plots
│
└── output/                              # Generated at runtime
    ├── csv/simulation/                  # Production CSVs
    ├── csv/geometry/                    # Geometry cache
    ├── csv/analysis/                    # Limit results
    ├── logs/simulation/                 # Production logs
    └── images/                          # Plots
```

## Setup

```bash
conda env create -f environment.yml
conda activate llpatcolliders
```

Pythia 8.315 is vendored. Build the C++ generator:

```bash
cd production/pythia_production
PYTHIA8=$(pwd)/pythia8315 make
```
MadGraph5 is also vendored (tarball + `madgraph_production/mg5/`) and electroweak production must run inside the Docker image in `production/madgraph_production`. Build the image if you don't already have it locally (internet needed once for Pythia/LHAPDF downloads; MG5 itself comes from the tarball):

```bash
cd production/madgraph_production
docker build -t hnl-madgraph .
# Launch container with repo mounted at /work
docker run --rm -it -v "$(pwd)/../..:/work" hnl-madgraph bash
# Inside the container, switch to the repo
cd /work
```
If the `hnl-madgraph` image is already on your machine, skip the build and just `docker run` it.

## Run pipeline

The pipeline has four stages. Stage 2 (EW) is optional; Stage 3 is recommended when multiple production channels exist at the same mass/flavour. Stage 4 needs whatever production CSVs you produced in Stage 1/2.

### Stage 1 — Meson/baryon production (Pythia)

Generates HNL events from K, D, B meson and baryon decays via Pythia. Runs 12 parallel jobs, ~100k events per mass point.

```bash
cd production/pythia_production
./run_parallel_production.sh              # all flavours, direct + fromTau
./run_parallel_production.sh tau          # tau only
./run_parallel_production.sh muon direct  # muon, direct mode only
```

Monitor with:

```bash
watch -n 10 ./monitor_production.sh
```

Output: `output/csv/simulation/HNL_<mass>GeV_<flavour>_<regime>[_fromTau].csv`

### Stage 2 — Electroweak production (MadGraph, optional)

Generates W/Z → ℓN events for high-mass HNLs. **Run this stage from inside the Docker container launched above** (MadGraph + LHAPDF live in `/opt` inside the image; running on the host will fail).

```bash
# Inside container, working directory /work/production/madgraph_production
cd /work/production/madgraph_production
python3 scripts/run_hnl_scan.py --flavour muon --min-mass 3
python3 scripts/run_hnl_scan.py --test    # quick smoke test
```

Output: `output/csv/simulation/HNL_<mass>GeV_<flavour>_ew.csv`

### Stage 3 — Combine production channels (optional but recommended with multiple sources)

Merges kaon/charm/beauty/ew CSVs at each mass point into a single `..._combined.csv`. Prefers form-factor (`_ff`) files when available. `analysis_pbc/limits/run.py` can also read separate files directly, but combining keeps geometry caching and bookkeeping simpler.

```bash
cd analysis_pbc
python limits/combine_production_channels.py --dry-run   # preview
python limits/combine_production_channels.py              # merge
```

### Stage 4 — Limits and plotting

Runs the full analysis: geometry preprocessing, decay acceptance, eps2 scan, limit extraction.

```bash
cd analysis_pbc
python limits/run.py --parallel --workers 12

# or single flavour:
python limits/run.py --parallel --flavour tau

# or single mass:
python limits/run.py --mass 2.6 --flavour muon
```

Key flags:

| Flag | Description |
|------|-------------|
| `--parallel` | Multi-core processing |
| `--workers N` | Number of cores |
| `--separation-mm 1.0` | Charged track separation cut (mm) |
| `--dirac` | 2x yield for Dirac HNL |
| `--decay-seed 12345` | RNG seed for decay sampling |
| `--timing` | Per-mass timing breakdown |

Output: `output/csv/analysis/HNL_U2_limits_summary.csv`

Generate exclusion plots:

```bash
python money_plot/plot_money_island.py
```

Output: `output/images/HNL_moneyplot_island.png`

## Directories not to modify

These contain vendored or external code and must not be edited:

- `production/pythia_production/pythia8315/`
- `production/madgraph_production/mg5/`
- `analysis_pbc/decay/external/`
- `analysis_pbc/HNLCalc/`
