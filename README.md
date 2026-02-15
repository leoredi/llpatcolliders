# llpatcolliders

HNL sensitivity projections for a transverse LLP detector near CMS at the HL-LHC.

## Run

### Bash

```bash
conda activate llpatcolliders

# Pythia production (auto + hardBc only)
cd production/pythia_production
PYTHIA8=$(pwd)/pythia8315 make main_hnl_production
./run_parallel_production.sh all both auto
./run_parallel_production.sh all direct hardBc 15
cd ../..

# EW production (inside Docker image: mg5-hnl)
docker run --rm -v "$(pwd):/work" mg5-hnl bash -c \
  "cd /work/production/madgraph_production && \
   python3 scripts/run_hnl_scan.py --flavour electron --min-mass 3 && \
   python3 scripts/run_hnl_scan.py --flavour muon --min-mass 3 && \
   python3 scripts/run_hnl_scan.py --flavour tau --min-mass 3"

# Analysis
cd analysis_pbc
python limits/combine_production_channels.py --allow-variant-drop
python limits/run.py --parallel --workers 12
cd ..

# Plotting
python money_plot/plot_money_island.py
```

### Fish

```fish
set -gx REPO /path/to/llpatcolliders
cd $REPO

source (conda info --base)/etc/fish/conf.d/conda.fish
conda activate llpatcolliders

# Optional: silence matplotlib cache warnings
set -gx MPLCONFIGDIR /tmp/mpl-cache
mkdir -p $MPLCONFIGDIR

# Pythia production
cd $REPO/production/pythia_production; or exit 1
set -lx PYTHIA8 (pwd)/pythia8315
make

cd $REPO/production/pythia_production; and begin
    ./run_parallel_production.sh all both auto
    and ./run_parallel_production.sh all direct hardBc 15
end

# EW production (inside Docker image: mg5-hnl)
docker run --rm -v "$REPO:/work" mg5-hnl:latest /bin/bash -lc "
set -euo pipefail
cd /work/production/madgraph_production
python3 scripts/run_hnl_scan.py --flavour electron --min-mass 3
python3 scripts/run_hnl_scan.py --flavour muon --min-mass 3
python3 scripts/run_hnl_scan.py --flavour tau --min-mass 3
"

# Analysis
cd $REPO/analysis_pbc
python limits/combine_production_channels.py --allow-variant-drop
python limits/run.py --parallel --workers 12

# Plotting
cd $REPO
python money_plot/plot_money_island.py
```

## Decay Library Overlay Workflow

Decay files in `analysis_pbc/decay/external/` are treated as read-only.
Generated decay libraries should be written under:

- `output/decay/generated/`

Hybrid routing policy:

- low-mass analytical regime (`mass <= low-mass threshold`): external analytical files.
- hadronized region below `5 GeV`: external MATHUSLA files.
- hadronized region at/above `5 GeV`: generated overlay files.

Global mass grid (`config_mass_grid.py`) is `116` points from `0.20` to `10.00 GeV`.

Generate overlay libraries from `>= 4 GeV` (validation overlap + high-mass coverage)
with flavour-pure couplings (`|U|^2=1e-3` by default, one non-zero flavour at a time):

```bash
python tools/decay/precompute_decay_library_overlay.py \
  --flavours electron,muon,tau \
  --from-mass-grid \
  --overlay-min-mass 4.0 \
  --u2-norm 1e-3 \
  --nevents 20000
```

If Docker runs with `--user` and you see `/.config` or matplotlib cache
permission errors, use this container-safe wrapper:

```bash
docker run --rm -it \
  --user "$(id -u):$(id -g)" \
  -e HOME=/tmp \
  -e XDG_CONFIG_HOME=/tmp/.config \
  -e XDG_CACHE_HOME=/tmp/.cache \
  -e MPLCONFIGDIR=/tmp/.config/matplotlib \
  -v "/Users/fredi/sandbox-offline/llpatcolliders:/work" \
  mg5-hnl:latest \
  /bin/bash -lc '
    set -euo pipefail
    mkdir -p "$XDG_CONFIG_HOME" "$XDG_CACHE_HOME" "$MPLCONFIGDIR"
    cd /work
    python3 tools/decay/precompute_decay_library_overlay.py \
      --flavours electron,muon,tau \
      --from-mass-grid \
      --overlay-min-mass 4.0 \
      --u2-norm 1e-3 \
      --nevents 20000
  '
```

Notes:

- precompute overwrite is enabled by default; use `--no-overwrite` to keep existing files.
- external decay files may encode PID values as integral floats (for example `16.0`); `load_decay_events()` normalizes these to integer PDG IDs.
- generated overlay libraries are kept physics-complete; no neutrino species are intentionally stripped to force overlap agreement.
- MG5 warnings about missing `fastjet-config`, `lhapdf-config`, or
  `eMELA-config` are expected for this decay-at-rest workflow and are not
  fatal by themselves.
- because limits discovery is file-driven, manually archive/remove stale `>10 GeV`
  simulation and generated decay files from old grids before production runs.

Validate the 4-5 GeV overlap between generated and external libraries:

```bash
python tools/decay/validate_decay_overlap.py \
  --flavours electron,muon,tau \
  --from-mass-grid \
  --min-mass 4.0 \
  --max-mass 5.0 \
  --out output/decay/overlap_validation.csv
```

Audit coverage and strict mass matching:

```bash
python tools/decay/audit_decay_coverage.py \
  --flavours electron,muon,tau \
  --from-mass-grid \
  --overlay-switch-mass 5.0 \
  --out output/decay/coverage_report.csv
```

## Calibrated Fast Decay Mode (`brvis-kappa`)

For faster limits scans, a calibrated analytical surrogate is available:

- `P_decay` from geometry + lifetime (same as baseline),
- multiplied by `BR_vis` from HNLCalc,
- multiplied by calibrated `kappa(mass, flavour)`.

Calibration convention is fixed to:

- `p_min = 0.6 GeV`
- `separation = 1.0 mm`

Build dense-mass calibration and report:

```bash
python tools/decay/calibrate_brvis_kappa.py \
  --flavours electron,muon,tau \
  --from-mass-grid \
  --p-min-gev 0.6 \
  --separation-mm 1.0 \
  --out output/csv/analysis/decay_kappa_table.csv \
  --report output/csv/analysis/decay_kappa_validation.csv
```

Validate calibrated mode against library mode:

```bash
python tools/decay/validate_brvis_kappa.py \
  --flavours electron,muon,tau \
  --from-mass-grid \
  --p-min-gev 0.6 \
  --separation-mm 1.0 \
  --kappa-table output/csv/analysis/decay_kappa_table.csv \
  --out output/csv/analysis/decay_kappa_validation_check.csv
```

Run limits with the calibrated mode:

```bash
cd analysis_pbc
python limits/run.py --parallel --workers 12 \
  --decay-mode brvis-kappa \
  --kappa-table ../output/csv/analysis/decay_kappa_table.csv \
  --p-min-gev 0.6 \
  --separation-mm 1.0
```

By default, decay-file selection fails when `|m_requested - m_file| > 0.5 GeV`.
For diagnostics only, downgrade failures to warnings with:

```bash
# bash/zsh
export HNL_ALLOW_DECAY_MASS_MISMATCH=1
unset HNL_ALLOW_DECAY_MASS_MISMATCH
```

```fish
# fish
set -gx HNL_ALLOW_DECAY_MASS_MISMATCH 1
set -e HNL_ALLOW_DECAY_MASS_MISMATCH
```

## Decay Physics QA Gates

Use these three checks as the decay sign-off gate before final limits:

1. Generated vs external overlap in `4 <= m < 5 GeV`:

```bash
python tools/decay/validate_decay_overlap.py \
  --flavours electron,muon,tau \
  --from-mass-grid \
  --min-mass 4.0 \
  --max-mass 5.0 \
  --out output/decay/overlap_check_now.csv
```

2. Hybrid source-policy + strict mass-match coverage:

```bash
python tools/decay/audit_decay_coverage.py \
  --flavours electron,muon,tau \
  --from-mass-grid \
  --overlay-switch-mass 5.0 \
  --out output/decay/coverage_check_now.csv
```

3. Fast surrogate consistency (`brvis-kappa`) against `library` mode:

```bash
python tools/decay/validate_brvis_kappa.py \
  --flavours electron,muon,tau \
  --from-mass-grid \
  --p-min-gev 0.6 \
  --separation-mm 1.0 \
  --kappa-table output/csv/analysis/decay_kappa_table.csv \
  --out output/csv/analysis/decay_kappa_validation_check.csv
```

For a PBC-style projection, run limits with an explicit reconstruction
efficiency (for example `--reco-efficiency 0.5`) rather than relying on the
default idealized `1.0`.

## Read next

- Physics assumptions and formulas: `PHYSICS.md`
- Code structure and execution details: `CODING.md`
