# llpatcolliders

HNL sensitivity projections for a transverse LLP detector near CMS at the HL-LHC.

## Run

### Bash

```bash
conda activate llpatcolliders

# Pythia production
cd production/pythia_production
PYTHIA8=$(pwd)/pythia8315 make main_hnl_production
./run_parallel_production.sh all both
./run_parallel_production.sh all direct hardccbar 10
./run_parallel_production.sh all direct hardbbbar 10
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
    and ./run_parallel_production.sh all direct hardccbar 10
    and ./run_parallel_production.sh all direct hardbbbar 10
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

## Read next

- Physics assumptions and formulas: `PHYSICS.md`
- Code structure and execution details: `CODING.md`
