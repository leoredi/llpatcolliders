# Full HNL Production Guide

Step-by-step bash commands to run the complete production pipeline.

---

## Prerequisites

```bash
# Activate conda environment
conda activate llpatcolliders
```

---

## Part 1: Pythia Production (Meson Decays)

Generates HNL events from meson decays (K, D, B mesons).

```bash
# Navigate to pythia production directory
cd production/pythia_production
```

```bash
# Check that the executable exists
ls -lh main_hnl_production
```

```bash
# Run full production (all flavours: electron, muon, tau)
./run_parallel_production.sh all
```

Or run individual flavours:

```bash
# Electron only
./run_parallel_production.sh electron
```

```bash
# Muon only
./run_parallel_production.sh muon
```

```bash
# Tau only
./run_parallel_production.sh tau
```

### Monitor Progress (in separate terminal)

```bash
# Watch running jobs
watch -n 5 "ps aux | grep main_hnl_production | grep -v grep | wc -l"
```

```bash
# Watch CSV files being created
watch -n 10 "ls ../../output/csv/simulation/HNL_*.csv 2>/dev/null | wc -l"
```

### Verify Pythia Output

```bash
# Count output files
ls ../../output/csv/simulation/HNL_*.csv | wc -l
```

```bash
# Check for failures
grep -l "FAILED\|ERROR" ../../output/logs/simulation/*.log
```

---

## Part 2: MadGraph Production (Electroweak)

Generates HNL events from W/Z boson decays (high mass regime 5-80 GeV).

### Step 2a: Build Docker Image

```bash
# Navigate to madgraph production directory
cd production/madgraph_production
```

```bash
# Build the Docker image (first time only, takes ~10 min)
docker build -t mg5-hnl:latest .
```

### Step 2b: Run MadGraph Production

```bash
# Start interactive Docker container with repo mounted
docker run -it --rm -v "$(pwd)/../..:/work" mg5-hnl:latest bash
```

Inside the Docker container:

```bash
# Navigate to madgraph production directory
cd /work/production/madgraph_production
```

```bash
# Test run (single mass point, 1000 events)
python3 scripts/run_hnl_scan.py --test
```

```bash
# Full production - all flavours (takes several hours)
python3 scripts/run_hnl_scan.py
```

Or run individual flavours:

```bash
# Electron only
python3 scripts/run_hnl_scan.py --flavour electron
```

```bash
# Muon only
python3 scripts/run_hnl_scan.py --flavour muon
```

```bash
# Tau only
python3 scripts/run_hnl_scan.py --flavour tau
```

```bash
# Custom mass points
python3 scripts/run_hnl_scan.py --masses 10 15 20 --flavour muon
```

```bash
# Exit Docker when done
exit
```

### Verify MadGraph Output

```bash
# Back on host machine
ls output/csv/simulation/HNL_*_ew.csv | wc -l
```

```bash
# Check summary file
cat output/csv/simulation/summary_HNL_ew_production.csv
```

---

## Part 3: Post-Production

### Combine Overlapping Channels

If both meson and EW files exist at overlapping masses (4-8 GeV):

```bash
cd analysis_pbc
```

```bash
# Dry run to see what would be combined
python limits/combine_production_channels.py --dry-run
```

```bash
# Actually combine
python limits/combine_production_channels.py
```

### Run Analysis

```bash
cd analysis_pbc
```

```bash
python limits/run.py
```

### Generate Plots

```bash
cd money_plot
```

```bash
python plot_money_island.py
```

---

## Quick Reference

| Step | Directory | Command |
|------|-----------|---------|
| Pythia | `production/pythia_production` | `./run_parallel_production.sh all` |
| MadGraph | `production/madgraph_production` (Docker) | `python3 scripts/run_hnl_scan.py` |
| Combine | `analysis_pbc` | `python limits/combine_production_channels.py` |
| Analysis | `analysis_pbc` | `python limits/run.py` |
| Plots | `money_plot` | `python plot_money_island.py` |

---

## Troubleshooting

### Pythia: Library not found

```bash
export DYLD_LIBRARY_PATH="$(pwd)/pythia8315/lib:$DYLD_LIBRARY_PATH"
```

### MadGraph: Docker not running

```bash
# Check Docker is running
docker info
```

### Kill stuck jobs

```bash
# Pythia jobs
killall main_hnl_production

# Docker container
docker ps
docker stop <container_id>
```
