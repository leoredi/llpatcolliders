# HNL Production Pipeline

Full chain from production to money plot.

---

## Prerequisites

Activate conda environment first:
```fish
conda activate llpatcolliders
```

---

## 1. Pythia Production (Meson Decays)

```fish
cd production/pythia_production
./run_parallel_production.sh all
```

The script will show a summary and ask for confirmation before starting.

For single flavour:
```fish
./run_parallel_production.sh electron
./run_parallel_production.sh muon
./run_parallel_production.sh tau
```

Monitor progress (separate terminal):
```fish
watch -n 10 'ls ../../output/csv/simulation/HNL_*.csv | wc -l'
```

---

## 2. MadGraph Production (Electroweak)

Start Docker:
```fish
cd production/madgraph_production
docker run -it --rm -v (pwd)/../..:/work mg5-hnl:latest bash
```

Inside Docker:
```bash
cd /work/production/madgraph_production

# All flavours
python3 scripts/run_hnl_scan.py

# Or single flavour
python3 scripts/run_hnl_scan.py --flavour electron
python3 scripts/run_hnl_scan.py --flavour muon
python3 scripts/run_hnl_scan.py --flavour tau

exit
```

---

## 3. Combine Production Channels

```fish
cd analysis_pbc
python limits/combine_production_channels.py
```

For single flavour:
```fish
python limits/combine_production_channels.py --flavour tau
```

Check output:
```fish
ls ../output/csv/simulation/*_combined.csv | wc -l
```

---

## 4. Run Analysis

```fish
cd analysis_pbc
python limits/run.py --parallel
```

Or tau only (preserves existing e/mu results):
```fish
python run_tau_only.py --parallel
```

For batch jobs (disable progress bars):
```fish
python limits/run.py --parallel --no-progress
```

---

## 5. Money Plot

```fish
cd money_plot
python plot_money_island.py
open ../output/images/HNL_moneyplot_island.png
```

---

## Quick Reference

| Step | Directory | Command |
|------|-----------|---------|
| Pythia | `production/pythia_production` | `./run_parallel_production.sh all` |
| MadGraph | `production/madgraph_production` (Docker) | `python3 scripts/run_hnl_scan.py` |
| Combine | `analysis_pbc` | `python limits/combine_production_channels.py` |
| Analysis | `analysis_pbc` | `python limits/run.py --parallel` |
| Analysis (batch) | `analysis_pbc` | `python limits/run.py --parallel --no-progress` |
| Plot | `money_plot` | `python plot_money_island.py` |

---

## Troubleshooting

Pythia library error:
```fish
set -x DYLD_LIBRARY_PATH (pwd)/pythia8315/lib $DYLD_LIBRARY_PATH
```

Kill stuck Pythia jobs:
```fish
killall main_hnl_production
```

Check/stop Docker:
```fish
docker ps
docker stop (docker ps -q)
```
