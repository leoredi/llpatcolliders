# MadGraph HNL Production - Docker Workflow

## Overview

MadGraph HNL production runs inside a Docker container with all dependencies (MadGraph5, LHAPDF6, Pythia8) pre-installed.

## Prerequisites

- Docker installed and running
- MadGraph tarball: `MG5_aMC_v3.6.6.tar.gz` (should be in this directory)

## Build Docker Image

```bash
cd /path/to/llpatcolliders/production/madgraph

# Build the image (first time only, ~5-10 minutes)
docker build -t mg5-hnl .
```

The Dockerfile:
- Installs Ubuntu 22.04 base
- Installs MadGraph5 v3.6.6
- Installs LHAPDF6 and Pythia8 via MadGraph
- Sets up Python environment

## Run Production

### Interactive Shell (Recommended for First Run)

```bash
# Start Docker container with current directory mounted at /work
docker run --rm -it -v "$(pwd)":/work mg5-hnl bash

# Inside container, you're in /work (which is your madgraph directory)
# Run test production
python3 scripts/run_hnl_scan.py --test

# Or run full muon production
python3 scripts/run_hnl_scan.py --flavour muon

# Exit when done
exit
```

### One-Line Production (Non-Interactive)

```bash
# Test run (15 GeV muon, 1000 events)
docker run --rm -v "$(pwd)":/work mg5-hnl python3 scripts/run_hnl_scan.py --test

# Full muon production (32 mass points, 50k events each)
docker run --rm -v "$(pwd)":/work mg5-hnl python3 scripts/run_hnl_scan.py --flavour muon

# Custom mass points
docker run --rm -v "$(pwd)":/work mg5-hnl python3 scripts/run_hnl_scan.py --flavour muon --masses 10 15 20
```

## Output

All output is written to the mounted directory, so it persists after the container exits:

```
/work/
├── csv/
│   ├── HNL_15p0GeV_muon_EW.csv       # Event data (Pythia CSV format)
│   ├── HNL_20p0GeV_muon_EW.csv
│   └── summary_HNL_EW_production.csv  # Cross-sections + metadata
├── work/
│   ├── hnl_muon_15.0GeV/              # MadGraph process directory
│   │   ├── Events/                     # LHE files
│   │   ├── Cards/                      # run_card, param_card
│   │   └── generate_events.log
│   └── ...
└── test_run.log                        # Optional log file
```

## CSV Output Format

The output CSVs match the Pythia meson production format:

```csv
event,weight,hnl_id,parent_pdg,pt,eta,phi,p,E,mass,prod_x_mm,prod_y_mm,prod_z_mm,boost_gamma
1,1,9900012,24,5.23,1.45,-0.82,15.6,16.4,15.0,0,0,0,1.04
2,1,9900012,-24,8.91,-2.31,2.15,42.3,43.9,15.0,0,0,0,2.82
...
```

Where:
- `hnl_id`: Always 9900012 (N1 PDG code)
- `parent_pdg`: 24 (W+), -24 (W-), 23 (Z), or 0 (if W/Z off-shell and not in LHE)
- `prod_x_mm`, `prod_y_mm`, `prod_z_mm`: Always (0,0,0) for EW production (at IP)
- `boost_gamma`: β γ = p / mass

## Troubleshooting

### Docker image not found
```bash
docker build -t mg5-hnl .
```

### Permission denied accessing /work
```bash
# Add your user to docker group (Linux)
sudo usermod -aG docker $USER
# Then log out and back in
```

### Container runs but no output
Check that you're mounting the correct directory:
```bash
pwd  # Should be .../production/madgraph
docker run --rm -v "$(pwd)":/work mg5-hnl ls -la /work
# Should show scripts/, cards/, mg5/, etc.
```

### LHAPDF errors inside container
This shouldn't happen if the image was built correctly. Rebuild:
```bash
docker rmi mg5-hnl
docker build -t mg5-hnl .
```

## Performance Notes

- **Test run (15 GeV, 1k events):** ~2-5 minutes
- **Single mass point (50k events):** ~10-20 minutes
- **Full muon production (32 points):** ~6-10 hours

The bottleneck is MadGraph's matrix element calculation and integration.

## Validation

After production, validate cross-sections:

```bash
# Inside or outside container
python3 scripts/validate_xsec.py csv/summary_HNL_EW_production.csv
```

Expected: σ ~ 10,000-18,000 pb for m=15 GeV with |U|²=1.0

## Environment Variables

Inside the container:
- `MG5_DIR=/opt/MG5_aMC_v3_6_6`
- `MG5_PATH=/opt/MG5_aMC_v3_6_6/bin/mg5_aMC`
- `PATH` includes MadGraph bin directory

## Next Steps

After generating EW production CSVs:
1. Combine with meson production CSVs (< 5 GeV)
2. Run geometry analysis to compute acceptance
3. Run limit calculation with HNLCalc
4. Generate exclusion plots

See main `CLAUDE.md` for full pipeline documentation.
