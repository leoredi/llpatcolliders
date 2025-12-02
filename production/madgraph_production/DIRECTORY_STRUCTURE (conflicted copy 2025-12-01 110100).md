# Directory Structure

## Current Layout (After Cleanup)

```
production/madgraph_production/
├── README.md                         # Main documentation
├── QUICKSTART.md                     # Quick user guide
├── EW_PRODUCTION_SUCCESS.md          # Production completion report
├── DOCKER_README.md                  # Docker-specific instructions
├── INSTALLATION.md                   # Setup guide
├── Dockerfile                        # Docker image definition
├── .gitignore                        # Git ignore rules
│
├── cards/                            # MadGraph input cards
│   ├── proc_card_electron.dat        # Process definition (e channel)
│   ├── proc_card_muon.dat            # Process definition (μ channel)
│   ├── proc_card_tau.dat             # Process definition (τ channel)
│   ├── run_card_template.dat         # Run parameters template
│   └── param_card_template.dat       # Model parameters template
│
├── scripts/                          # Production scripts
│   ├── run_hnl_scan.py               # Main driver (Docker/local)
│   └── lhe_to_csv.py                 # LHE → CSV converter
│
├── mg5/                              # MadGraph installation (147 MB)
│   └── bin/mg5_aMC                   # MadGraph executable
│
├── MG5_aMC_v3.6.6.tar.gz            # MadGraph source (30 MB)
│
├── work/                             # Process directories (cleaned)
├── lhe/                              # LHE output (cleaned)
│
└── archive/                          # Old files (402 MB)
    ├── csv_backup/                   # Old CSV location (now in output/)
    ├── *_production.log              # Completed run logs
    └── test_*.log                    # Test run logs
```

## Output Files (External)

Production output is saved to the central directory:

```
../../output/csv/simulation_new/
├── HNL_5p0GeV_electron_EW.csv       # EW production (5-80 GeV)
├── HNL_5p0GeV_muon_EW.csv
├── HNL_5p0GeV_tau_EW.csv
├── ...
├── HNL_80p0GeV_electron_EW.csv
├── HNL_80p0GeV_muon_EW.csv
├── HNL_80p0GeV_tau_EW.csv           # 96 files total
├── summary_HNL_EW_production.csv    # Cross-section metadata
│
├── HNL_0p20GeV_electron_kaon.csv    # Meson production (0.2-5 GeV)
├── ...                               # 152 files from Pythia
└── (267 total CSV files covering 0.2-80 GeV)
```

## Disk Usage

| Directory | Size | Purpose |
|-----------|------|---------|
| `mg5/` | 147 MB | MadGraph installation |
| `MG5_aMC_v3.6.6.tar.gz` | 30 MB | Source tarball |
| `archive/` | 402 MB | Old files (can be deleted) |
| `cards/` | 24 KB | Input templates |
| `scripts/` | 56 KB | Python scripts |
| `docs/` | ~100 KB | Documentation |
| **Total** | ~580 MB | (without archive: ~180 MB) |

## What Can Be Deleted?

### Safe to delete:
- `archive/` - Old production logs and CSV backup (402 MB)
- `work/` - Process directories (regenerated on each run)
- `lhe/` - LHE files (converted to CSV, regenerated if needed)

### Keep for future runs:
- `cards/` - Input templates (required)
- `scripts/` - Production scripts (required)
- `mg5/` - MadGraph installation (required if running locally)
- `MG5_aMC_v3.6.6.tar.gz` - Source (optional, for reproducibility)
- `Dockerfile` - Docker image definition (required for Docker runs)
- All `.md` files - Documentation

## Future Runs

To generate new mass points:

```bash
# Using Docker (recommended)
docker run --rm -v "$(pwd)":/work mg5-hnl \
  python3 scripts/run_hnl_scan.py --masses 85 90 --flavour muon

# Files will be saved to:
# ../../output/csv/simulation_new/HNL_85p0GeV_muon_EW.csv
# ../../output/csv/simulation_new/HNL_90p0GeV_muon_EW.csv
```

Work directories are automatically created/cleaned during each run.
