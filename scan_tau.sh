#!/bin/bash

# Mass Scan Script for TAUS
# Special handling: Skip intermediate masses (3.6-5.2 GeV) where yield is zero

# 1. SAFETY SETTINGS
set -e
cd "$(dirname "$0")"

# 2. CONFIGURATION
PYTHIA_DIR="/Users/fredi/cernbox/Physics/llpatcolliders/pythia-install"
LEPTON="tau"

# Tau mass range - TWO REGIMES:
# Low Mass (B-meson decays): 0.5 - 3.4 GeV (stop before zero-yield region)
# High Mass (W-decays): 10.0+ GeV (resume after zero-yield region)
# Skip: 3.6 4.0 4.4 4.8 5.2 (Yield will be zero due to tau mass threshold)

MASSES_LOW="0.5 1.0 1.5 2.0 2.5 2.8 3.1 3.4"
MASSES_HIGH="10.0 15.0 20.0 40.0 60.0"
MASSES="$MASSES_LOW $MASSES_HIGH"

# Number of CPU cores to use
CORES=10

# 3. COMPILATION
echo "======================================================="
echo "TAU HNL MASS SCAN"
echo "======================================================="
echo "Step 1: Compiling main_hnl_single.cc..."
echo "-------------------------------------------------------"

g++ main_hnl_single.cc -o main_hnl_single \
    -I${PYTHIA_DIR}/include \
    -L${PYTHIA_DIR}/lib -Wl,-rpath,${PYTHIA_DIR}/lib -lpythia8 \
    -std=c++17 \
    -O2

echo "Compilation successful."

# 4. DIRECTORY SETUP
echo "-------------------------------------------------------"
echo "Step 2: Checking directories..."
echo "-------------------------------------------------------"

mkdir -p csv
mkdir -p logs
echo "Directories 'csv' and 'logs' are ready."

# 5. PARALLEL EXECUTION
echo "-------------------------------------------------------"
echo "Step 3: Launching Parallel Scan for TAUS on $CORES cores..."
echo "-------------------------------------------------------"
echo "Low Mass Regime:  $MASSES_LOW"
echo "High Mass Regime: $MASSES_HIGH"
echo "NOTE: Skipping 3.6-5.2 GeV (zero yield due to tau mass)"
echo "-------------------------------------------------------"

# Run the scan with tau parameter
echo $MASSES | xargs -n 1 -P $CORES -I {} sh -c "./main_hnl_single {} $LEPTON > logs/log_tau_{}.txt 2>&1"

echo "======================================================="
echo "TAU SCAN COMPLETE."
echo "Results saved to csv/HNL_mass_*_tau.csv"
echo "======================================================="
