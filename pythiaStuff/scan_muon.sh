#!/bin/bash

# Mass Scan Script for MUONS
# Standard mass range starting from muon mass threshold

# 1. SAFETY SETTINGS
set -e
cd "$(dirname "$0")"

# 2. CONFIGURATION
PYTHIA_DIR="/Users/fredi/cernbox/Physics/llpatcolliders/pythia-install"
LEPTON="muon"

# Muon mass range: Starting from ~0.5 GeV (above muon mass)
MASSES="0.5 1.0 1.5 2.0 2.5 2.8 3.1 3.3 3.6 4.0 4.4 4.8 5.2 10.0 15.0 20.0 40.0 80.0"

# Number of CPU cores to use
CORES=10

# 3. COMPILATION
echo "======================================================="
echo "MUON HNL MASS SCAN"
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
echo "Step 3: Launching Parallel Scan for MUONS on $CORES cores..."
echo "Mass points: $MASSES"
echo "-------------------------------------------------------"

# Run the scan with muon parameter
echo $MASSES | xargs -n 1 -P $CORES -I {} sh -c "./main_hnl_single {} $LEPTON > logs/log_muon_{}.txt 2>&1"

echo "======================================================="
echo "MUON SCAN COMPLETE."
echo "Results saved to csv/HNL_mass_*_muon.csv"
echo "======================================================="
