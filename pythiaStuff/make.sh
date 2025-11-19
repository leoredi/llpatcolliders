#!/bin/bash

# HNL Mass Scan - Universal Script
# Usage: ./make.sh [lepton_flavor]
# lepton_flavor: electron, muon (default), tau, or 'all' for all three

# 1. SAFETY SETTINGS
set -e
cd "$(dirname "$0")"

# 2. CONFIGURATION
PYTHIA_DIR="/Users/fredi/cernbox/Physics/llpatcolliders/pythia-install"
CORES=10

# Parse command line argument
LEPTON="${1:-muon}"  # Default to muon if no argument

# 3. COMPILATION (ONCE)
echo "======================================================="
echo "HNL MASS SCAN"
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

# 5. RUN APPROPRIATE SCAN(S)
echo "-------------------------------------------------------"
echo "Step 3: Running Mass Scan..."
echo "-------------------------------------------------------"

case "$LEPTON" in
    electron|e)
        echo "Running ELECTRON scan..."
        LEPTON_NAME="electron"
        MASSES="0.2 0.3 0.4 0.5 1.0 1.5 2.0 2.5 2.8 3.1 3.3 3.6 4.0 4.4 4.8 5.2 10.0 15.0 20.0 40.0 80.0"
        echo "Mass points: $MASSES"
        echo $MASSES | xargs -n 1 -P $CORES -I {} sh -c "./main_hnl_single {} $LEPTON_NAME > logs/log_${LEPTON_NAME}_{}.txt 2>&1"
        echo "ELECTRON scan complete. Results: csv/HNL_mass_*_electron.csv"
        ;;

    muon|mu)
        echo "Running MUON scan..."
        LEPTON_NAME="muon"
        MASSES="0.5 1.0 1.5 2.0 2.5 2.8 3.1 3.3 3.6 4.0 4.4 4.8 5.2 10.0 15.0 20.0 40.0 80.0"
        echo "Mass points: $MASSES"
        echo $MASSES | xargs -n 1 -P $CORES -I {} sh -c "./main_hnl_single {} $LEPTON_NAME > logs/log_${LEPTON_NAME}_{}.txt 2>&1"
        echo "MUON scan complete. Results: csv/HNL_mass_*_muon.csv"
        ;;

    tau)
        echo "Running TAU scan..."
        LEPTON_NAME="tau"
        MASSES_LOW="0.5 1.0 1.5 2.0 2.5 2.8 3.1 3.4"
        MASSES_HIGH="10.0 15.0 20.0 40.0 80.0"
        MASSES="$MASSES_LOW $MASSES_HIGH"
        echo "Low Mass:  $MASSES_LOW"
        echo "High Mass: $MASSES_HIGH"
        echo "NOTE: Skipping 3.6-5.2 GeV (zero yield)"
        echo $MASSES | xargs -n 1 -P $CORES -I {} sh -c "./main_hnl_single {} $LEPTON_NAME > logs/log_${LEPTON_NAME}_{}.txt 2>&1"
        echo "TAU scan complete. Results: csv/HNL_mass_*_tau.csv"
        ;;

    all)
        echo "Running ALL lepton flavor scans (electron, muon, tau)..."
        echo ""

        # Electron
        echo "--- ELECTRON ---"
        LEPTON_NAME="electron"
        MASSES="0.2 0.3 0.4 0.5 1.0 1.5 2.0 2.5 2.8 3.1 3.3 3.6 4.0 4.4 4.8 5.2 10.0 15.0 20.0 40.0 80.0"
        echo $MASSES | xargs -n 1 -P $CORES -I {} sh -c "./main_hnl_single {} $LEPTON_NAME > logs/log_${LEPTON_NAME}_{}.txt 2>&1"
        echo "Electron complete."
        echo ""

        # Muon
        echo "--- MUON ---"
        LEPTON_NAME="muon"
        MASSES="0.5 1.0 1.5 2.0 2.5 2.8 3.1 3.3 3.6 4.0 4.4 4.8 5.2 10.0 15.0 20.0 40.0 80.0"
        echo $MASSES | xargs -n 1 -P $CORES -I {} sh -c "./main_hnl_single {} $LEPTON_NAME > logs/log_${LEPTON_NAME}_{}.txt 2>&1"
        echo "Muon complete."
        echo ""

        # Tau
        echo "--- TAU ---"
        LEPTON_NAME="tau"
        MASSES_LOW="0.5 1.0 1.5 2.0 2.5 2.8 3.1 3.4"
        MASSES_HIGH="10.0 15.0 20.0 40.0 80.0"
        MASSES="$MASSES_LOW $MASSES_HIGH"
        echo $MASSES | xargs -n 1 -P $CORES -I {} sh -c "./main_hnl_single {} $LEPTON_NAME > logs/log_${LEPTON_NAME}_{}.txt 2>&1"
        echo "Tau complete."
        echo ""

        echo "ALL scans complete."
        ;;

    *)
        echo "ERROR: Unknown lepton flavor '$LEPTON'"
        echo "Usage: ./make.sh [electron|muon|tau|all]"
        echo ""
        echo "Or use the dedicated scripts:"
        echo "  ./scan_electron.sh"
        echo "  ./scan_muon.sh"
        echo "  ./scan_tau.sh"
        exit 1
        ;;
esac

echo "======================================================="
echo "SCAN COMPLETE."
echo "======================================================="