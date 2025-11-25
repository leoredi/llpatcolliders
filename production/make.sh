#!/bin/bash

# HNL Mass Scan - Universal Script
# Usage: ./make.sh [lepton_flavor]
# lepton_flavor: electron, muon (default), tau, or 'all' for all three

# 1. SAFETY SETTINGS
set -e
cd "$(dirname "$0")"

# 2. CONFIGURATION
PYTHIA_RELDIR="pythia/pythia8315"
CORES=10

# Parse command line argument
LEPTON="${1:-muon}"  # Default to muon if no argument

# 2.5 AUTO-COMPILE PYTHIA IF NEEDED
if [ ! -f "$PYTHIA_RELDIR/lib/libpythia8.a" ]; then
    echo "======================================================="
    echo "PYTHIA NOT COMPILED - Building Pythia 8.315..."
    echo "======================================================="
    cd "$PYTHIA_RELDIR"
    ./configure --prefix=$(pwd)
    make -j4
    cd ../..
    echo "Pythia compilation complete."
    echo ""
fi

# Export absolute path for compiler
export PYTHIA_DIR="$(pwd)/$PYTHIA_RELDIR"

# 3. COMPILATION (ONCE)
echo "======================================================="
echo "HNL MASS SCAN"
echo "======================================================="
echo "Step 1: Compiling main_hnl_single.cc..."
echo "-------------------------------------------------------"

g++ main_hnl_single.cc -o main_hnl_single \
    -I$PYTHIA_DIR/include \
    -L$PYTHIA_DIR/lib -Wl,-rpath,$PYTHIA_DIR/lib -lpythia8 \
    -std=c++17 \
    -O2

echo "Compilation successful."

# 4. DIRECTORY SETUP
echo "-------------------------------------------------------"
echo "Step 2: Checking directories..."
echo "-------------------------------------------------------"

mkdir -p ../output/csv/simulation
mkdir -p ../output/logs/simulation
echo "Directories '../output/csv/simulation' and '../output/logs/simulation' are ready."

# 5. RUN APPROPRIATE SCAN(S)
echo "-------------------------------------------------------"
echo "Step 3: Running Mass Scan..."
echo "-------------------------------------------------------"

case "$LEPTON" in
    electron|e)
        echo "Running ELECTRON scan..."
        LEPTON_NAME="electron"
        MASSES="0.2 0.25 0.3 0.35 0.4 0.45 0.5 \
                0.6 0.7 0.8 0.9 1.0 \
                1.2 1.4 1.6 1.8 \
                2.0 2.3 2.6 3.0 3.4 \
                3.8 4.2 4.6 \
                5.0 6.0 7.0 8.0 9.0 \
                10.0 12.0 15.0 20.0 30.0 40.0 50.0 60.0 80.0"
        echo "Mass points: $MASSES"
        echo $MASSES | xargs -n 1 -P $CORES -I {} sh -c "./main_hnl_single {} $LEPTON_NAME > ../output/logs/simulation/log_${LEPTON_NAME}_{}.txt 2>&1"
        echo "ELECTRON scan complete. Results: ../output/csv/simulation/HNL_mass_*_electron_*.csv"
        ;;

    muon|mu)
        echo "Running MUON scan..."
        LEPTON_NAME="muon"
        MASSES="0.2 0.25 0.3 0.35 0.4 0.45 0.5 \
                0.6 0.7 0.8 0.9 1.0 \
                1.2 1.4 1.6 1.8 \
                2.0 2.3 2.6 3.0 3.4 \
                3.8 4.2 4.6 \
                5.0 6.0 7.0 8.0 9.0 \
                10.0 12.0 15.0 20.0 30.0 40.0 50.0 60.0 80.0"
        echo "Mass points: $MASSES"
        echo $MASSES | xargs -n 1 -P $CORES -I {} sh -c "./main_hnl_single {} $LEPTON_NAME > ../output/logs/simulation/log_${LEPTON_NAME}_{}.txt 2>&1"
        echo "MUON scan complete. Results: ../output/csv/simulation/HNL_mass_*_muon_*.csv"
        ;;

    tau)
        echo "Running TAU scan..."
        LEPTON_NAME="tau"
        MASSES="0.5 0.7 1.0 1.3 1.6 2.0 2.4 2.8 3.2 3.6 4.0 4.5 \
                5.0 6.0 7.0 8.0 9.0 \
                10.0 12.0 15.0 20.0 30.0 40.0 50.0 60.0 80.0"
        echo "Mass points: $MASSES"
        echo $MASSES | xargs -n 1 -P $CORES -I {} sh -c "./main_hnl_single {} $LEPTON_NAME > ../output/logs/simulation/log_${LEPTON_NAME}_{}.txt 2>&1"
        echo "TAU scan complete. Results: ../output/csv/simulation/HNL_mass_*_tau_*.csv"
        ;;

    all)
        echo "Running ALL lepton flavor scans (electron, muon, tau)..."
        echo ""

        # Electron
        echo "--- ELECTRON ---"
        LEPTON_NAME="electron"
        MASSES="0.2 0.25 0.3 0.35 0.4 0.45 0.5 \
                0.6 0.7 0.8 0.9 1.0 \
                1.2 1.4 1.6 1.8 \
                2.0 2.3 2.6 3.0 3.4 \
                3.8 4.2 4.6 \
                5.0 6.0 7.0 8.0 9.0 \
                10.0 12.0 15.0 20.0 30.0 40.0 50.0 60.0 80.0"
        echo $MASSES | xargs -n 1 -P $CORES -I {} sh -c "./main_hnl_single {} $LEPTON_NAME > ../output/logs/simulation/log_${LEPTON_NAME}_{}.txt 2>&1"
        echo "Electron complete."
        echo ""

        # Muon
        echo "--- MUON ---"
        LEPTON_NAME="muon"
        MASSES="0.2 0.25 0.3 0.35 0.4 0.45 0.5 \
                0.6 0.7 0.8 0.9 1.0 \
                1.2 1.4 1.6 1.8 \
                2.0 2.3 2.6 3.0 3.4 \
                3.8 4.2 4.6 \
                5.0 6.0 7.0 8.0 9.0 \
                10.0 12.0 15.0 20.0 30.0 40.0 50.0 60.0 80.0"
        echo $MASSES | xargs -n 1 -P $CORES -I {} sh -c "./main_hnl_single {} $LEPTON_NAME > ../output/logs/simulation/log_${LEPTON_NAME}_{}.txt 2>&1"
        echo "Muon complete."
        echo ""

        # Tau
        echo "--- TAU ---"
        LEPTON_NAME="tau"
        MASSES="0.5 0.7 1.0 1.3 1.6 2.0 2.4 2.8 3.2 3.6 4.0 4.5 \
                5.0 6.0 7.0 8.0 9.0 \
                10.0 12.0 15.0 20.0 30.0 40.0 50.0 60.0 80.0"
        echo $MASSES | xargs -n 1 -P $CORES -I {} sh -c "./main_hnl_single {} $LEPTON_NAME > ../output/logs/simulation/log_${LEPTON_NAME}_{}.txt 2>&1"
        echo "Tau complete."
        echo ""

        echo "ALL scans complete."
        ;;

    *)
        echo "ERROR: Unknown lepton flavor '$LEPTON'"
        echo "Usage: ./make.sh [electron|muon|tau|all]"
        echo ""
        exit 1
        ;;
esac

echo "======================================================="
echo "SCAN COMPLETE."
echo "======================================================="