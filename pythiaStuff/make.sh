#!/bin/bash

# This line ensures the script will stop if any command fails.
set -e

# Change to the script's directory.
cd "$(dirname "$0")"

# --- Define path to your Pythia installation ---
PYTHIA_DIR="/Users/fredi/cernbox/Physics/llpatcolliders/pythia-install"

#-------------------------------------------------------
# STEP 1: Compile the main program and link everything.
#-------------------------------------------------------
echo "Step 1: Compiling main_hnl_scan.cc and linking..."
g++ main_hnl_scan.cc -o main_hnl_scan \
    -I${PYTHIA_DIR}/include \
    -L${PYTHIA_DIR}/lib -Wl,-rpath,${PYTHIA_DIR}/lib -lpythia8

echo "---"
echo "Done! Executable 'main_hnl_scan' created successfully."
