#!/bin/bash

# This line ensures the script will stop if any command fails.
set -e

# Change to the script's directory.
cd "$(dirname "$0")"

# --- Define path to your Pythia installation ---
PYTHIA_DIR="/tmp/pythia8315/pythia8315"

#-------------------------------------------------------
# STEP 1: Compile the main program and link everything.
#-------------------------------------------------------
echo "Step 1: Compiling main144.cc and linking..."
g++ main144.cc -o main144 \
    -I${PYTHIA_DIR}/include \
    -L${PYTHIA_DIR}/lib -Wl,-rpath,${PYTHIA_DIR}/lib -lpythia8

echo "---"
echo "Done! Executable 'main144' created successfully."
