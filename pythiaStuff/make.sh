#!/bin/bash

# This line ensures the script will stop if any command fails.
set -e

# --- Define path to your Pythia installation ---
PYTHIA_DIR="../../pythia8315"

#-------------------------------------------------------
# STEP 1: Generate ROOT dictionary C++ code.
#-------------------------------------------------------
echo "Step 1: Generating ROOT dictionary with rootcling..."
rootcling -f main144Dct.cc -I${PYTHIA_DIR}/include main144Dct.h

#-------------------------------------------------------
# STEP 2: Compile dictionary into a shared library.
#-------------------------------------------------------
echo "Step 2: Compiling dictionary into main144Dct.so..."
g++ -shared -fPIC -o main144Dct.so main144Dct.cc \
    -I${PYTHIA_DIR}/include \
    `root-config --cflags --libs`

#-------------------------------------------------------
# STEP 3: Compile the main program and link everything.
#-------------------------------------------------------
echo "Step 3: Compiling main144.cc and linking..."
g++ main144.cc -o main144 \
    -I${PYTHIA_DIR}/include \
    -L${PYTHIA_DIR}/lib -Wl,-rpath,${PYTHIA_DIR}/lib -lpythia8 \
    main144Dct.so \
    `root-config --cflags --libs` \
    -DPY8ROOT

echo "---"
echo "Done! Executable 'main144' created successfully."