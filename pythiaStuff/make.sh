#!/bin/bash

# This line ensures the script will stop if any command fails.
set -e

# Change to the script's directory.
cd "$(dirname "$0")"

echo "-------------------------------------------------------"
echo "--- Pythia8 Not Found: Manual Installation Required ---"
echo "-------------------------------------------------------"
echo "The Pythia8 library is not installed in the expected location, and I am unable to install it automatically."
echo "Please follow these steps to install Pythia8 and compile the simulation code:"
echo ""
echo "1. Download Pythia8 (e.g., version 8.315):"
echo "   wget https://pythia.org/download/pythia83/pythia8315.tgz"
echo ""
echo "2. Extract the archive:"
echo "   tar xvfz pythia8315.tgz"
echo ""
echo "3. Configure and build Pythia8:"
echo "   cd pythia8315"
echo "   ./configure"
echo "   make"
echo ""
echo "4. Update the PYTHIA_DIR variable in this script (pythiaStuff/make.sh) to point to your Pythia8 installation directory."
echo ""
echo "5. Run this script again to compile the simulation."
echo ""
exit 1
