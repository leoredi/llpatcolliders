#!/bin/bash
# Complete HNL Analysis Pipeline
# This script runs the full workflow: PYTHIA simulation → ROOT to CSV → Analysis

set -e  # Exit on error

echo "=========================================="
echo "HNL Analysis Pipeline"
echo "=========================================="
echo ""

# Check if main144 exists
if [ ! -f "pythiaStuff/main144" ]; then
    echo "ERROR: pythiaStuff/main144 not found!"
    echo "Please compile PYTHIA code first:"
    echo "  cd pythiaStuff && ./make.sh"
    exit 1
fi

# Step 1: Generate HNL events with PYTHIA
echo "Step 1: Generating HNL events with PYTHIA..."
echo "Command: cd pythiaStuff && ./main144 -c hnlLL.cmnd"
echo ""
cd pythiaStuff
./main144 -c hnlLL.cmnd
echo ""
echo "✓ PYTHIA simulation complete"
echo ""

# Check if ROOT file was created
if [ ! -f "main144.root" ]; then
    echo "ERROR: main144.root not created!"
    exit 1
fi

# Step 2: Convert ROOT to CSV
echo "Step 2: Converting ROOT output to CSV..."
echo "Command: python convertRootToCsv_HNL.py main144.root ../LLP.csv"
echo ""
python convertRootToCsv_HNL.py main144.root ../LLP.csv
echo ""
echo "✓ ROOT to CSV conversion complete"
echo ""

# Return to main directory
cd ..

# Check if CSV file was created
if [ ! -f "LLP.csv" ]; then
    echo "ERROR: LLP.csv not created!"
    exit 1
fi

# Step 3: Run decay probability analysis
echo "Step 3: Running decay probability analysis..."
echo "Command: python decayProbPerEvent.py"
echo ""
python decayProbPerEvent.py
echo ""
echo "✓ Analysis complete"
echo ""

# Summary
echo "=========================================="
echo "HNL Analysis Complete!"
echo "=========================================="
echo ""
echo "Output files:"
echo "  - HNL_exclusion_vs_lifetime.png  (main exclusion plot)"
echo "  - HNL_particle_decay_results.csv (particle-level data)"
echo "  - HNL_event_decay_statistics.csv (event-level statistics)"
echo ""
echo "To view the exclusion plot:"
echo "  - Linux: xdg-open HNL_exclusion_vs_lifetime.png"
echo "  - macOS: open HNL_exclusion_vs_lifetime.png"
echo ""
