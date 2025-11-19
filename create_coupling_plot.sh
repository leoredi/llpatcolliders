#!/bin/bash
#
# Script to create HNL coupling vs mass "money plot"
# This runs the complete pipeline from BR vs lifetime data to coupling limits
#
# Usage:
#   bash create_coupling_plot.sh [mu|tau]

SCENARIO=${1:-mu}  # Default to muon-coupled

echo "=========================================="
echo "HNL Coupling Limit Analysis"
echo "Scenario: ${SCENARIO}-coupled"
echo "=========================================="
echo ""

# Mass points to analyze (GeV)
# Only include masses with meaningful sensitivity (BR_limit < 1)
# For m >= 47 GeV, BR_limits are >1, indicating essentially zero acceptance
MASSES=(15 23 31 39)

# Step 1: Run decayProbPerEvent.py for each mass (if not already done)
echo "Step 1: Checking for exclusion data..."
for mass in "${MASSES[@]}"; do
    if [ "$SCENARIO" = "mu" ]; then
        CSV_FILE="output/csv/hnlLL_m${mass}GeVLLP.csv"
        EXCL_FILE="output/csv/hnlLL_m${mass}GeVLLP_exclusion_data.csv"
    else
        CSV_FILE="output/csv/hnlTauLL_m${mass}GeVLLP.csv"
        EXCL_FILE="output/csv/hnlTauLL_m${mass}GeVLLP_exclusion_data.csv"
    fi

    if [ ! -f "$EXCL_FILE" ]; then
        echo "  Missing: $EXCL_FILE"
        if [ -f "$CSV_FILE" ]; then
            echo "  Running analysis for m=${mass} GeV..."
            conda run -n llpatcolliders python decayProbPerEvent.py "$CSV_FILE"
        else
            echo "  ERROR: CSV file not found: $CSV_FILE"
            echo "  Run simulation first: cd pythiaStuff && python run_mass_scan.py --scenario $SCENARIO"
        fi
    else
        echo "  ✓ Found: $EXCL_FILE"
    fi
done

echo ""
echo "Step 2: Computing coupling limits from BR vs lifetime data..."
conda run -n llpatcolliders python hnl_coupling_limit.py --scenario "$SCENARIO"

echo ""
echo "=========================================="
echo "Analysis complete!"
echo "Check: output/images/hnl_coupling_vs_mass_${SCENARIO}.png"
echo "=========================================="
