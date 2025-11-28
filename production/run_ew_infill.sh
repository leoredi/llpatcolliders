#!/bin/bash
# Generate only the missing EW regime mass points (11-40 GeV)
# Adds 13 new points for better granularity

set -e

NEVENTS=200000
OUTPUT_DIR="../output/csv/simulation_new"
LOG_DIR="../output/logs/simulation_new"

mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_DIR"

# Only the NEW mass points we're adding
NEW_EW_MASSES=(11.0 13.0 14.0 16.0 17.0 18.0 19.0 22.0 25.0 28.0 32.0 35.0 38.0 40.0)

echo "=========================================="
echo "EW Regime Infill Production"
echo "=========================================="
echo "Adding ${#NEW_EW_MASSES[@]} new mass points per flavor"
echo "Total simulations: $((${#NEW_EW_MASSES[@]} * 3)) (electron + muon + tau)"
echo "Events per mass: $NEVENTS"
echo ""

TOTAL_SIMS=$((${#NEW_EW_MASSES[@]} * 3))
CURRENT=0

# Electron
for MASS in "${NEW_EW_MASSES[@]}"; do
    CURRENT=$((CURRENT + 1))
    MASS_STR=$(printf "%.2f" $MASS | sed 's/\./p/g')
    echo "[$CURRENT/$TOTAL_SIMS] Generating electron, m=$MASS GeV..."
    ./main_hnl_production $MASS electron $NEVENTS > "${LOG_DIR}/HNL_${MASS_STR}GeV_electron_ew.log" 2>&1
done

# Muon
for MASS in "${NEW_EW_MASSES[@]}"; do
    CURRENT=$((CURRENT + 1))
    MASS_STR=$(printf "%.2f" $MASS | sed 's/\./p/g')
    echo "[$CURRENT/$TOTAL_SIMS] Generating muon, m=$MASS GeV..."
    ./main_hnl_production $MASS muon $NEVENTS > "${LOG_DIR}/HNL_${MASS_STR}GeV_muon_ew.log" 2>&1
done

# Tau
for MASS in "${NEW_EW_MASSES[@]}"; do
    CURRENT=$((CURRENT + 1))
    MASS_STR=$(printf "%.2f" $MASS | sed 's/\./p/g')
    echo "[$CURRENT/$TOTAL_SIMS] Generating tau, m=$MASS GeV..."
    ./main_hnl_production $MASS tau $NEVENTS > "${LOG_DIR}/HNL_${MASS_STR}GeV_tau_ew.log" 2>&1
done

echo ""
echo "EW infill complete!"
echo "Generated $TOTAL_SIMS new simulation files"
