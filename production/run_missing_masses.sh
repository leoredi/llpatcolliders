#!/bin/bash
# Run missing mass points from dense grid scan

set -e
cd "$(dirname "$0")"

echo "======================================================="
echo "RUNNING MISSING MASS POINTS"
echo "======================================================="
echo ""

# Missing masses: 0.25, 0.35, 0.45 for electron and muon
MISSING_MASSES="0.25 0.35 0.45"

echo "Missing masses identified: $MISSING_MASSES"
echo "Flavors: electron, muon"
echo ""

# Electron
echo "--- ELECTRON ---"
for mass in $MISSING_MASSES; do
    echo "Running m = $mass GeV (electron)..."
    ./main_hnl_single $mass electron > ../output/logs/simulation/log_electron_${mass}.txt 2>&1
    echo "  Done: HNL_mass_${mass}_electron_*.csv"
done

echo ""

# Muon
echo "--- MUON ---"
for mass in $MISSING_MASSES; do
    echo "Running m = $mass GeV (muon)..."
    ./main_hnl_single $mass muon > ../output/logs/simulation/log_muon_${mass}.txt 2>&1
    echo "  Done: HNL_mass_${mass}_muon_*.csv"
done

echo ""
echo "======================================================="
echo "MISSING MASSES COMPLETE"
echo "======================================================="
echo ""
echo "Summary:"
echo "  - Added 6 mass points (3 electron + 3 muon)"
echo "  - New total: 102 mass points"
echo ""
echo "Verify:"
ls -1 ../output/csv/simulation/*_electron_*.csv 2>/dev/null | wc -l | xargs -I {} echo "  Electron files: {}"
ls -1 ../output/csv/simulation/*_muon_*.csv 2>/dev/null | wc -l | xargs -I {} echo "  Muon files: {}"
ls -1 ../output/csv/simulation/*_tau_*.csv 2>/dev/null | wc -l | xargs -I {} echo "  Tau files: {}"
