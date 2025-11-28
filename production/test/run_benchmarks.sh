#!/bin/bash
# ==========================================================================
# run_benchmarks.sh
#
# Run HNL production for Physics Beyond Colliders benchmark points
# BC6 (electron), BC7 (muon), BC8 (tau)
#
# Usage: ./run_benchmarks.sh [nEvents]
# ==========================================================================

set -e

# Number of events (default 100k, can override)
NEVENTS=${1:-100000}

# Mass points to scan (in GeV)
# Following MATHUSLA methodology
MASSES_LOW="0.1 0.2 0.3 0.4"                          # Kaon regime
MASSES_CHARM="0.5 0.75 1.0 1.25 1.5 1.75"            # D meson regime  
MASSES_BEAUTY="2.0 2.5 3.0 3.5 4.0 4.5"              # B meson regime
MASSES_EW="5.0 7.5 10.0 15.0 20.0 30.0 50.0"         # EW regime

# Combine all masses
ALL_MASSES="$MASSES_LOW $MASSES_CHARM $MASSES_BEAUTY $MASSES_EW"

# Flavors (PBC benchmarks)
FLAVORS="electron muon tau"

echo "============================================"
echo "HNL Production Benchmark Run"
echo "============================================"
echo "Events per point: $NEVENTS"
echo "Mass points: $(echo $ALL_MASSES | wc -w)"
echo "Flavors: $FLAVORS"
echo "Total jobs: $(( $(echo $ALL_MASSES | wc -w) * $(echo $FLAVORS | wc -w) ))"
echo "============================================"
echo ""

# Check if executable exists
if [ ! -f ./main_hnl_production ]; then
    echo "Error: main_hnl_production not found. Run 'make' first."
    exit 1
fi

# Create output directory
mkdir -p output
cd output

# Counter
NJOBS=0
NFAILED=0

# Loop over flavors
for FLAVOR in $FLAVORS; do
    echo ""
    echo "=== Running $FLAVOR (BC$([ "$FLAVOR" = "electron" ] && echo 6 || ([ "$FLAVOR" = "muon" ] && echo 7 || echo 8))) ==="
    
    # Loop over masses
    for MASS in $ALL_MASSES; do
        NJOBS=$((NJOBS + 1))
        echo -n "  Mass $MASS GeV... "
        
        # Run simulation
        if ../main_hnl_production $MASS $FLAVOR $NEVENTS > /dev/null 2>&1; then
            echo "OK"
        else
            echo "FAILED"
            NFAILED=$((NFAILED + 1))
        fi
    done
done

echo ""
echo "============================================"
echo "Summary"
echo "============================================"
echo "Total jobs: $NJOBS"
echo "Successful: $((NJOBS - NFAILED))"
echo "Failed: $NFAILED"
echo ""
echo "Output files in: $(pwd)"
ls -la *.csv 2>/dev/null | wc -l
echo " CSV files created"
echo "============================================"
