#!/bin/bash
# Run decay probability analysis for multiple mass points in parallel
# Usage: ./run_decay_analysis.sh [muon|electron|both]

LEPTON="${1:-muon}"

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate llpatcolliders

# Select mass points based on lepton
case "$LEPTON" in
    muon)
        echo "Running decay analysis for MUON mass points..."
        MASSES="10.0 15.0 20.0 40.0"
        LEPTON_NAME="muon"
        ;;
    electron)
        echo "Running decay analysis for ELECTRON mass points..."
        MASSES="10.0 15.0 20.0 40.0"
        LEPTON_NAME="electron"
        ;;
    both)
        echo "Running decay analysis for BOTH muon and electron..."
        # Run muon first
        $0 muon
        # Then electron
        $0 electron
        exit 0
        ;;
    *)
        echo "Usage: $0 [muon|electron|both]"
        exit 1
        ;;
esac

# Create output directories
mkdir -p output/csv output/images

# Run analyses in parallel (2 at a time to avoid overload)
echo "Processing masses: $MASSES"
echo $MASSES | xargs -n 1 -P 2 -I {} bash -c "
    echo 'Starting analysis for mass={} GeV (${LEPTON_NAME})...'
    python decayProbPerEvent.py csv/HNL_mass_{}_${LEPTON_NAME}.csv > logs/decay_analysis_{}GeV_${LEPTON_NAME}.log 2>&1
    echo 'Completed analysis for mass={} GeV (${LEPTON_NAME})'
"

echo "All decay analyses complete for ${LEPTON_NAME}!"
echo "Results saved to output/csv/"
