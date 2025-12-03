#!/bin/bash
# Test parallel production with a few mass points
# Run this BEFORE full production to verify everything works

set -e

echo "============================================"
echo "TEST: Parallel Production (3 mass points)"
echo "============================================"
echo ""

# Configuration
NEVENTS=10000  # Reduced for testing
MAX_PARALLEL=3

# Pythia library path
PYTHIA_ROOT="$(cd "$(dirname "$0")/../pythia/pythia8315" && pwd)"
export DYLD_LIBRARY_PATH="$PYTHIA_ROOT/lib:${DYLD_LIBRARY_PATH:-}"
export LD_LIBRARY_PATH="$PYTHIA_ROOT/lib:${LD_LIBRARY_PATH:-}"
OUTPUT_DIR="../../output/csv/simulation_test"
LOG_DIR="../../output/logs/simulation_test"

# Create directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_DIR"

# Test masses (one from each regime)
TEST_MASSES=(
    "0.5:muon:kaon"      # Kaon regime
    "2.0:muon:charm"     # D-meson regime
    "4.5:muon:beauty"    # B-meson regime
)

echo "Test configuration:"
echo "  Events per mass: $NEVENTS (reduced for testing)"
echo "  Parallel jobs: $MAX_PARALLEL"
echo "  Test masses: ${#TEST_MASSES[@]}"
echo ""
echo "Test points:"
for test in "${TEST_MASSES[@]}"; do
    IFS=':' read -r mass flavour regime <<< "$test"
    echo "  - $mass GeV $flavour ($regime)"
done
echo ""

read -p "Run test? (y/n) " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Test cancelled."
    exit 0
fi

# Job control
count_jobs() {
    jobs -r | wc -l | tr -d ' '
}

wait_for_slot() {
    while [ $(count_jobs) -ge $MAX_PARALLEL ]; do
        sleep 1
    done
}

# Run jobs
echo ""
echo "Starting test jobs..."
echo ""

start_time=$(date +%s)

for test in "${TEST_MASSES[@]}"; do
    IFS=':' read -r mass flavour regime <<< "$test"

    wait_for_slot

    {
        echo "[$(date +%H:%M:%S)] START: $mass GeV $flavour"
        ./main_hnl_production ${mass} ${flavour} $NEVENTS 2>&1
        exit_code=$?

        if [ $exit_code -eq 0 ]; then
            # Wait for filesystem sync and find the generated file
            sleep 0.5
            mass_label=$(echo "$mass" | sed 's/\.0*$//' | sed 's/\./p/')
            output_file=$(find . -maxdepth 1 -name "HNL_${mass_label}GeV_${flavour}_*.csv" -type f -mmin -2 2>/dev/null | head -1)

            if [ -n "$output_file" ] && [ -f "$output_file" ]; then
                mv "$output_file" "$OUTPUT_DIR/"
                echo "[$(date +%H:%M:%S)] SUCCESS: $(basename $output_file) → simulation_test/"
            else
                echo "[$(date +%H:%M:%S)] WARNING: No CSV found"
            fi
        else
            echo "[$(date +%H:%M:%S)] FAILED: exit code $exit_code"
        fi
    } > "${LOG_DIR}/test_${mass}GeV_${flavour}.log" 2>&1 &

    echo "Queued: $mass GeV $flavour"
done

echo ""
echo "Waiting for jobs to complete..."
wait

end_time=$(date +%s)
elapsed=$((end_time - start_time))

echo ""
echo "============================================"
echo "TEST COMPLETE"
echo "============================================"
echo "Time: ${elapsed}s"
echo ""

# Check results
csv_count=$(find "$OUTPUT_DIR" -name "*.csv" -type f 2>/dev/null | wc -l | tr -d ' ')
failed_count=$(grep -l "FAILED" "${LOG_DIR}"/test_*.log 2>/dev/null | wc -l | tr -d ' ')

echo "Results:"
echo "  CSV files created: $csv_count / ${#TEST_MASSES[@]}"
echo "  Failed jobs: $failed_count"
echo ""

if [ $csv_count -eq ${#TEST_MASSES[@]} ] && [ $failed_count -eq 0 ]; then
    echo "✓ TEST PASSED - Ready for full production"
    echo ""
    echo "Run full production with:"
    echo "  ./run_parallel_production.sh"
else
    echo "✗ TEST FAILED - Check logs:"
    echo "  ls -lh $LOG_DIR/test_*.log"
    exit 1
fi
