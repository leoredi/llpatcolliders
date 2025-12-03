#!/bin/bash
# Parallel HNL Production Run
# Generates all meson (Pythia) mass points using parallel job execution

set -e  # Exit on error

echo "============================================"
echo "HNL Production - Parallel Execution"
echo "============================================"
echo ""

# Configuration
NEVENTS=100000
MAX_PARALLEL=8  # Adjust based on your CPU cores (leave 1-2 for system)

# Pythia library path (needed on macOS for libpythia8.dylib lookup)
PYTHIA_ROOT="$(cd "$(dirname "$0")/../pythia/pythia8315" && pwd)"
export DYLD_LIBRARY_PATH="$PYTHIA_ROOT/lib:${DYLD_LIBRARY_PATH:-}"
export LD_LIBRARY_PATH="$PYTHIA_ROOT/lib:${LD_LIBRARY_PATH:-}"
OUTPUT_DIR="../../output/csv/simulation_new"
LOG_DIR="../../output/logs/simulation_new"

# Create directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_DIR"

# Timestamp for this run
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOGFILE="${LOG_DIR}/production_run_${TIMESTAMP}.log"

echo "Configuration:" | tee "$LOGFILE"
echo "  Events per mass: $NEVENTS" | tee -a "$LOGFILE"
echo "  Parallel jobs: $MAX_PARALLEL" | tee -a "$LOGFILE"
echo "  Output directory: $OUTPUT_DIR" | tee -a "$LOGFILE"
echo "  Log file: $LOGFILE" | tee -a "$LOGFILE"
echo ""

# ===========================================================================
# Load Mass Grids from Central Configuration
# ===========================================================================
source ./load_mass_grid.sh

# Use meson (Pythia) mass arrays
ELECTRON_MASSES=("${ELECTRON_MASSES_MESON[@]}")
MUON_MASSES=("${MUON_MASSES_MESON[@]}")
TAU_MASSES=("${TAU_MASSES_MESON[@]}")

# ===========================================================================
# Job Control Functions
# ===========================================================================

# Count running background jobs
count_jobs() {
    jobs -r | wc -l | tr -d ' '
}

# Wait for a job slot to become available
wait_for_slot() {
    while [ $(count_jobs) -ge $MAX_PARALLEL ]; do
        sleep 1
    done
}

# Run a single production job in background
run_production_job() {
    local mass=$1
    local flavour=$2
    local mode=${3:-direct}  # direct or fromTau

    local log_file="${LOG_DIR}/HNL_${mass}GeV_${flavour}_${mode}_${TIMESTAMP}.log"

    {
        echo "[$(date +%H:%M:%S)] Starting: $mass GeV $flavour ($mode)"

        if [ "$mode" = "fromTau" ]; then
            ./main_hnl_production ${mass} ${flavour} $NEVENTS fromTau 2>&1
        else
            ./main_hnl_production ${mass} ${flavour} $NEVENTS 2>&1
        fi

        local exit_code=$?

        if [ $exit_code -eq 0 ]; then
            # Find the generated CSV - wait a moment for filesystem sync
            sleep 0.5

            # Convert mass to filename format (e.g., 2.6 → 2p6)
            local mass_label=$(echo "$mass" | sed 's/\.0*$//' | sed 's/\./p/')

            # Search for CSV file matching this job
            local output_file=$(find . -maxdepth 1 -name "HNL_${mass_label}GeV_${flavour}_*.csv" -type f -mmin -2 2>/dev/null | head -1)

            if [ -n "$output_file" ] && [ -f "$output_file" ]; then
                mv "$output_file" "$OUTPUT_DIR/" 2>/dev/null
                local basename=$(basename "$output_file")
                echo "[$(date +%H:%M:%S)] SUCCESS: $basename → simulation_new/"
            else
                echo "[$(date +%H:%M:%S)] WARNING: No CSV found for $mass GeV $flavour (searched: HNL_${mass_label}GeV_${flavour}_*.csv)"
            fi
        else
            echo "[$(date +%H:%M:%S)] FAILED: $mass GeV $flavour (exit code $exit_code)"
        fi
    } > "$log_file" 2>&1
}

# ===========================================================================
# Calculate Total Jobs
# ===========================================================================

total_electron=${#ELECTRON_MASSES[@]}
total_muon=${#MUON_MASSES[@]}
total_tau=0

# Count tau jobs (some masses have dual mode)
for mass in "${TAU_MASSES[@]}"; do
    total_tau=$((total_tau + 1))
    # Check if this mass needs fromTau mode (m < 1.64 GeV)
    if (( $(echo "$mass < 1.64" | bc -l) )); then
        total_tau=$((total_tau + 1))
    fi
done

total_jobs=$((total_electron + total_muon + total_tau))

echo "Total mass points:" | tee -a "$LOGFILE"
echo "  Electron: $total_electron points" | tee -a "$LOGFILE"
echo "  Muon: $total_muon points" | tee -a "$LOGFILE"
echo "  Tau: $total_tau runs (including dual-mode)" | tee -a "$LOGFILE"
echo "  TOTAL: $total_jobs simulation runs" | tee -a "$LOGFILE"
echo "" | tee -a "$LOGFILE"

# Estimate time
single_core_hours=$((total_jobs / 6))  # ~6 jobs per hour on average
parallel_hours=$((single_core_hours / MAX_PARALLEL + 1))

echo "Estimated time:" | tee -a "$LOGFILE"
echo "  Single core: ~$single_core_hours hours" | tee -a "$LOGFILE"
echo "  $MAX_PARALLEL cores: ~$parallel_hours hours" | tee -a "$LOGFILE"
echo "" | tee -a "$LOGFILE"

read -p "Continue with parallel production? (y/n) " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Production cancelled." | tee -a "$LOGFILE"
    exit 0
fi
echo "" | tee -a "$LOGFILE"

start_time=$(date +%s)
completed_jobs=0

# ===========================================================================
# Electron Production
# ===========================================================================

echo "============================================" | tee -a "$LOGFILE"
echo "ELECTRON COUPLING (BC6)" | tee -a "$LOGFILE"
echo "============================================" | tee -a "$LOGFILE"

for mass in "${ELECTRON_MASSES[@]}"; do
    wait_for_slot
    run_production_job "$mass" "electron" "direct" &
    completed_jobs=$((completed_jobs + 1))
    echo "[$completed_jobs/$total_jobs] Queued: $mass GeV electron" | tee -a "$LOGFILE"
done

# ===========================================================================
# Muon Production
# ===========================================================================

echo "" | tee -a "$LOGFILE"
echo "============================================" | tee -a "$LOGFILE"
echo "MUON COUPLING (BC7)" | tee -a "$LOGFILE"
echo "============================================" | tee -a "$LOGFILE"

for mass in "${MUON_MASSES[@]}"; do
    wait_for_slot
    run_production_job "$mass" "muon" "direct" &
    completed_jobs=$((completed_jobs + 1))
    echo "[$completed_jobs/$total_jobs] Queued: $mass GeV muon" | tee -a "$LOGFILE"
done

# ===========================================================================
# Tau Production (DUAL MODE)
# ===========================================================================

echo "" | tee -a "$LOGFILE"
echo "============================================" | tee -a "$LOGFILE"
echo "TAU COUPLING (BC8) - DUAL MODE" | tee -a "$LOGFILE"
echo "============================================" | tee -a "$LOGFILE"

for mass in "${TAU_MASSES[@]}"; do
    # MODE A: Direct production (all masses)
    wait_for_slot
    run_production_job "$mass" "tau" "direct" &
    completed_jobs=$((completed_jobs + 1))
    echo "[$completed_jobs/$total_jobs] Queued: $mass GeV tau (direct)" | tee -a "$LOGFILE"

    # MODE B: fromTau cascade (only m < 1.64 GeV)
    if (( $(echo "$mass < 1.64" | bc -l) )); then
        wait_for_slot
        run_production_job "$mass" "tau" "fromTau" &
        completed_jobs=$((completed_jobs + 1))
        echo "[$completed_jobs/$total_jobs] Queued: $mass GeV tau (fromTau)" | tee -a "$LOGFILE"
    fi
done

# ===========================================================================
# Wait for All Jobs to Complete
# ===========================================================================

echo "" | tee -a "$LOGFILE"
echo "============================================" | tee -a "$LOGFILE"
echo "Waiting for all jobs to complete..." | tee -a "$LOGFILE"
echo "============================================" | tee -a "$LOGFILE"

wait

end_time=$(date +%s)
elapsed=$((end_time - start_time))
elapsed_hours=$((elapsed / 3600))
elapsed_mins=$(((elapsed % 3600) / 60))

echo "" | tee -a "$LOGFILE"
echo "============================================" | tee -a "$LOGFILE"
echo "PRODUCTION COMPLETE" | tee -a "$LOGFILE"
echo "============================================" | tee -a "$LOGFILE"
echo "Total time: ${elapsed_hours}h ${elapsed_mins}m" | tee -a "$LOGFILE"
echo "" | tee -a "$LOGFILE"

# ===========================================================================
# Summary Statistics
# ===========================================================================

total_csv=$(find "$OUTPUT_DIR" -name "HNL_*.csv" -type f | wc -l | tr -d ' ')
total_size=$(du -sh "$OUTPUT_DIR" | cut -f1)

echo "Output summary:" | tee -a "$LOGFILE"
echo "  CSV files: $total_csv" | tee -a "$LOGFILE"
echo "  Total size: $total_size" | tee -a "$LOGFILE"
echo "  Location: $OUTPUT_DIR" | tee -a "$LOGFILE"
echo "" | tee -a "$LOGFILE"

# Check for failures
failed_logs=$(find "$LOG_DIR" -name "*${TIMESTAMP}*.log" -type f -exec grep -l "FAILED\|ERROR" {} \; | wc -l | tr -d ' ')

if [ $failed_logs -gt 0 ]; then
    echo "⚠ Warning: $failed_logs job(s) may have failed" | tee -a "$LOGFILE"
    echo "  Check logs in: $LOG_DIR" | tee -a "$LOGFILE"
else
    echo "✓ All jobs completed successfully" | tee -a "$LOGFILE"
fi

echo "" | tee -a "$LOGFILE"
echo "Next steps:" | tee -a "$LOGFILE"
echo "1. Verify CSV files: ls -lh $OUTPUT_DIR/*.csv" | tee -a "$LOGFILE"
echo "2. Check for overlaps: python analysis_pbc/limits/combine_production_channels.py --dry-run" | tee -a "$LOGFILE"
echo "3. Run analysis: cd analysis_pbc && python limits/run_serial.py" | tee -a "$LOGFILE"
echo "============================================" | tee -a "$LOGFILE"
