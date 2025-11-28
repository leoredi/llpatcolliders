#!/bin/bash
# Parallel HNL Production - Enhanced Mass Grid
# Runs 10 jobs in parallel to reduce runtime from 150h → ~1-2h

set -e

echo "============================================"
echo "HNL Parallel Production - Enhanced Mass Grid"
echo "============================================"
echo ""

# Configuration
NEVENTS=200000
OUTPUT_DIR="../output/csv/simulation_new"
LOG_DIR="../output/logs/simulation_new"
MAX_PARALLEL=10  # 10 parallel jobs

# Create directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_DIR"

# Timestamp
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
SUMMARY_LOG="${LOG_DIR}/parallel_run_${TIMESTAMP}.log"

echo "Configuration:" | tee "$SUMMARY_LOG"
echo "  Events per mass: $NEVENTS" | tee -a "$SUMMARY_LOG"
echo "  Parallel jobs: $MAX_PARALLEL" | tee -a "$SUMMARY_LOG"
echo "  Output dir: $OUTPUT_DIR" | tee -a "$SUMMARY_LOG"
echo "  Log dir: $LOG_DIR" | tee -a "$SUMMARY_LOG"
echo "" | tee -a "$SUMMARY_LOG"

# Enhanced mass grids
ELECTRON_MASSES=(
    0.2 0.22 0.25 0.28 0.3 0.32 0.35 0.38 0.40 0.42 0.45 0.48 0.50 0.52 0.55
    0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.75 1.8 1.82 1.85 1.9
    2.0 2.3 2.6 3.0 3.4 3.8 4.2 4.6 4.8 5.0 5.2
    6.0 7.0 8.0 9.0 10.0 12.0 15.0 20.0 30.0 40.0 50.0 60.0 80.0
)

MUON_MASSES=(
    0.2 0.22 0.25 0.28 0.3 0.32 0.35 0.37 0.38 0.39 0.40 0.42 0.45 0.48 0.50 0.55
    0.6 0.7 0.8 0.9 1.0 1.2 1.4 1.6 1.65 1.70 1.75 1.76 1.78 1.80 1.85 1.90
    2.0 2.3 2.6 3.0 3.4 3.8 4.2 4.6 4.8 5.0 5.2
    6.0 8.0 10.0 12.0 15.0 20.0 30.0 40.0 50.0 60.0 80.0
)

TAU_MASSES=(
    0.5 0.55 0.6 0.65 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.45 1.5 1.55
    1.6 1.62 1.64 1.66 1.70 1.74 1.78 1.80 1.85 1.9
    2.0 2.4 2.8 3.0 3.2 3.6 4.0 4.5
    6.0 7.0 8.0 10.0 15.0 20.0 30.0 40.0 50.0 60.0 80.0
)

# Function to wait for available slot
wait_for_slot() {
    while [ $(jobs -r | wc -l) -ge $MAX_PARALLEL ]; do
        sleep 2
    done
}

# Function to run a single mass point
run_mass_point() {
    local mass=$1
    local flavor=$2
    local mode=${3:-""}

    # Create temp log
    local temp_log="${LOG_DIR}/temp_${mass}_${flavor}_${mode}_$$.log"

    # Run simulation
    if [ -z "$mode" ]; then
        ./main_hnl_production $mass $flavor $NEVENTS > "$temp_log" 2>&1
    else
        ./main_hnl_production $mass $flavor $NEVENTS $mode > "$temp_log" 2>&1
    fi

    # Find output CSV
    local pattern="HNL_*${flavor}*"
    if [ -n "$mode" ]; then
        pattern="HNL_*${flavor}_*_${mode}.csv"
    else
        pattern="HNL_*${flavor}_*.csv"
    fi

    local output_file=$(ls -t $pattern 2>/dev/null | head -1)

    if [ -n "$output_file" ] && [ -f "$output_file" ]; then
        # Success - rename log to match CSV
        local log_name="${output_file%.csv}.log"
        mv "$temp_log" "${LOG_DIR}/${log_name}"
        mv "$output_file" "$OUTPUT_DIR/"
        echo "[SUCCESS] $mass GeV $flavor $mode → $output_file" | tee -a "$SUMMARY_LOG"
    else
        # Failed - keep temp log with FAILED prefix
        local failed_log="FAILED_${mass}GeV_${flavor}"
        [ -n "$mode" ] && failed_log="${failed_log}_${mode}"
        mv "$temp_log" "${LOG_DIR}/${failed_log}.log"
        echo "[FAILED]  $mass GeV $flavor $mode → check ${failed_log}.log" | tee -a "$SUMMARY_LOG"
    fi
}

echo "Starting parallel production..." | tee -a "$SUMMARY_LOG"
echo "" | tee -a "$SUMMARY_LOG"
start_time=$(date +%s)

# ============================================
# Electron Production
# ============================================
echo "Electron coupling (${#ELECTRON_MASSES[@]} points)..." | tee -a "$SUMMARY_LOG"
for mass in "${ELECTRON_MASSES[@]}"; do
    wait_for_slot
    run_mass_point $mass electron &
done

# ============================================
# Muon Production
# ============================================
echo "Muon coupling (${#MUON_MASSES[@]} points)..." | tee -a "$SUMMARY_LOG"
for mass in "${MUON_MASSES[@]}"; do
    wait_for_slot
    run_mass_point $mass muon &
done

# ============================================
# Tau Production (Dual Mode)
# ============================================
echo "Tau coupling (${#TAU_MASSES[@]} points, dual-mode where applicable)..." | tee -a "$SUMMARY_LOG"
for mass in "${TAU_MASSES[@]}"; do
    # Always run direct mode
    wait_for_slot
    run_mass_point $mass tau direct &

    # Run fromTau mode if kinematically allowed (m < 1.64 GeV)
    if (( $(echo "${mass} < 1.64" | bc -l) )); then
        wait_for_slot
        run_mass_point $mass tau fromTau &
    fi
done

# Wait for all background jobs to complete
echo "" | tee -a "$SUMMARY_LOG"
echo "Waiting for all jobs to complete..." | tee -a "$SUMMARY_LOG"
wait

# ============================================
# Summary
# ============================================
end_time=$(date +%s)
duration=$((end_time - start_time))
hours=$((duration / 3600))
minutes=$(((duration % 3600) / 60))

echo "" | tee -a "$SUMMARY_LOG"
echo "============================================" | tee -a "$SUMMARY_LOG"
echo "PRODUCTION COMPLETE!" | tee -a "$SUMMARY_LOG"
echo "============================================" | tee -a "$SUMMARY_LOG"
echo "" | tee -a "$SUMMARY_LOG"
echo "Runtime: ${hours}h ${minutes}m" | tee -a "$SUMMARY_LOG"
echo "" | tee -a "$SUMMARY_LOG"

# Count results
total_csv=$(ls "$OUTPUT_DIR"/*.csv 2>/dev/null | wc -l)
total_failed=$(ls "$LOG_DIR"/FAILED_*.log 2>/dev/null | wc -l)

echo "Results:" | tee -a "$SUMMARY_LOG"
echo "  Successful: $total_csv CSV files" | tee -a "$SUMMARY_LOG"
echo "  Failed: $total_failed runs" | tee -a "$SUMMARY_LOG"
echo "" | tee -a "$SUMMARY_LOG"

if [ $total_failed -gt 0 ]; then
    echo "Failed runs:" | tee -a "$SUMMARY_LOG"
    ls -1 "$LOG_DIR"/FAILED_*.log | sed 's|.*/||' | tee -a "$SUMMARY_LOG"
    echo "" | tee -a "$SUMMARY_LOG"
fi

echo "Output directory: $OUTPUT_DIR" | tee -a "$SUMMARY_LOG"
echo "Total disk usage: $(du -sh $OUTPUT_DIR | cut -f1)" | tee -a "$SUMMARY_LOG"
echo "" | tee -a "$SUMMARY_LOG"
echo "Summary log: $SUMMARY_LOG" | tee -a "$SUMMARY_LOG"
echo "" | tee -a "$SUMMARY_LOG"
