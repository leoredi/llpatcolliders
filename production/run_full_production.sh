#!/bin/bash
# Full HNL Production Run - Uniform 0.1 GeV Grid
# Mass range: 0.2 - 10.0 GeV with 0.1 GeV spacing
# Standard 200k events per point

set -e  # Exit on error

echo "============================================"
echo "HNL Full Production - Uniform Mass Grid"
echo "============================================"
echo ""
echo "Mass grid: 0.1 GeV spacing from 0.2 to 10.0 GeV"
echo "Total mass points: 99 per flavor"
echo "Total events: 200k per mass point"
echo ""

# Configuration
NEVENTS=200000
OUTPUT_DIR="../output/csv/simulation_new"
LOG_DIR="../output/logs/simulation_new"

# Create directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_DIR"

# Timestamp for this run
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOGFILE="${LOG_DIR}/production_run_${TIMESTAMP}.log"

echo "Output directory: $OUTPUT_DIR"
echo "Log file: $LOGFILE"
echo ""

# ===========================================================================
# Uniform Mass Grid - Electron Coupling (BC6)
# ===========================================================================
# 0.1 GeV spacing from 0.2 to 10.0 GeV: 99 points

ELECTRON_MASSES=($(seq 0.2 0.1 10.0))

# ===========================================================================
# Uniform Mass Grid - Muon Coupling (BC7)
# ===========================================================================
# 0.1 GeV spacing from 0.2 to 10.0 GeV: 99 points

MUON_MASSES=($(seq 0.2 0.1 10.0))

# ===========================================================================
# Uniform Mass Grid - Tau Coupling (BC8)
# ===========================================================================
# 0.1 GeV spacing from 0.2 to 10.0 GeV: 99 points
#
# DUAL MODE: direct + fromTau where kinematically allowed
# fromTau available when: m_HNL + m_pi < m_tau (m < 1.64 GeV)

TAU_MASSES=($(seq 0.2 0.1 10.0))

# ===========================================================================
# Execution
# ===========================================================================

start_time=$(date +%s)
total_jobs=$((${#ELECTRON_MASSES[@]} + ${#MUON_MASSES[@]} + ${#TAU_MASSES[@]} * 2))
current_job=0

echo "Total mass points:"
echo "  Electron: ${#ELECTRON_MASSES[@]}"
echo "  Muon: ${#MUON_MASSES[@]}"
echo "  Tau: ${#TAU_MASSES[@]} (some with dual mode for m < 1.64 GeV)"
echo "  TOTAL: ~312 simulation runs"
echo ""
echo "Estimated time: ~10 hours (single core), ~1 hour (10 cores)"
echo "Recommendation: Run in parallel on multi-core system"
echo ""
read -p "Continue with full production? (y/n) " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Production cancelled."
    exit 0
fi
echo ""

# ===========================================================================
# Electron Production
# ===========================================================================

echo "============================================" | tee -a "$LOGFILE"
echo "ELECTRON COUPLING (BC6)" | tee -a "$LOGFILE"
echo "============================================" | tee -a "$LOGFILE"
echo "" | tee -a "$LOGFILE"

for mass in "${ELECTRON_MASSES[@]}"; do
    current_job=$((current_job + 1))
    echo "[$current_job/$total_jobs] Running ${mass} GeV electron..." | tee -a "$LOGFILE"

    # Run and capture output
    ./main_hnl_production ${mass} electron $NEVENTS > "${LOG_DIR}/temp_run.log" 2>&1
    tail -5 "${LOG_DIR}/temp_run.log" | tee -a "$LOGFILE"

    # Move output and create matching log
    output_file=$(ls -t HNL_*GeV_electron_*.csv 2>/dev/null | head -1)
    if [ -n "$output_file" ]; then
        # Derive log name from CSV name (replace .csv with .log)
        log_name="${output_file%.csv}.log"
        mv "${LOG_DIR}/temp_run.log" "${LOG_DIR}/${log_name}"
        mv "$output_file" "$OUTPUT_DIR/"
        echo "  → CSV: $OUTPUT_DIR/$output_file" | tee -a "$LOGFILE"
        echo "  → LOG: $LOG_DIR/$log_name" | tee -a "$LOGFILE"
    else
        # No CSV created - keep temp log for debugging
        mv "${LOG_DIR}/temp_run.log" "${LOG_DIR}/FAILED_${mass}GeV_electron.log"
        echo "  → FAILED - check $LOG_DIR/FAILED_${mass}GeV_electron.log" | tee -a "$LOGFILE"
    fi
    echo "" | tee -a "$LOGFILE"
done

# ===========================================================================
# Muon Production
# ===========================================================================

echo "============================================" | tee -a "$LOGFILE"
echo "MUON COUPLING (BC7)" | tee -a "$LOGFILE"
echo "============================================" | tee -a "$LOGFILE"
echo "" | tee -a "$LOGFILE"

for mass in "${MUON_MASSES[@]}"; do
    current_job=$((current_job + 1))
    echo "[$current_job/$total_jobs] Running ${mass} GeV muon..." | tee -a "$LOGFILE"

    # Run and capture output
    ./main_hnl_production ${mass} muon $NEVENTS > "${LOG_DIR}/temp_run.log" 2>&1
    tail -5 "${LOG_DIR}/temp_run.log" | tee -a "$LOGFILE"

    # Move output and create matching log
    output_file=$(ls -t HNL_*GeV_muon_*.csv 2>/dev/null | head -1)
    if [ -n "$output_file" ]; then
        # Derive log name from CSV name (replace .csv with .log)
        log_name="${output_file%.csv}.log"
        mv "${LOG_DIR}/temp_run.log" "${LOG_DIR}/${log_name}"
        mv "$output_file" "$OUTPUT_DIR/"
        echo "  → CSV: $OUTPUT_DIR/$output_file" | tee -a "$LOGFILE"
        echo "  → LOG: $LOG_DIR/$log_name" | tee -a "$LOGFILE"
    else
        # No CSV created - keep temp log for debugging
        mv "${LOG_DIR}/temp_run.log" "${LOG_DIR}/FAILED_${mass}GeV_muon.log"
        echo "  → FAILED - check $LOG_DIR/FAILED_${mass}GeV_muon.log" | tee -a "$LOGFILE"
    fi
    echo "" | tee -a "$LOGFILE"
done

# ===========================================================================
# Tau Production (DUAL MODE)
# ===========================================================================

echo "============================================" | tee -a "$LOGFILE"
echo "TAU COUPLING (BC8) - DUAL MODE" | tee -a "$LOGFILE"
echo "============================================" | tee -a "$LOGFILE"
echo "" | tee -a "$LOGFILE"

for mass in "${TAU_MASSES[@]}"; do
    # Always run direct mode
    current_job=$((current_job + 1))
    echo "[$current_job/$total_jobs] Running ${mass} GeV tau (direct mode)..." | tee -a "$LOGFILE"

    # Run and capture output
    ./main_hnl_production ${mass} tau $NEVENTS direct > "${LOG_DIR}/temp_run.log" 2>&1
    tail -5 "${LOG_DIR}/temp_run.log" | tee -a "$LOGFILE"

    # Move output and create matching log
    output_file=$(ls -t HNL_*GeV_tau_*_direct.csv 2>/dev/null | head -1)
    if [ -n "$output_file" ]; then
        # Derive log name from CSV name (replace .csv with .log)
        log_name="${output_file%.csv}.log"
        mv "${LOG_DIR}/temp_run.log" "${LOG_DIR}/${log_name}"
        mv "$output_file" "$OUTPUT_DIR/"
        echo "  → CSV: $OUTPUT_DIR/$output_file" | tee -a "$LOGFILE"
        echo "  → LOG: $LOG_DIR/$log_name" | tee -a "$LOGFILE"
    else
        # No CSV created - keep temp log for debugging
        mv "${LOG_DIR}/temp_run.log" "${LOG_DIR}/FAILED_${mass}GeV_tau_direct.log"
        echo "  → FAILED - check $LOG_DIR/FAILED_${mass}GeV_tau_direct.log" | tee -a "$LOGFILE"
    fi
    echo "" | tee -a "$LOGFILE"

    # Run fromTau mode if kinematically allowed
    # Condition: m_HNL + m_pi < m_tau
    # m_pi = 0.14 GeV, m_tau = 1.777 GeV
    # Therefore: m_HNL < 1.64 GeV
    if (( $(echo "${mass} < 1.64" | bc -l) )); then
        current_job=$((current_job + 1))
        echo "[$current_job/$total_jobs] Running ${mass} GeV tau (fromTau mode)..." | tee -a "$LOGFILE"

        # Run and capture output
        ./main_hnl_production ${mass} tau $NEVENTS fromTau > "${LOG_DIR}/temp_run.log" 2>&1
        tail -5 "${LOG_DIR}/temp_run.log" | tee -a "$LOGFILE"

        # Move output and create matching log
        output_file=$(ls -t HNL_*GeV_tau_*_fromTau.csv 2>/dev/null | head -1)
        if [ -n "$output_file" ]; then
            # Derive log name from CSV name (replace .csv with .log)
            log_name="${output_file%.csv}.log"
            mv "${LOG_DIR}/temp_run.log" "${LOG_DIR}/${log_name}"
            mv "$output_file" "$OUTPUT_DIR/"
            echo "  → CSV: $OUTPUT_DIR/$output_file" | tee -a "$LOGFILE"
            echo "  → LOG: $LOG_DIR/$log_name" | tee -a "$LOGFILE"
        else
            # No CSV created - keep temp log for debugging
            mv "${LOG_DIR}/temp_run.log" "${LOG_DIR}/FAILED_${mass}GeV_tau_fromTau.log"
            echo "  → FAILED - check $LOG_DIR/FAILED_${mass}GeV_tau_fromTau.log" | tee -a "$LOGFILE"
        fi
        echo "" | tee -a "$LOGFILE"
    else
        echo "  (fromTau mode not kinematically allowed at this mass)" | tee -a "$LOGFILE"
        echo "" | tee -a "$LOGFILE"
    fi
done

# ===========================================================================
# Summary
# ===========================================================================

end_time=$(date +%s)
duration=$((end_time - start_time))
hours=$((duration / 3600))
minutes=$(((duration % 3600) / 60))

echo "============================================" | tee -a "$LOGFILE"
echo "PRODUCTION COMPLETE!" | tee -a "$LOGFILE"
echo "============================================" | tee -a "$LOGFILE"
echo "" | tee -a "$LOGFILE"
echo "Total runtime: ${hours}h ${minutes}m" | tee -a "$LOGFILE"
echo "" | tee -a "$LOGFILE"
echo "Output files in: $OUTPUT_DIR" | tee -a "$LOGFILE"
ls -lh "$OUTPUT_DIR" | tee -a "$LOGFILE"
echo "" | tee -a "$LOGFILE"
echo "Total disk usage:" | tee -a "$LOGFILE"
du -sh "$OUTPUT_DIR" | tee -a "$LOGFILE"
echo "" | tee -a "$LOGFILE"
echo "Log file: $LOGFILE" | tee -a "$LOGFILE"
echo "" | tee -a "$LOGFILE"
echo "Next steps:" | tee -a "$LOGFILE"
echo "  1. Run geometry preprocessing" | tee -a "$LOGFILE"
echo "  2. Calculate limits" | tee -a "$LOGFILE"
echo "  3. Generate exclusion plots" | tee -a "$LOGFILE"
echo "" | tee -a "$LOGFILE"
