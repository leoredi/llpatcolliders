#!/bin/bash
# Full HNL Production Run - Enhanced Mass Grid
# Increased sampling near kinematic thresholds
# Standard 200k events per point (no increase in simulation size)

set -e  # Exit on error

echo "============================================"
echo "HNL Full Production - Enhanced Mass Grid"
echo "============================================"
echo ""
echo "Enhanced sampling near thresholds:"
echo "  - K threshold: m_K - m_lepton"
echo "  - D threshold: m_D - m_lepton"
echo "  - B threshold: m_B - m_lepton"
echo "  - Tau threshold: m_tau = 1.777 GeV"
echo "  - W threshold: m_W = 80 GeV"
echo ""
echo "Total events: 200k per mass point (unchanged)"
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
# Enhanced Mass Grid - Electron Coupling (BC6)
# ===========================================================================
# Standard: 38 points
# Enhanced: 50 points (added threshold regions)

ELECTRON_MASSES=(
    # Ultra-low mass (kaon regime, below K threshold)
    0.2 0.22 0.25 0.28 0.3 0.32 0.35
    # K threshold region (m_K - m_e ≈ 0.493 GeV)
    0.38 0.40 0.42 0.45 0.48 0.50 0.52
    # Charm regime
    0.55 0.6 0.7 0.8 0.9 1.0 1.1 1.2
    # D threshold region (m_D - m_e ≈ 1.87 GeV)
    1.3 1.4 1.5 1.6 1.7 1.75 1.8 1.82 1.85 1.9
    # Intermediate
    2.0 2.3 2.6 3.0 3.4 3.8 4.2 4.6
    # B threshold approach (m_B - m_e ≈ 5.28 GeV)
    4.8 5.0 5.2
    # High mass (EW regime) - Enhanced granularity up to 40 GeV
    6.0 7.0 8.0 9.0 10.0 11.0 12.0 13.0 14.0 15.0 16.0 17.0 18.0 19.0 20.0 22.0 25.0 28.0 30.0 32.0 35.0 38.0 40.0
)

# ===========================================================================
# Enhanced Mass Grid - Muon Coupling (BC7)
# ===========================================================================
# Standard: 38 points
# Enhanced: 48 points (added threshold regions)

MUON_MASSES=(
    # Low mass (kaon regime, below K threshold)
    0.2 0.22 0.25 0.28 0.3 0.32 0.35
    # K threshold region (m_K - m_mu ≈ 0.388 GeV)
    0.37 0.38 0.39 0.40 0.42 0.45 0.48 0.50
    # Charm regime
    0.55 0.6 0.7 0.8 0.9 1.0 1.2 1.4 1.6
    # D threshold region (m_D - m_mu ≈ 1.76 GeV)
    1.65 1.70 1.75 1.76 1.78 1.8 1.85 1.9
    # Intermediate
    2.0 2.3 2.6 3.0 3.4 3.8 4.2 4.6
    # B threshold approach (m_B - m_mu ≈ 5.17 GeV)
    4.8 5.0 5.2
    # High mass (EW regime) - Enhanced granularity up to 40 GeV
    6.0 7.0 8.0 9.0 10.0 11.0 12.0 13.0 14.0 15.0 16.0 17.0 18.0 19.0 20.0 22.0 25.0 28.0 30.0 32.0 35.0 38.0 40.0
)

# ===========================================================================
# Enhanced Mass Grid - Tau Coupling (BC8)
# ===========================================================================
# Standard: 26 points (starts at 0.5 GeV)
# Enhanced: 35 points (added tau threshold region)
#
# DUAL MODE: direct + fromTau where kinematically allowed
# fromTau available when: m_HNL + m_pi < m_tau (m < 1.64 GeV)

TAU_MASSES=(
    # Charm regime
    0.5 0.55 0.6 0.65 0.7 0.8 0.9 1.0 1.1 1.2 1.3
    # Below tau mass (both modes available)
    1.4 1.45 1.5 1.55 1.6 1.62 1.64
    # Tau threshold region (m_tau = 1.777 GeV)
    1.66 1.70 1.74 1.777 1.8 1.85 1.9
    # Beauty regime
    2.0 2.4 2.8 3.2 3.6 4.0 4.5
    # High mass (EW regime)
    5.0 6.0 7.0 8.0 10.0 12.0 15.0 20.0 30.0 40.0 50.0 60.0 80.0
)

# ===========================================================================
# Execution
# ===========================================================================

start_time=$(date +%s)
total_jobs=$((${#ELECTRON_MASSES[@]} + ${#MUON_MASSES[@]} + ${#TAU_MASSES[@]} * 2))
current_job=0

echo "Total mass points:"
echo "  Electron: ${#ELECTRON_MASSES[@]}"
echo "  Muon: ${#MUON_MASSES[@]}"
echo "  Tau: ${#TAU_MASSES[@]} (x2 for dual mode where applicable)"
echo "  TOTAL: ~$total_jobs simulation runs"
echo ""
echo "Estimated time: 120-150 hours (single core)"
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
