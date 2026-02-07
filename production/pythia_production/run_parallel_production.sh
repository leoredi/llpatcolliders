#!/bin/bash

set -e

FLAVOUR="${1:-all}"
FLAVOUR=$(echo "$FLAVOUR" | tr '[:upper:]' '[:lower:]')
MODE="${2:-both}"
MODE=$(echo "$MODE" | tr '[:upper:]' '[:lower:]')

if [[ ! "$FLAVOUR" =~ ^(electron|muon|tau|all)$ ]]; then
    echo "Error: Invalid flavour '$FLAVOUR'"
    echo "Usage: $0 [electron|muon|tau|all] [direct|fromTau|both]"
    exit 1
fi

if [[ ! "$MODE" =~ ^(direct|fromtau|both)$ ]]; then
    echo "Error: Invalid mode '$MODE'"
    echo "Usage: $0 [electron|muon|tau|all] [direct|fromTau|both]"
    exit 1
fi

if [[ "$MODE" == "fromtau" && "$FLAVOUR" != "tau" && "$FLAVOUR" != "all" ]]; then
    echo "Error: 'fromTau' mode only valid for tau flavour"
    echo "  fromTau uses τ → N X decay chain, which requires tau coupling"
    exit 1
fi

echo "============================================"
echo "HNL Production - Parallel Execution"
echo "Flavour: $FLAVOUR"
echo "Mode: $MODE"
echo "============================================"
echo ""

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
NEVENTS=100000
NEVENTS_FROMTAU=$NEVENTS
MAX_PARALLEL=12

FROMTAU_MASS_THRESHOLD=1.77

PYTHIA_ROOT="$SCRIPT_DIR/pythia8315"
export DYLD_LIBRARY_PATH="$PYTHIA_ROOT/lib:${DYLD_LIBRARY_PATH:-}"
export LD_LIBRARY_PATH="$PYTHIA_ROOT/lib:${LD_LIBRARY_PATH:-}"
export PYTHIA8DATA="$PYTHIA_ROOT/share/Pythia8/xmldoc"
OUTPUT_DIR="$SCRIPT_DIR/../../output/csv/simulation"
LOG_DIR="$SCRIPT_DIR/../../output/logs/simulation"

mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_DIR"

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOGFILE="${LOG_DIR}/production_run_${TIMESTAMP}.log"

echo "Configuration:" | tee "$LOGFILE"
echo "  Events per mass: $NEVENTS" | tee -a "$LOGFILE"
echo "  Parallel jobs: $MAX_PARALLEL" | tee -a "$LOGFILE"
echo "  Output directory: $OUTPUT_DIR" | tee -a "$LOGFILE"
echo "  Log file: $LOGFILE" | tee -a "$LOGFILE"
echo ""

source ./load_mass_grid.sh

ELECTRON_MASSES=("${MASS_GRID[@]}")
MUON_MASSES=("${MASS_GRID[@]}")
TAU_MASSES=("${MASS_GRID[@]}")

if [[ "$FLAVOUR" != "all" ]]; then
    [[ "$FLAVOUR" != "electron" ]] && ELECTRON_MASSES=()
    [[ "$FLAVOUR" != "muon" ]] && MUON_MASSES=()
    [[ "$FLAVOUR" != "tau" ]] && TAU_MASSES=()
fi


count_jobs() {
    jobs -r | wc -l | tr -d ' '
}

wait_for_slot() {
    while [ $(count_jobs) -ge $MAX_PARALLEL ]; do
        sleep 1
    done
}

run_production_job() {
    local mass=$1
    local flavour=$2
    local mode=${3:-direct}

    local log_file="${LOG_DIR}/HNL_${mass}GeV_${flavour}_${mode}_${TIMESTAMP}.log"
    local script_dir=$(pwd)

    {
        echo "[$(date +%H:%M:%S)] Starting: $mass GeV $flavour ($mode)"

        if [ "$mode" = "fromTau" ]; then
            ./main_hnl_production ${mass} ${flavour} $NEVENTS_FROMTAU fromTau 2>&1
        else
            ./main_hnl_production ${mass} ${flavour} $NEVENTS 2>&1
        fi

        local exit_code=$?

        if [ $exit_code -eq 0 ]; then
            sleep 0.5
            local mass_label=$(printf "%.2f" "$mass" | tr '.' 'p')
            local csv_file=$(ls HNL_${mass_label}GeV_${flavour}_*.csv 2>/dev/null | head -1)

            if [ -n "$csv_file" ]; then
                mv "$csv_file" "$OUTPUT_DIR/"
                echo "[$(date +%H:%M:%S)] SUCCESS: $csv_file → simulation/"
            else
                echo "[$(date +%H:%M:%S)] WARNING: CSV not found"
            fi
        else
            echo "[$(date +%H:%M:%S)] FAILED: $mass GeV $flavour (exit $exit_code)"
        fi
    } > "$log_file" 2>&1
}


count_tau_runs() {
    local count=0
    for mass in "${TAU_MASSES[@]}"; do
        if [[ "$MODE" == "direct" || "$MODE" == "both" ]]; then
            count=$((count + 1))
        fi
        if [[ "$MODE" == "fromtau" || "$MODE" == "both" ]]; then
            if (( $(echo "$mass < $FROMTAU_MASS_THRESHOLD" | bc -l) )); then
                count=$((count + 1))
            fi
        fi
    done
    echo $count
}

total_electron=${#ELECTRON_MASSES[@]}
total_muon=${#MUON_MASSES[@]}
total_tau=$(count_tau_runs)
total_jobs=$((total_electron + total_muon + total_tau))

echo "Total mass points:" | tee -a "$LOGFILE"
echo "  Electron: $total_electron points" | tee -a "$LOGFILE"
echo "  Muon: $total_muon points" | tee -a "$LOGFILE"
echo "  Tau: $total_tau runs (including dual-mode)" | tee -a "$LOGFILE"
echo "  TOTAL: $total_jobs simulation runs" | tee -a "$LOGFILE"
echo "" | tee -a "$LOGFILE"

single_core_hours=$((total_jobs / 6))
parallel_hours=$((single_core_hours / MAX_PARALLEL + 1))

echo "Estimated time:" | tee -a "$LOGFILE"
echo "  Single core: ~$single_core_hours hours" | tee -a "$LOGFILE"
echo "  $MAX_PARALLEL cores: ~$parallel_hours hours" | tee -a "$LOGFILE"
echo "" | tee -a "$LOGFILE"

echo "Starting production..." | tee -a "$LOGFILE"
echo "" | tee -a "$LOGFILE"

start_time=$(date +%s)
completed_jobs=0


if [[ "$FLAVOUR" == "electron" || "$FLAVOUR" == "all" ]]; then
    echo "============================================" | tee -a "$LOGFILE"
    echo "ELECTRON COUPLING (BC6)" | tee -a "$LOGFILE"
    echo "============================================" | tee -a "$LOGFILE"

    for mass in "${ELECTRON_MASSES[@]}"; do
        wait_for_slot
        run_production_job "$mass" "electron" "direct" &
        completed_jobs=$((completed_jobs + 1))
        echo "[$completed_jobs/$total_jobs] Queued: $mass GeV electron" | tee -a "$LOGFILE"
    done
fi


if [[ "$FLAVOUR" == "muon" || "$FLAVOUR" == "all" ]]; then
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
fi


if [[ "$FLAVOUR" == "tau" || "$FLAVOUR" == "all" ]]; then
    echo "" | tee -a "$LOGFILE"
    echo "============================================" | tee -a "$LOGFILE"
    echo "TAU COUPLING (BC8) - MODE: $MODE" | tee -a "$LOGFILE"
    echo "============================================" | tee -a "$LOGFILE"

    for mass in "${TAU_MASSES[@]}"; do
        if [[ "$MODE" == "direct" || "$MODE" == "both" ]]; then
            wait_for_slot
            run_production_job "$mass" "tau" "direct" &
            completed_jobs=$((completed_jobs + 1))
            echo "[$completed_jobs/$total_jobs] Queued: $mass GeV tau (direct)" | tee -a "$LOGFILE"
        fi

        if [[ "$MODE" == "fromtau" || "$MODE" == "both" ]]; then
            if (( $(echo "$mass < $FROMTAU_MASS_THRESHOLD" | bc -l) )); then
                wait_for_slot
                run_production_job "$mass" "tau" "fromTau" &
                completed_jobs=$((completed_jobs + 1))
                echo "[$completed_jobs/$total_jobs] Queued: $mass GeV tau (fromTau)" | tee -a "$LOGFILE"
            fi
        fi
    done
fi


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


total_csv=$(find "$OUTPUT_DIR" -name "HNL_*.csv" -type f | wc -l | tr -d ' ')
total_size=$(du -sh "$OUTPUT_DIR" | cut -f1)

echo "Output summary:" | tee -a "$LOGFILE"
echo "  CSV files: $total_csv" | tee -a "$LOGFILE"
echo "  Total size: $total_size" | tee -a "$LOGFILE"
echo "  Location: $OUTPUT_DIR" | tee -a "$LOGFILE"
echo "" | tee -a "$LOGFILE"

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
echo "3. Run analysis: cd analysis_pbc && python limits/run.py" | tee -a "$LOGFILE"
echo "============================================" | tee -a "$LOGFILE"
