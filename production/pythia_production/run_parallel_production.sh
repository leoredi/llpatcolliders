#!/bin/bash

set -e

FLAVOUR="${1:-all}"
FLAVOUR=$(echo "$FLAVOUR" | tr '[:upper:]' '[:lower:]')
MODE="${2:-both}"
MODE=$(echo "$MODE" | tr '[:upper:]' '[:lower:]')
QCD_MODE_RAW="${3:-auto}"
PTHAT_MIN="${4:-}"

if [[ ! "$FLAVOUR" =~ ^(electron|muon|tau|all)$ ]]; then
    echo "Error: Invalid flavour '$FLAVOUR'"
    echo "Usage: $0 [electron|muon|tau|all] [direct|fromTau|both] [auto|hardBc|hardccbar|hardbbbar] [pTHatMin]"
    exit 1
fi

if [[ ! "$MODE" =~ ^(direct|fromtau|both)$ ]]; then
    echo "Error: Invalid mode '$MODE'"
    echo "Usage: $0 [electron|muon|tau|all] [direct|fromTau|both] [auto|hardBc|hardccbar|hardbbbar] [pTHatMin]"
    exit 1
fi

if [[ "$MODE" == "fromtau" && "$FLAVOUR" != "tau" && "$FLAVOUR" != "all" ]]; then
    echo "Error: 'fromTau' mode only valid for tau flavour"
    echo "  fromTau uses τ → N X decay chain, which requires tau coupling"
    exit 1
fi

QCD_MODE_LOWER=$(echo "$QCD_MODE_RAW" | tr '[:upper:]' '[:lower:]')
case "$QCD_MODE_LOWER" in
    auto) QCD_MODE="auto" ;;
    hardbc) QCD_MODE="hardBc" ;;
    hardccbar) QCD_MODE="hardccbar" ;;
    hardbbbar) QCD_MODE="hardbbbar" ;;
    *)
        echo "Error: Invalid qcd mode '$QCD_MODE_RAW'"
        echo "Valid options: auto, hardBc, hardccbar, hardbbbar"
        exit 1
        ;;
esac

if [[ -n "$PTHAT_MIN" ]] && ! [[ "$PTHAT_MIN" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
    echo "Error: pTHatMin must be numeric, got '$PTHAT_MIN'"
    exit 1
fi

if [[ "$QCD_MODE" == "auto" && -n "$PTHAT_MIN" ]]; then
    echo "Error: pTHatMin is only valid with hard QCD modes (hardBc/hardccbar/hardbbbar)."
    exit 1
fi

echo "============================================"
echo "HNL Production - Parallel Execution"
echo "Flavour: $FLAVOUR"
echo "Mode: $MODE"
echo "QCD mode: $QCD_MODE"
if [[ -n "$PTHAT_MIN" ]]; then
    echo "pTHatMin: $PTHAT_MIN GeV"
fi
echo "============================================"
echo ""

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "$SCRIPT_DIR/load_mass_grid.sh"
NEVENTS=$N_EVENTS_DEFAULT
NEVENTS_FROMTAU=$NEVENTS
MAX_PARALLEL=12

FROMTAU_MASS_THRESHOLD=1.77
M_BC=6.275
M_ELECTRON=0.000511
M_MUON=0.10566
M_TAU=1.777

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
echo "  Events per mass: $NEVENTS (pp collisions)" | tee -a "$LOGFILE"
echo "  Max signal events: $MAX_SIGNAL_EVENTS (HNL cap)" | tee -a "$LOGFILE"
echo "  Parallel jobs: $MAX_PARALLEL" | tee -a "$LOGFILE"
echo "  QCD mode: $QCD_MODE" | tee -a "$LOGFILE"
if [[ -n "$PTHAT_MIN" ]]; then
    echo "  pTHatMin: $PTHAT_MIN GeV" | tee -a "$LOGFILE"
fi
echo "  Kinematic prefilter: direct mHNL < (mBc - m_l), fromTau mHNL < ${FROMTAU_MASS_THRESHOLD} GeV" | tee -a "$LOGFILE"
echo "  Output directory: $OUTPUT_DIR" | tee -a "$LOGFILE"
echo "  Log file: $LOGFILE" | tee -a "$LOGFILE"
echo ""

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

direct_mass_ceiling() {
    local flavour="$1"
    case "$flavour" in
        electron)
            echo "$M_BC - $M_ELECTRON" | bc -l
            ;;
        muon)
            echo "$M_BC - $M_MUON" | bc -l
            ;;
        tau)
            echo "$M_BC - $M_TAU" | bc -l
            ;;
        *)
            return 1
            ;;
    esac
}

can_run_point() {
    local mass="$1"
    local flavour="$2"
    local mode="$3"
    local mode_lower
    mode_lower=$(echo "$mode" | tr '[:upper:]' '[:lower:]')

    if [[ "$mode_lower" == "fromtau" ]]; then
        if (( $(echo "$mass < $FROMTAU_MASS_THRESHOLD" | bc -l) )); then
            return 0
        fi
        return 1
    fi

    local ceiling
    ceiling=$(direct_mass_ceiling "$flavour") || return 1
    if (( $(echo "$mass < $ceiling" | bc -l) )); then
        return 0
    fi
    return 1
}

skip_reason() {
    local flavour="$1"
    local mode="$2"
    local mode_lower
    mode_lower=$(echo "$mode" | tr '[:upper:]' '[:lower:]')

    if [[ "$mode_lower" == "fromtau" ]]; then
        echo "kinematically forbidden: fromTau requires mHNL < ${FROMTAU_MASS_THRESHOLD} GeV"
        return
    fi

    local ceiling
    ceiling=$(direct_mass_ceiling "$flavour")
    echo "kinematically forbidden: direct mode requires mHNL < ${ceiling} GeV"
}

mass_to_label() {
    local mass="$1"
    printf "%.2f" "$mass" | tr '.' 'p'
}

determine_regime() {
    local mass="$1"
    local flavour="$2"
    local qcd_mode="$3"

    if [[ "$qcd_mode" == "hardBc" ]]; then
        echo "Bc"
        return
    fi
    if [[ "$qcd_mode" == "hardccbar" ]]; then
        echo "charm"
        return
    fi
    if [[ "$qcd_mode" == "hardbbbar" ]]; then
        echo "beauty"
        return
    fi

    if [[ "$flavour" == "tau" ]]; then
        if (( $(echo "$mass < 2.0" | bc -l) )); then
            echo "charm"
        else
            echo "beauty"
        fi
        return
    fi

    if (( $(echo "$mass < 0.5" | bc -l) )); then
        echo "kaon"
    elif (( $(echo "$mass < 2.0" | bc -l) )); then
        echo "charm"
    else
        echo "beauty"
    fi
}

effective_pthat_min() {
    local qcd_mode="$1"
    local pthat_user="$2"

    if [[ "$qcd_mode" == "hardBc" ]]; then
        if [[ -n "$pthat_user" ]]; then
            echo "$pthat_user"
        else
            echo "15.0"
        fi
        return
    fi
    if [[ "$qcd_mode" == "hardccbar" || "$qcd_mode" == "hardbbbar" ]]; then
        if [[ -n "$pthat_user" ]]; then
            echo "$pthat_user"
        else
            echo "10.0"
        fi
        return
    fi
}

expected_output_csv() {
    local mass="$1"
    local flavour="$2"
    local mode="$3"
    local qcd_mode="$4"
    local pthat_user="$5"

    local mass_label
    mass_label=$(mass_to_label "$mass")
    local regime
    regime=$(determine_regime "$mass" "$flavour" "$qcd_mode")

    local name="HNL_${mass_label}GeV_${flavour}_${regime}"
    if [[ "$flavour" == "tau" ]]; then
        name+="_${mode}"
    fi

    if [[ "$qcd_mode" != "auto" ]]; then
        local pthat_eff
        pthat_eff=$(effective_pthat_min "$qcd_mode" "$pthat_user")
        name+="_${qcd_mode}"
        if [[ -n "$pthat_eff" ]]; then
            local pthat_label
            pthat_label=$(mass_to_label "$pthat_eff")
            name+="_pTHat${pthat_label}"
        fi
    fi

    echo "${name}.csv"
}

run_production_job() {
    local mass=$1
    local flavour=$2
    local mode=${3:-direct}

    local log_file="${LOG_DIR}/HNL_${mass}GeV_${flavour}_${mode}_${TIMESTAMP}.log"

    {
        echo "[$(date +%H:%M:%S)] Starting: $mass GeV $flavour ($mode)"

        local events=$NEVENTS
        if [ "$mode" = "fromTau" ]; then
            events=$NEVENTS_FROMTAU
        fi

        local cmd=(./main_hnl_production "${mass}" "${flavour}" "${events}" "${mode}")
        # Always pass qcdMode + pTHatMin so we can append maxSignalEvents
        cmd+=("$QCD_MODE")
        cmd+=("${PTHAT_MIN:--1}")
        cmd+=("$MAX_SIGNAL_EVENTS")
        "${cmd[@]}" 2>&1

        local exit_code=$?

        if [ $exit_code -eq 0 ]; then
            sleep 0.5
            local csv_file
            csv_file=$(expected_output_csv "$mass" "$flavour" "$mode" "$QCD_MODE" "$PTHAT_MIN")

            if [ -f "$csv_file" ]; then
                mv "$csv_file" "$OUTPUT_DIR/"
                local meta_file="${csv_file}.meta.json"
                if [ -f "$meta_file" ]; then
                    mv "$meta_file" "$OUTPUT_DIR/"
                fi
                echo "[$(date +%H:%M:%S)] SUCCESS: $csv_file → simulation/"
            else
                local mass_label
                mass_label=$(mass_to_label "$mass")
                local candidates
                candidates=$(ls HNL_${mass_label}GeV_${flavour}_*.csv 2>/dev/null | tr '\n' ' ')
                echo "[$(date +%H:%M:%S)] WARNING: Expected CSV not found: $csv_file"
                if [ -n "$candidates" ]; then
                    echo "[$(date +%H:%M:%S)] WARNING: Candidate files: $candidates"
                fi
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
            if can_run_point "$mass" "tau" "direct"; then
                count=$((count + 1))
            fi
        fi
        if [[ "$MODE" == "fromtau" || "$MODE" == "both" ]]; then
            if can_run_point "$mass" "tau" "fromTau"; then
                count=$((count + 1))
            fi
        fi
    done
    echo $count
}

count_direct_runs() {
    local flavour="$1"
    shift
    local count=0
    local mass
    for mass in "$@"; do
        if can_run_point "$mass" "$flavour" "direct"; then
            count=$((count + 1))
        fi
    done
    echo "$count"
}

total_electron=$(count_direct_runs "electron" "${ELECTRON_MASSES[@]}")
total_muon=$(count_direct_runs "muon" "${MUON_MASSES[@]}")
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
        if can_run_point "$mass" "electron" "direct"; then
            wait_for_slot
            run_production_job "$mass" "electron" "direct" &
            completed_jobs=$((completed_jobs + 1))
            echo "[$completed_jobs/$total_jobs] Queued: $mass GeV electron" | tee -a "$LOGFILE"
        else
            echo "[SKIP] $mass GeV electron (direct): $(skip_reason "electron" "direct")" | tee -a "$LOGFILE"
        fi
    done
fi


if [[ "$FLAVOUR" == "muon" || "$FLAVOUR" == "all" ]]; then
    echo "" | tee -a "$LOGFILE"
    echo "============================================" | tee -a "$LOGFILE"
    echo "MUON COUPLING (BC7)" | tee -a "$LOGFILE"
    echo "============================================" | tee -a "$LOGFILE"

    for mass in "${MUON_MASSES[@]}"; do
        if can_run_point "$mass" "muon" "direct"; then
            wait_for_slot
            run_production_job "$mass" "muon" "direct" &
            completed_jobs=$((completed_jobs + 1))
            echo "[$completed_jobs/$total_jobs] Queued: $mass GeV muon" | tee -a "$LOGFILE"
        else
            echo "[SKIP] $mass GeV muon (direct): $(skip_reason "muon" "direct")" | tee -a "$LOGFILE"
        fi
    done
fi


if [[ "$FLAVOUR" == "tau" || "$FLAVOUR" == "all" ]]; then
    echo "" | tee -a "$LOGFILE"
    echo "============================================" | tee -a "$LOGFILE"
    echo "TAU COUPLING (BC8) - MODE: $MODE" | tee -a "$LOGFILE"
    echo "============================================" | tee -a "$LOGFILE"

    for mass in "${TAU_MASSES[@]}"; do
        if [[ "$MODE" == "direct" || "$MODE" == "both" ]]; then
            if can_run_point "$mass" "tau" "direct"; then
                wait_for_slot
                run_production_job "$mass" "tau" "direct" &
                completed_jobs=$((completed_jobs + 1))
                echo "[$completed_jobs/$total_jobs] Queued: $mass GeV tau (direct)" | tee -a "$LOGFILE"
            else
                echo "[SKIP] $mass GeV tau (direct): $(skip_reason "tau" "direct")" | tee -a "$LOGFILE"
            fi
        fi

        if [[ "$MODE" == "fromtau" || "$MODE" == "both" ]]; then
            if can_run_point "$mass" "tau" "fromTau"; then
                wait_for_slot
                run_production_job "$mass" "tau" "fromTau" &
                completed_jobs=$((completed_jobs + 1))
                echo "[$completed_jobs/$total_jobs] Queued: $mass GeV tau (fromTau)" | tee -a "$LOGFILE"
            else
                echo "[SKIP] $mass GeV tau (fromTau): $(skip_reason "tau" "fromTau")" | tee -a "$LOGFILE"
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
