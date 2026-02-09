#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

OUTPUT_DIR="$REPO_ROOT/output/csv/simulation"
LOG_DIR="$REPO_ROOT/output/logs/simulation"

echo "============================================"
echo "PRODUCTION MONITORING"
echo "============================================"
echo ""

n_running=$(pgrep -fc main_hnl_production 2>/dev/null || echo 0)
echo "Running jobs: $n_running"
echo ""

n_total=$(find "$OUTPUT_DIR" -maxdepth 1 -name "HNL_*.csv" 2>/dev/null | wc -l | tr -d ' ')
echo "CSV files: $n_total"
echo ""

n_recent=$(find "$OUTPUT_DIR" -maxdepth 1 -name "HNL_*.csv" -mmin -5 2>/dev/null | wc -l | tr -d ' ')
echo "Recent (last 5 min): $n_recent files"
echo ""

latest_log=$(ls -t "$LOG_DIR"/production_run_*.log 2>/dev/null | head -1)
if [ -n "$latest_log" ]; then
    echo "Latest progress:"
    tail -15 "$latest_log" | grep -E "Queued:|COMPLETE" | tail -5
    echo ""
    echo "Full log: $latest_log"
fi

echo ""
echo "============================================"
echo "To watch live (from repo root): watch -n 10 ./tools/pythia/monitor_production.sh"
echo "============================================"
