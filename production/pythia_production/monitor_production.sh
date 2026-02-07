#!/bin/bash

OUTPUT_DIR="../../output/csv/simulation"
LOG_DIR="../../output/logs/simulation"

echo "============================================"
echo "PRODUCTION MONITORING"
echo "============================================"
echo ""

n_running=$(ps aux | grep "main_hnl_production" | grep -v grep | wc -l | tr -d ' ')
echo "Running jobs: $n_running"
echo ""

n_total=$(ls $OUTPUT_DIR/HNL_*.csv 2>/dev/null | wc -l | tr -d ' ')
echo "CSV files: $n_total"
echo ""

n_recent=$(find $OUTPUT_DIR -name "HNL_*.csv" -mmin -5 2>/dev/null | wc -l | tr -d ' ')
echo "Recent (last 5 min): $n_recent files"
echo ""

latest_log=$(ls -t $LOG_DIR/production_run_*.log 2>/dev/null | head -1)
if [ -n "$latest_log" ]; then
    echo "Latest progress:"
    tail -15 "$latest_log" | grep -E "Queued:|COMPLETE" | tail -5
    echo ""
    echo "Full log: $latest_log"
fi

echo ""
echo "============================================"
echo "To watch live: watch -n 10 ./monitor_production.sh"
echo "============================================"
