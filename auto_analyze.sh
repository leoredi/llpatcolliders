#!/bin/bash
# Auto-run analysis after mass scan completes

LOG_FILE="/tmp/auto_analysis.log"
echo "=== Auto-analysis started at $(date) ===" > $LOG_FILE

# Wait for mass scan to complete
echo "Waiting for mass scan to complete..." >> $LOG_FILE
while ps aux | grep -q "[r]un_mass_scan.py"; do
    sleep 60
done

echo "Mass scan completed at $(date)" >> $LOG_FILE
echo "Starting analysis on all mass points..." >> $LOG_FILE

# Activate conda environment and run analysis on all mass points
cd /Users/fredi/cernbox/Physics/llpatcolliders/llpatcolliders

MASSES=(15 23 31 39 47 55 63 71)

for mass in "${MASSES[@]}"; do
    CSV_FILE="output/csv/hnlLL_m${mass}GeVLLP.csv"

    if [ -f "$CSV_FILE" ]; then
        echo "[$(date '+%H:%M:%S')] Analyzing m=${mass} GeV..." >> $LOG_FILE
        conda run -n llpatcolliders python decayProbPerEvent.py "$CSV_FILE" >> $LOG_FILE 2>&1
        echo "[$(date '+%H:%M:%S')] Completed m=${mass} GeV" >> $LOG_FILE
    else
        echo "[$(date '+%H:%M:%S')] WARNING: $CSV_FILE not found, skipping" >> $LOG_FILE
    fi
done

echo "" >> $LOG_FILE
echo "=== All analysis completed at $(date) ===" >> $LOG_FILE
echo "" >> $LOG_FILE
echo "SUMMARY:" >> $LOG_FILE
echo "--------" >> $LOG_FILE
ls -lh output/images/hnlLL_m*GeVLLP_*.png 2>/dev/null | tail -20 >> $LOG_FILE

# Create completion marker
touch /tmp/analysis_complete.marker
echo "Analysis complete! Check /tmp/auto_analysis.log for details"
