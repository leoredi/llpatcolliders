#!/usr/bin/env python3
"""
Run decay probability analysis for all mass points
"""

import subprocess
import os
from datetime import datetime

csv_dir = "pythiaStuff/mass_scan_hnl"
csv_files = sorted([f for f in os.listdir(csv_dir) if f.endswith("LLP.csv")])

print("=" * 70)
print("DECAY PROBABILITY ANALYSIS FOR ALL MASS POINTS")
print("=" * 70)
print(f"Found {len(csv_files)} CSV files to process\n")

start_time = datetime.now()

for i, csv_file in enumerate(csv_files, 1):
    filepath = os.path.join(csv_dir, csv_file)
    mass = csv_file.split('_m')[1].split('GeV')[0]

    print(f"\n[{datetime.now().strftime('%H:%M:%S')}] Processing {i}/{len(csv_files)}: m = {mass} GeV")
    print(f"  File: {csv_file}")

    file_start = datetime.now()

    try:
        # Run the decay probability analysis
        result = subprocess.run(
            ['python', 'decayProbPerEvent.py', filepath],
            capture_output=True,
            text=True,
            timeout=600  # 10 minute timeout
        )

        file_end = datetime.now()
        duration = (file_end - file_start).total_seconds()

        if result.returncode == 0:
            print(f"  ✓ Completed in {duration:.1f}s")

            # Look for output plot
            expected_plot = filepath.replace('.csv', '_exclusion.png')
            if os.path.exists(expected_plot):
                print(f"  Plot saved: {os.path.basename(expected_plot)}")
        else:
            print(f"  ✗ Failed ({duration:.1f}s)")
            print(f"  Error: {result.stderr[:200]}")

    except subprocess.TimeoutExpired:
        print(f"  ✗ Timeout after 10 minutes")
    except Exception as e:
        print(f"  ✗ Error: {e}")

end_time = datetime.now()
total_duration = (end_time - start_time).total_seconds()

print(f"\n{'=' * 70}")
print(f"ANALYSIS COMPLETE")
print(f"{'=' * 70}")
print(f"Total time: {total_duration/60:.1f} minutes")
print(f"Average per file: {total_duration/len(csv_files):.1f}s")
print(f"{'=' * 70}\n")
