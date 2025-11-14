#!/usr/bin/env python3
"""
Quick script to check for duplicate HNL entries in mass scan CSV files
"""

import pandas as pd
import os

csv_dir = "pythiaStuff/mass_scan_hnl"
csv_files = sorted([f for f in os.listdir(csv_dir) if f.endswith("LLP.csv")])

print("=" * 70)
print("DUPLICATE CHECK FOR MASS SCAN CSV FILES")
print("=" * 70)

for csv_file in csv_files:
    filepath = os.path.join(csv_dir, csv_file)
    df = pd.read_csv(filepath, sep=',')

    # Strip whitespace from column names
    df.columns = df.columns.str.strip()

    # Extract mass from filename
    mass = csv_file.split('_m')[1].split('GeV')[0]

    # Check total entries
    total_entries = len(df)

    # Check for exact duplicates across all columns
    duplicates = df.duplicated()
    n_duplicates = duplicates.sum()

    # Check for kinematic duplicates within same event
    kinematic_cols = ['event', 'pt', 'eta', 'phi', 'momentum']
    kinematic_duplicates = df.duplicated(subset=kinematic_cols)
    n_kinematic_duplicates = kinematic_duplicates.sum()

    # Count unique events
    n_events = df['event'].nunique()

    # Count particles per event
    particles_per_event = df.groupby('event').size()
    events_with_1 = (particles_per_event == 1).sum()
    events_with_2 = (particles_per_event == 2).sum()
    events_with_more = (particles_per_event > 2).sum()

    print(f"\n{csv_file} (m = {mass} GeV):")
    print(f"  Total entries: {total_entries:,}")
    print(f"  Unique events: {n_events:,}")
    print(f"  Exact duplicates: {n_duplicates:,}")
    print(f"  Kinematic duplicates: {n_kinematic_duplicates:,}")
    print(f"  Events with 1 HNL: {events_with_1:,} ({events_with_1/n_events*100:.1f}%)")
    print(f"  Events with 2 HNLs: {events_with_2:,} ({events_with_2/n_events*100:.1f}%)")
    if events_with_more > 0:
        print(f"  Events with >2 HNLs: {events_with_more:,} ({events_with_more/n_events*100:.1f}%)")
        print(f"    WARNING: Unexpected duplicates detected!")

print("\n" + "=" * 70)
print("SUMMARY: All CSV files checked")
print("=" * 70 + "\n")
