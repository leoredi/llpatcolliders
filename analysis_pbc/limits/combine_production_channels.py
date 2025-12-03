#!/usr/bin/env python3
"""
Combine Pythia (meson) and MadGraph (EW) production channels at overlapping masses.

At masses like 4-8 GeV, HNLs can be produced from BOTH:
1. B/D-meson decays (Pythia) → parent_pdg = 511, 521, 421, 431, etc.
2. W/Z-boson decays (MadGraph) → parent_pdg = 23, 24

These are DIFFERENT production mechanisms that should be ADDED, not double-counted.

This script:
- Identifies masses with multiple production files (e.g., beauty + ew)
- Combines CSV files at each mass into a single unified file
- Preserves all parent PDG codes for proper per-parent counting

Usage:
    python combine_production_channels.py [--flavour electron|muon|tau] [--dry-run]
"""

import pandas as pd
import re
from pathlib import Path
import argparse
from collections import defaultdict


def find_production_files(sim_dir, flavour=None):
    """
    Find all production CSV files and group by mass.

    Returns:
        dict: {mass_val: [(regime, filepath), ...]}
    """
    # Accept both 1 and 2 decimal formats during transition: 5p0 or 5p00
    pattern = re.compile(r"HNL_([0-9]+p[0-9]{1,2})GeV_([a-z]+)_(kaon|charm|beauty|ew)(?:_direct|_fromTau)?\.csv")

    files_by_mass = defaultdict(list)

    for f in sim_dir.glob("HNL_*.csv"):
        if f.stat().st_size < 1000:  # Skip empty
            continue

        match = pattern.search(f.name)
        if match:
            mass_str = match.group(1)
            file_flavour = match.group(2)
            regime = match.group(3)

            if flavour and file_flavour != flavour:
                continue

            mass_val = float(mass_str.replace('p', '.'))
            files_by_mass[(mass_val, file_flavour)].append((regime, f))

    return files_by_mass


def combine_csvs(csv_paths, output_path):
    """
    Combine multiple CSV files into one, preserving all columns.

    All files should have the same column structure (Pythia/MadGraph format):
    event, weight, hnl_id, parent_pdg, pt, eta, phi, p, E, mass, prod_x_mm, prod_y_mm, prod_z_mm, boost_gamma
    """
    dfs = []

    for regime, path in csv_paths:
        df = pd.read_csv(path)
        # Add regime column for tracking
        df['source_regime'] = regime
        dfs.append(df)
        print(f"    {regime:10s}: {len(df):6d} HNLs")

    combined = pd.concat(dfs, ignore_index=True)

    # Renumber events to avoid conflicts
    combined['event'] = range(len(combined))

    combined.to_csv(output_path, index=False)
    print(f"    → Combined: {len(combined):6d} HNLs → {output_path.name}")

    return len(combined)


def main():
    parser = argparse.ArgumentParser(description="Combine production channels at overlapping masses")
    parser.add_argument("--flavour", choices=["electron", "muon", "tau"], help="Process specific flavour only")
    parser.add_argument("--dry-run", action="store_true", help="Show what would be done without writing files")
    parser.add_argument("--no-cleanup", action="store_true", help="Don't move/backup files, just create combined versions")
    args = parser.parse_args()

    # Paths
    repo_root = Path(__file__).parent.parent.parent
    sim_dir = repo_root / "output" / "csv" / "simulation_new"
    combined_dir = sim_dir / "combined"
    backup_dir = sim_dir / "originals_backup"

    if not args.dry_run:
        combined_dir.mkdir(exist_ok=True)
        if not args.no_cleanup:
            backup_dir.mkdir(exist_ok=True)

    print("=" * 70)
    print("PRODUCTION CHANNEL COMBINATION")
    print("=" * 70)
    print(f"Input directory: {sim_dir}")
    if args.flavour:
        print(f"Flavour filter: {args.flavour}")
    if args.dry_run:
        print("DRY RUN MODE - no files will be written")
    print()

    # Find all files grouped by mass
    files_by_mass = find_production_files(sim_dir, args.flavour)

    # Find masses with multiple production channels
    multi_channel_masses = {k: v for k, v in files_by_mass.items() if len(v) > 1}

    if not multi_channel_masses:
        print("✓ No overlapping production channels found")
        print("  (Each mass has only one production file)")
        return

    print(f"Found {len(multi_channel_masses)} masses with multiple production channels:\n")

    for (mass_val, flavour), regimes in sorted(multi_channel_masses.items()):
        regime_names = [r for r, _ in regimes]
        print(f"  m = {mass_val:5.1f} GeV ({flavour}): {', '.join(regime_names)}")

    print(f"\n{'-' * 70}\n")

    # Combine files
    total_combined = 0
    files_to_backup = []

    for (mass_val, flavour), csv_list in sorted(multi_channel_masses.items()):
        # Use 2 decimal places for consistency (5.0 → 5p00, not 5p0)
        mass_str = f"{mass_val:.2f}".replace('.', 'p')
        output_path = combined_dir / f"HNL_{mass_str}GeV_{flavour}_combined.csv"

        print(f"Mass {mass_val} GeV ({flavour}):")

        if args.dry_run:
            print(f"    [DRY RUN] Would combine {len(csv_list)} files")
            for regime, fpath in csv_list:
                print(f"              - {fpath.name}")
            continue

        n_total = combine_csvs(csv_list, output_path)
        total_combined += 1

        # Track original files for backup
        for regime, fpath in csv_list:
            files_to_backup.append(fpath)

        print()

    # Cleanup: move combined files and backup originals
    if not args.dry_run and not args.no_cleanup:
        print(f"\n{'-' * 70}")
        print("CLEANUP: Moving files")
        print(f"{'-' * 70}\n")

        # Move combined files to main directory
        combined_files = list(combined_dir.glob("*.csv"))
        for f in combined_files:
            dest = sim_dir / f.name
            f.rename(dest)
            print(f"  ✓ {f.name} → {sim_dir.name}/")

        # Move original overlapping files to backup
        print()
        for f in files_to_backup:
            if f.exists():  # Check still exists (not already moved)
                dest = backup_dir / f.name
                f.rename(dest)
                print(f"  ✓ {f.name} → originals_backup/")

        # Remove empty combined directory
        if combined_dir.exists() and not any(combined_dir.iterdir()):
            combined_dir.rmdir()

    print("\n" + "=" * 70)
    if args.dry_run:
        print(f"DRY RUN: Would create {len(multi_channel_masses)} combined files")
        print(f"         Would backup {len(set(files_to_backup))} original files")
    elif args.no_cleanup:
        print(f"✓ Created {total_combined} combined files in: {combined_dir}")
        print(f"\nManual steps needed:")
        print(f"1. Move combined files to {sim_dir}")
        print(f"2. Backup original files to avoid double-counting")
    else:
        print(f"✓ Created {total_combined} combined files")
        print(f"✓ Moved to {sim_dir}")
        print(f"✓ Backed up {len(set(files_to_backup))} original files to {backup_dir}")
        print(f"\n✅ Ready to run analysis: python limits/run_serial.py")
    print("=" * 70)


if __name__ == "__main__":
    main()
