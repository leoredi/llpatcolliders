#!/usr/bin/env python3
"""
Combine Pythia (meson) and MadGraph (EW) production channels at overlapping masses.

**IMPORTANT: Run this BEFORE analysis to avoid undercounting signal!**

At masses like 4-8 GeV, HNLs can be produced from BOTH:
1. B/D-meson decays (Pythia) → parent_pdg = 511, 521, 421, 431, etc.
2. W/Z-boson decays (MadGraph) → parent_pdg = 23, 24

These are DIFFERENT production mechanisms that should be ADDED, not double-counted.

This script:
- Identifies masses with multiple production files (e.g., beauty + ew)
- Combines CSV files at each mass into a single unified file
- Preserves all parent PDG codes for proper per-parent counting
- DELETES original separate files (data preserved in combined files)
- Saves ~2 GB of disk space by removing duplicates

Why this is critical:
- WITHOUT combining: Analysis uses only one production channel → undercounts signal
- WITH combining: Analysis includes all production mechanisms → correct sensitivity

Usage:
    # Combine all flavors (recommended)
    python combine_production_channels.py

    # Single flavor only
    python combine_production_channels.py --flavour electron

    # Preview changes without writing files
    python combine_production_channels.py --dry-run

    # Create combined files but don't delete originals
    python combine_production_channels.py --no-cleanup

Typical workflow:
    1. Run Pythia production (kaon/charm/beauty regimes)
    2. Run MadGraph production (EW regime)
    3. **Run this script** ← YOU ARE HERE
    4. Run analysis: python limits/run.py --parallel
    5. Generate plots: python ../money_plot/plot_money_island.py
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
        dict: {(mass_val, flavour): [(base_regime, mode, is_ff, filepath), ...]}
    """
    # Accept both 1 and 2 decimal formats during transition: 5p0 or 5p00
    pattern = re.compile(
        r"HNL_([0-9]+p[0-9]{1,2})GeV_([a-z]+)_((?:kaon|charm|beauty|ew)(?:_ff)?)(?:_(direct|fromTau))?\.csv"
    )

    files_by_mass = defaultdict(list)

    for f in sim_dir.glob("HNL_*.csv"):
        if f.stat().st_size < 1000:  # Skip empty
            continue

        match = pattern.search(f.name)
        if match:
            mass_str = match.group(1)
            file_flavour = match.group(2)
            regime = match.group(3)  # e.g., beauty or beauty_ff
            mode = match.group(4)    # direct/fromTau (tau only) or None

            if flavour and file_flavour != flavour:
                continue

            mass_val = float(mass_str.replace('p', '.'))
            is_ff = regime.endswith("_ff")
            base_regime = regime.replace("_ff", "")
            files_by_mass[(mass_val, file_flavour)].append((base_regime, mode, is_ff, f))

    return files_by_mass


def prefer_ff(regime_files):
    """
    If both base and *_ff versions exist for the same regime, keep only the *_ff version.
    regime_files: list of (base_regime, mode, is_ff, path)
    Returns filtered list.
    """
    chosen = {}
    for base_regime, mode, is_ff, path in regime_files:
        key = (base_regime, mode)
        keep_current = key not in chosen or is_ff
        if keep_current:
            chosen[key] = (base_regime, mode, is_ff, path)
    return list(chosen.values())


def normalize_boost_column(df):
    """
    Normalize boost factor column naming for backward compatibility.

    Production CSVs generated before Dec 2025 used 'boost_gamma'.
    Current production uses 'beta_gamma' (physically correct: βγ = p/m, not γ = E/m).
    This function ensures both formats are accepted.

    Args:
        df: DataFrame loaded from production CSV

    Returns:
        DataFrame with 'beta_gamma' column (renames 'boost_gamma' if present)
    """
    if 'boost_gamma' in df.columns and 'beta_gamma' not in df.columns:
        df = df.rename(columns={'boost_gamma': 'beta_gamma'})
    elif 'boost_gamma' in df.columns and 'beta_gamma' in df.columns:
        # Both present - drop the legacy column
        df = df.drop(columns=['boost_gamma'])
    # If only beta_gamma present, no action needed
    return df


def _format_source_label(base_regime: str, mode: str | None, is_ff: bool) -> str:
    label = base_regime
    if is_ff:
        label += "_ff"
    if mode:
        label += f"_{mode}"
    return label


def combine_csvs(csv_paths, output_path):
    """
    Combine multiple CSV files into one, preserving all columns.

    All files should have the same column structure (Pythia/MadGraph format):
    event, weight, hnl_id, parent_pdg, pt, eta, phi, p, E, mass, prod_x_mm, prod_y_mm, prod_z_mm, beta_gamma

    Note: Legacy files using 'boost_gamma' column name are automatically converted to 'beta_gamma'.
    """
    dfs = []

    for base_regime, mode, is_ff, path in csv_paths:
        df = pd.read_csv(path)
        # Normalize column naming for backward compatibility
        df = normalize_boost_column(df)
        # Add provenance columns for tracking
        df["source_regime"] = base_regime
        df["source_mode"] = mode if mode is not None else ""
        df["source_is_ff"] = bool(is_ff)
        dfs.append(df)
        print(f"    {_format_source_label(base_regime, mode, is_ff):16s}: {len(df):6d} HNLs")

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
    # Current simulation output directory (Pythia + MadGraph)
    sim_dir = repo_root / "output" / "csv" / "simulation"
    combined_dir = sim_dir / "combined"

    if not args.dry_run:
        combined_dir.mkdir(exist_ok=True)

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
    # Prefer *_ff replacements over base regimes to avoid double counting
    files_by_mass = {
        key: prefer_ff(regimes) for key, regimes in files_by_mass.items()
    }

    # Find masses with multiple production channels
    multi_channel_masses = {k: v for k, v in files_by_mass.items() if len(v) > 1}

    if not multi_channel_masses:
        print("✓ No overlapping production channels found")
        print("  (Each mass has only one production file)")
        return

    print(f"Found {len(multi_channel_masses)} masses with multiple production channels:\n")

    for (mass_val, flavour), regimes in sorted(multi_channel_masses.items()):
        regime_names = [_format_source_label(r, m, ff) for r, m, ff, _ in regimes]
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
            for base_regime, mode, is_ff, fpath in csv_list:
                print(f"              - {_format_source_label(base_regime, mode, is_ff)}: {fpath.name}")
            continue

        n_total = combine_csvs(csv_list, output_path)
        total_combined += 1

        # Track original files for deletion
        for _, _, _, fpath in csv_list:
            files_to_backup.append(fpath)

        print()

    # Cleanup: move combined files and delete originals
    if not args.dry_run and not args.no_cleanup:
        print(f"\n{'-' * 70}")
        print("CLEANUP: Moving files and removing duplicates")
        print(f"{'-' * 70}\n")

        # Move combined files to main directory
        combined_files = list(combined_dir.glob("*.csv"))
        for f in combined_files:
            dest = sim_dir / f.name
            f.rename(dest)
            print(f"  ✓ {f.name} → {sim_dir.name}/")

        # Delete original overlapping files (they're now combined)
        print()
        for f in files_to_backup:
            if f.exists():  # Check still exists (not already moved)
                f.unlink()
                print(f"  ✓ Deleted: {f.name}")

        # Remove empty combined directory
        if combined_dir.exists() and not any(combined_dir.iterdir()):
            combined_dir.rmdir()

    print("\n" + "=" * 70)
    if args.dry_run:
        print(f"DRY RUN: Would create {len(multi_channel_masses)} combined files")
        print(f"         Would delete {len(set(files_to_backup))} original files")
    elif args.no_cleanup:
        print(f"✓ Created {total_combined} combined files in: {combined_dir}")
        print(f"\nManual steps needed:")
        print(f"1. Move combined files to {sim_dir}")
        print(f"2. Delete original files to avoid double-counting")
    else:
        print(f"✓ Created {total_combined} combined files")
        print(f"✓ Moved to {sim_dir}")
        print(f"✓ Deleted {len(set(files_to_backup))} original files (data preserved in combined files)")
        print(f"\n✅ Ready to run analysis: python limits/run.py")
    print("=" * 70)


if __name__ == "__main__":
    main()
