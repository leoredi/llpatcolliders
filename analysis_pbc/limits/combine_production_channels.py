#!/usr/bin/env python3

import pandas as pd
import re
from pathlib import Path
import argparse
from collections import defaultdict


def find_production_files(sim_dir, flavour=None):
    pattern = re.compile(
        r"HNL_([0-9]+p[0-9]{1,2})GeV_([a-z]+)_((?:kaon|charm|beauty|ew)(?:_ff)?)(?:_(direct|fromTau))?\.csv"
    )

    files_by_mass = defaultdict(list)

    for f in sim_dir.glob("HNL_*.csv"):
        if f.stat().st_size < 1000:
            continue

        match = pattern.search(f.name)
        if match:
            mass_str = match.group(1)
            file_flavour = match.group(2)
            regime = match.group(3)
            mode = match.group(4)

            if flavour and file_flavour != flavour:
                continue

            mass_val = float(mass_str.replace('p', '.'))
            is_ff = regime.endswith("_ff")
            base_regime = regime.replace("_ff", "")
            files_by_mass[(mass_val, file_flavour)].append((base_regime, mode, is_ff, f))

    return files_by_mass


def prefer_ff(regime_files):
    chosen = {}
    for base_regime, mode, is_ff, path in regime_files:
        key = (base_regime, mode)
        keep_current = key not in chosen or is_ff
        if keep_current:
            chosen[key] = (base_regime, mode, is_ff, path)
    return list(chosen.values())


def _format_source_label(base_regime: str, mode: str | None, is_ff: bool) -> str:
    label = base_regime
    if is_ff:
        label += "_ff"
    if mode:
        label += f"_{mode}"
    return label


def combine_csvs(csv_paths, output_path):
    dfs = []

    for base_regime, mode, is_ff, path in csv_paths:
        df = pd.read_csv(path)
        df["source_regime"] = base_regime
        df["source_mode"] = mode if mode is not None else "direct"
        df["source_is_ff"] = bool(is_ff)
        dfs.append(df)
        print(f"    {_format_source_label(base_regime, mode, is_ff):16s}: {len(df):6d} HNLs")

    combined = pd.concat(dfs, ignore_index=True)

    combined['event'] = range(len(combined))

    combined.to_csv(output_path, index=False)
    print(f"    → Combined: {len(combined):6d} HNLs → {output_path.name}")

    return len(combined)


def main():
    parser = argparse.ArgumentParser(description="Combine production channels at overlapping masses")
    parser.add_argument("--flavour", choices=["electron", "muon", "tau"], help="Process specific flavour only")
    parser.add_argument("--dry-run", action="store_true", help="Show what would be done without writing files")
    parser.add_argument("--no-cleanup", action="store_true", help="Don't move combined files, just create them in a temporary folder")
    parser.add_argument("--delete-originals", action="store_true", help="Delete original files after combining (opt-in)")
    args = parser.parse_args()

    repo_root = Path(__file__).parent.parent.parent
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

    files_by_mass = find_production_files(sim_dir, args.flavour)
    files_by_mass = {
        key: prefer_ff(regimes) for key, regimes in files_by_mass.items()
    }

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

    total_combined = 0
    files_to_backup = []

    for (mass_val, flavour), csv_list in sorted(multi_channel_masses.items()):
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

        for _, _, _, fpath in csv_list:
            files_to_backup.append(fpath)

        print()

    if not args.dry_run and not args.no_cleanup:
        print(f"\n{'-' * 70}")
        print("CLEANUP: Moving combined files")
        print(f"{'-' * 70}\n")

        combined_files = list(combined_dir.glob("*.csv"))
        for f in combined_files:
            dest = sim_dir / f.name
            f.rename(dest)
            print(f"  ✓ {f.name} → {sim_dir.name}/")

        if args.delete_originals:
            print()
            for f in files_to_backup:
                if f.exists():
                    f.unlink()
                    print(f"  ✓ Deleted: {f.name}")

        if combined_dir.exists() and not any(combined_dir.iterdir()):
            combined_dir.rmdir()

    print("\n" + "=" * 70)
    if args.dry_run:
        print(f"DRY RUN: Would create {len(multi_channel_masses)} combined files")
        if args.delete_originals:
            print(f"         Would delete {len(set(files_to_backup))} original files")
    elif args.no_cleanup:
        print(f"✓ Created {total_combined} combined files in: {combined_dir}")
        print(f"\nManual steps needed:")
        print(f"1. Move combined files to {sim_dir}")
        print(f"2. (Optional) Delete original files to save space")
    else:
        print(f"✓ Created {total_combined} combined files")
        print(f"✓ Moved to {sim_dir}")
        if args.delete_originals:
            print(f"✓ Deleted {len(set(files_to_backup))} original files (data preserved in combined files)")
        else:
            print(f"• Kept {len(set(files_to_backup))} original files")
        print(f"\n Ready to run analysis: python limits/run.py")
    print("=" * 70)


if __name__ == "__main__":
    main()
