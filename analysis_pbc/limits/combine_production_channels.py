#!/usr/bin/env python3

import pandas as pd
import re
import json
from pathlib import Path
import argparse
from collections import defaultdict


def find_production_files(sim_dir, flavour=None):
    pattern = re.compile(
        r"HNL_([0-9]+p[0-9]{1,2})GeV_([a-z]+)_((?:kaon|charm|beauty|Bc|ew)(?:_ff)?)(?:_(direct|fromTau))?(?:_(hardBc|hardccbar|hardbbbar)(?:_pTHat([0-9]+(?:p[0-9]+)?))?)?\.csv"
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
            qcd_mode = match.group(5) or "auto"
            pthat_token = match.group(6)
            pthat_min = float(pthat_token.replace("p", ".")) if pthat_token else None

            if flavour and file_flavour != flavour:
                continue

            mass_val = float(mass_str.replace('p', '.'))
            is_ff = regime.endswith("_ff")
            base_regime = regime.replace("_ff", "")
            files_by_mass[(mass_val, file_flavour)].append(
                (base_regime, mode, is_ff, qcd_mode, pthat_min, f)
            )

    return files_by_mass


def _variant_priority(base_regime: str, is_ff: bool, qcd_mode: str, pthat_min: float | None) -> tuple[int, int, float]:
    if base_regime == "charm" and qcd_mode == "hardccbar":
        qcd_priority = 3
    elif base_regime in {"beauty", "Bc"} and qcd_mode in {"hardbbbar", "hardBc"}:
        qcd_priority = 3
    elif qcd_mode != "auto":
        qcd_priority = 2
    else:
        qcd_priority = 1
    ff_priority = 1 if is_ff else 0
    pthat_priority = float(pthat_min) if pthat_min is not None else -1.0
    return (qcd_priority, ff_priority, pthat_priority)


def prefer_best_variant(regime_files, allow_variant_drop=False):
    chosen = {}
    all_candidates = {}
    for base_regime, mode, is_ff, qcd_mode, pthat_min, path in regime_files:
        key = (base_regime, mode)
        if key not in all_candidates:
            all_candidates[key] = []
        all_candidates[key].append((base_regime, mode, is_ff, qcd_mode, pthat_min, path))
        if key not in chosen:
            chosen[key] = (base_regime, mode, is_ff, qcd_mode, pthat_min, path)
            continue

        current = chosen[key]
        new_priority = _variant_priority(base_regime, is_ff, qcd_mode, pthat_min)
        old_priority = _variant_priority(current[0], current[2], current[3], current[4])
        if new_priority > old_priority:
            chosen[key] = (base_regime, mode, is_ff, qcd_mode, pthat_min, path)

    # Detect dropped variants (compare by path, not identity)
    for key, candidates in all_candidates.items():
        if len(candidates) > 1:
            kept_path = str(chosen[key][5])
            dropped = [c for c in candidates if str(c[5]) != kept_path]
            if dropped:
                kept_label = _format_source_label(*chosen[key][:5])
                dropped_labels = [_format_source_label(*d[:5]) for d in dropped]
                if allow_variant_drop:
                    print(
                        f"[DROP] {key}: keeping {kept_label}, "
                        f"dropping {dropped_labels}"
                    )
                else:
                    raise ValueError(
                        f"Multiple variants for {key}: keeping {kept_label}, "
                        f"would drop {dropped_labels}. "
                        f"Pass --allow-variant-drop to override."
                    )

    return list(chosen.values())


def _format_source_label(
    base_regime: str,
    mode: str | None,
    is_ff: bool,
    qcd_mode: str,
    pthat_min: float | None,
) -> str:
    label = base_regime
    if is_ff:
        label += "_ff"
    if mode:
        label += f"_{mode}"
    if qcd_mode != "auto":
        label += f"_{qcd_mode}"
        if pthat_min is not None:
            label += f"_pTHat{f'{pthat_min:g}'.replace('.', 'p')}"
    return label


def _load_sim_metadata(sim_csv: Path) -> dict:
    meta_path = Path(f"{sim_csv}.meta.json")
    if not meta_path.exists():
        return {}
    try:
        with meta_path.open("r", encoding="utf-8") as f:
            data = json.load(f)
            if isinstance(data, dict):
                return data
    except Exception as exc:
        print(f"[WARN] Failed to read metadata sidecar {meta_path.name}: {exc}")
    return {}


def combine_csvs(csv_paths, output_path):
    dfs = []

    for base_regime, mode, is_ff, qcd_mode, pthat_min, path in csv_paths:
        df = pd.read_csv(path)
        if "qcd_mode" not in df.columns or "sigma_gen_pb" not in df.columns or "pthat_min_gev" not in df.columns:
            meta = _load_sim_metadata(path)
            if "qcd_mode" not in df.columns:
                df["qcd_mode"] = str(meta.get("qcd_mode", qcd_mode))
            if "sigma_gen_pb" not in df.columns:
                sigma_pb = pd.to_numeric(pd.Series([meta.get("sigma_gen_pb")]), errors="coerce").iloc[0]
                df["sigma_gen_pb"] = float(sigma_pb) if pd.notna(sigma_pb) else float("nan")
            if "pthat_min_gev" not in df.columns:
                pthat_from_meta = pd.to_numeric(pd.Series([meta.get("pthat_min_gev")]), errors="coerce").iloc[0]
                if pd.notna(pthat_from_meta):
                    df["pthat_min_gev"] = float(pthat_from_meta)
                elif pthat_min is not None:
                    df["pthat_min_gev"] = float(pthat_min)
                else:
                    df["pthat_min_gev"] = float("nan")

        df["source_regime"] = base_regime
        df["source_mode"] = mode if mode is not None else "direct"
        df["source_is_ff"] = bool(is_ff)
        df["source_qcd_mode"] = qcd_mode
        df["source_pthat_min_gev"] = float(pthat_min) if pthat_min is not None else float("nan")
        dfs.append(df)
        print(f"    {_format_source_label(base_regime, mode, is_ff, qcd_mode, pthat_min):30s}: {len(df):6d} HNLs")

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
    parser.add_argument("--allow-variant-drop", action="store_true", help="Allow silently dropping lower-priority pTHat/QCD variants (default: error).")
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
        key: prefer_best_variant(regimes, allow_variant_drop=args.allow_variant_drop)
        for key, regimes in files_by_mass.items()
    }

    multi_channel_masses = {k: v for k, v in files_by_mass.items() if len(v) > 1}

    if not multi_channel_masses:
        print("✓ No overlapping production channels found")
        print("  (Each mass has only one production file)")
        return

    print(f"Found {len(multi_channel_masses)} masses with multiple production channels:\n")

    for (mass_val, flavour), regimes in sorted(multi_channel_masses.items()):
        regime_names = [_format_source_label(r, m, ff, qm, pt) for r, m, ff, qm, pt, _ in regimes]
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
            for base_regime, mode, is_ff, qcd_mode, pthat_min, fpath in csv_list:
                print(f"              - {_format_source_label(base_regime, mode, is_ff, qcd_mode, pthat_min)}: {fpath.name}")
            continue

        n_total = combine_csvs(csv_list, output_path)
        total_combined += 1

        for _, _, _, _, _, fpath in csv_list:
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
