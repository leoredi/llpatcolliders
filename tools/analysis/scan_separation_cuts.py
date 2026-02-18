#!/usr/bin/env python3

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Iterable

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RUN_SCRIPT = REPO_ROOT / "analysis_pbc" / "limits" / "run.py"
DEFAULT_SUMMARY = REPO_ROOT / "output" / "csv" / "analysis" / "HNL_U2_limits_summary.csv"
DEFAULT_SCAN_DIR = REPO_ROOT / "output" / "csv" / "analysis" / "separation_scan_runs"
DEFAULT_CONSOLIDATED = REPO_ROOT / "output" / "csv" / "analysis" / "separation_scan_results.csv"
DEFAULT_MANIFEST = REPO_ROOT / "output" / "csv" / "analysis" / "separation_scan_manifest.csv"


def _parse_float_list(raw: str) -> list[float]:
    vals = []
    for token in (x.strip() for x in raw.split(",")):
        if not token:
            continue
        vals.append(float(token))
    if not vals:
        raise ValueError("Expected a comma-separated list with at least one numeric value.")
    return vals


def _format_mm(val: float) -> str:
    return f"{float(val):g}".replace(".", "p")


def _iter_scan_points(
    separation_vals: Iterable[float],
    max_separation_vals: Iterable[float] | None,
) -> list[tuple[float, float | None]]:
    points: list[tuple[float, float | None]] = []
    if max_separation_vals is None:
        for sep in separation_vals:
            points.append((float(sep), None))
        return points

    for sep in separation_vals:
        for max_sep in max_separation_vals:
            if float(max_sep) <= float(sep):
                continue
            points.append((float(sep), float(max_sep)))
    return points


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Scan separation cut points by invoking analysis_pbc/limits/run.py as subprocesses. "
            "Produces a consolidated results CSV and a per-run manifest."
        )
    )
    parser.add_argument(
        "--separation-mm-values",
        type=str,
        default="1,10,100",
        help="Comma-separated list of minimum separation cuts in mm (default: 1,10,100).",
    )
    parser.add_argument(
        "--max-separation-mm-values",
        type=str,
        default=None,
        help=(
            "Optional comma-separated list of exploratory max separation cuts in mm. "
            "When set, scans all valid (min,max) pairs with max > min."
        ),
    )
    parser.add_argument(
        "--run-script",
        type=str,
        default=str(DEFAULT_RUN_SCRIPT),
        help=f"Path to limits run script (default: {DEFAULT_RUN_SCRIPT}).",
    )
    parser.add_argument(
        "--summary-path",
        type=str,
        default=str(DEFAULT_SUMMARY),
        help="Path where run.py writes the summary CSV.",
    )
    parser.add_argument(
        "--scan-dir",
        type=str,
        default=str(DEFAULT_SCAN_DIR),
        help="Directory where per-run summary copies are stored.",
    )
    parser.add_argument("--out", type=str, default=str(DEFAULT_CONSOLIDATED))
    parser.add_argument("--manifest", type=str, default=str(DEFAULT_MANIFEST))
    parser.add_argument(
        "--python",
        type=str,
        default=sys.executable,
        help="Python executable used to launch run.py.",
    )
    parser.add_argument(
        "--run-args",
        nargs=argparse.REMAINDER,
        default=[],
        help=(
            "Additional args forwarded to run.py. Example: "
            "--run-args --parallel --workers 12 --flavour electron"
        ),
    )
    args = parser.parse_args()

    separation_vals = _parse_float_list(args.separation_mm_values)
    if any(v <= 0.0 for v in separation_vals):
        raise ValueError("All --separation-mm-values must be positive.")

    max_vals = None
    if args.max_separation_mm_values is not None:
        max_vals = _parse_float_list(args.max_separation_mm_values)
        if any(v <= 0.0 for v in max_vals):
            raise ValueError("All --max-separation-mm-values must be positive.")

    scan_points = _iter_scan_points(separation_vals, max_vals)
    if not scan_points:
        raise ValueError("No valid scan points found (check min/max separation lists).")

    run_script = Path(args.run_script).resolve()
    summary_path = Path(args.summary_path).resolve()
    scan_dir = Path(args.scan_dir).resolve()
    out_path = Path(args.out).resolve()
    manifest_path = Path(args.manifest).resolve()
    scan_dir.mkdir(parents=True, exist_ok=True)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)

    manifest_rows: list[dict] = []
    result_frames: list[pd.DataFrame] = []

    for idx, (sep_mm, max_sep_mm) in enumerate(scan_points, start=1):
        run_id = f"scan_{idx:03d}_sep{_format_mm(sep_mm)}"
        if max_sep_mm is not None:
            run_id += f"_max{_format_mm(max_sep_mm)}"

        cmd = [args.python, str(run_script), "--separation-mm", f"{sep_mm:g}"]
        if max_sep_mm is not None:
            cmd.extend(["--max-separation-mm", f"{max_sep_mm:g}"])
        if args.run_args:
            cmd.extend(args.run_args)

        start = datetime.now(timezone.utc)
        proc = subprocess.run(
            cmd,
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
        )
        end = datetime.now(timezone.utc)
        duration_s = (end - start).total_seconds()

        summary_copy = ""
        status = "ok" if proc.returncode == 0 else "fail"
        error_msg = ""

        if proc.returncode == 0:
            if not summary_path.exists():
                status = "fail"
                error_msg = f"Expected summary file not found: {summary_path}"
            else:
                summary_copy_path = scan_dir / f"HNL_U2_limits_summary_{run_id}.csv"
                shutil.copy2(summary_path, summary_copy_path)
                summary_copy = str(summary_copy_path)
                df = pd.read_csv(summary_copy_path)
                df["scan_run_id"] = run_id
                df["scan_separation_mm"] = sep_mm
                df["scan_max_separation_mm"] = max_sep_mm
                result_frames.append(df)
        else:
            output = f"{proc.stdout}\n{proc.stderr}".strip()
            error_msg = output[-4000:]

        manifest_rows.append(
            {
                "run_id": run_id,
                "separation_mm": sep_mm,
                "max_separation_mm": max_sep_mm,
                "status": status,
                "returncode": proc.returncode,
                "start_utc": start.isoformat(),
                "end_utc": end.isoformat(),
                "duration_s": duration_s,
                "command": " ".join(cmd),
                "summary_copy": summary_copy,
                "error": error_msg,
            }
        )
        print(
            f"[{run_id}] status={status} returncode={proc.returncode} "
            f"duration={duration_s:.1f}s"
        )

    manifest_df = pd.DataFrame(manifest_rows)
    manifest_df.to_csv(manifest_path, index=False)
    print(f"Manifest written: {manifest_path}")

    if result_frames:
        combined = pd.concat(result_frames, ignore_index=True)
        combined.to_csv(out_path, index=False)
        print(f"Consolidated results written: {out_path}")
    else:
        print("No successful runs; consolidated results not written.")

    n_fail = int((manifest_df["status"] != "ok").sum()) if len(manifest_df) > 0 else 0
    if n_fail > 0:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
