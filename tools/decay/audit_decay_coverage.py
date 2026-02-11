#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import List

THIS_FILE = Path(__file__).resolve()
REPO_ROOT = THIS_FILE.parents[2]
ANALYSIS_ROOT = REPO_ROOT / "analysis_pbc"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(ANALYSIS_ROOT) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_ROOT))

from config_mass_grid import MASS_GRID
from decay.rhn_decay_library import MAX_DECAY_FILE_DELTA_GEV, select_decay_file

ALL_FLAVOURS = ("electron", "muon", "tau")


def _parse_csv_list(raw: str) -> List[str]:
    return [x.strip() for x in raw.split(",") if x.strip()]


def _parse_mass_list(raw: str) -> List[float]:
    return [float(x) for x in _parse_csv_list(raw)]


def _resolve_flavours(raw: str) -> List[str]:
    flavours = _parse_csv_list(raw) if raw else list(ALL_FLAVOURS)
    if not flavours:
        return list(ALL_FLAVOURS)
    unknown = [f for f in flavours if f not in ALL_FLAVOURS]
    if unknown:
        raise ValueError(f"Unknown flavour(s): {unknown}. Allowed: {ALL_FLAVOURS}")
    return flavours


def _resolve_masses(args: argparse.Namespace) -> List[float]:
    if args.masses:
        masses = _parse_mass_list(args.masses)
    elif args.from_mass_grid:
        masses = list(MASS_GRID)
    else:
        raise ValueError("No masses selected. Provide --masses or keep --from-mass-grid enabled.")
    return sorted({float(m) for m in masses})


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Audit decay-library coverage and strict mass matching across masses/flavours."
    )
    parser.add_argument("--flavours", type=str, default="electron,muon,tau")
    parser.add_argument("--masses", type=str, default=None, help="Comma-separated mass list in GeV.")
    parser.add_argument(
        "--from-mass-grid",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Use MASS_GRID when --masses is not provided (default: true).",
    )
    parser.add_argument("--out", type=str, default=None, help="Optional CSV report path.")
    args = parser.parse_args()

    flavours = _resolve_flavours(args.flavours)
    masses = _resolve_masses(args)

    rows: List[dict] = []
    n_fail = 0

    for flavour in flavours:
        for mass in masses:
            row = {
                "flavour": flavour,
                "requested_mass_GeV": f"{mass:.6g}",
                "selected_mass_GeV": "",
                "delta_GeV": "",
                "source": "",
                "selected_path": "",
                "status": "",
                "error": "",
            }
            try:
                selected = select_decay_file(flavour, mass)
                delta = abs(float(selected.mass_GeV) - float(mass))
                row["selected_mass_GeV"] = f"{selected.mass_GeV:.6g}"
                row["delta_GeV"] = f"{delta:.6g}"
                row["source"] = selected.source
                row["selected_path"] = str(selected.path)
                if delta > MAX_DECAY_FILE_DELTA_GEV:
                    row["status"] = "strict_fail"
                    row["error"] = (
                        f"delta {delta:.6g} > threshold {MAX_DECAY_FILE_DELTA_GEV:.6g} "
                        "(selection allowed only via override env var)"
                    )
                    n_fail += 1
                else:
                    row["status"] = "ok"
            except Exception as exc:
                row["status"] = "missing_or_error"
                row["error"] = str(exc)
                n_fail += 1
            rows.append(row)

    if args.out:
        out_path = Path(args.out).resolve()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with out_path.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=[
                    "flavour",
                    "requested_mass_GeV",
                    "selected_mass_GeV",
                    "delta_GeV",
                    "source",
                    "selected_path",
                    "status",
                    "error",
                ],
            )
            writer.writeheader()
            writer.writerows(rows)
        print(f"Wrote coverage report: {out_path}")

    print(
        f"Audit complete: points={len(rows)} "
        f"ok={sum(r['status'] == 'ok' for r in rows)} "
        f"fail={n_fail} threshold={MAX_DECAY_FILE_DELTA_GEV:.3f} GeV"
    )
    if n_fail > 0:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

