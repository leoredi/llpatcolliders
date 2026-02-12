#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import os
import sys
import warnings
from pathlib import Path
from typing import List, Sequence, Tuple

THIS_FILE = Path(__file__).resolve()
REPO_ROOT = THIS_FILE.parents[2]
ANALYSIS_ROOT = REPO_ROOT / "analysis_pbc"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(ANALYSIS_ROOT) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_ROOT))

from config_mass_grid import MASS_GRID
from decay.rhn_decay_library import (
    FLAVOUR_CONFIG,
    MAX_DECAY_FILE_DELTA_GEV,
    DecayFileEntry,
    list_decay_files,
    load_decay_events,
)

ALL_FLAVOURS = ("electron", "muon", "tau")
INVISIBLE_NEUTRINO_PIDS = {12, 14, 16}


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
    masses = sorted({float(m) for m in masses})
    return [m for m in masses if args.min_mass <= m < args.max_mass]


def _allow_large_mismatch_from_env() -> bool:
    raw = os.environ.get("HNL_ALLOW_DECAY_MASS_MISMATCH", "").strip().lower()
    return raw in {"1", "true", "yes", "on"}


def _check_delta(source: str, flavour: str, requested_mass: float, selected: DecayFileEntry) -> Tuple[bool, float]:
    delta = abs(float(selected.mass_GeV) - float(requested_mass))
    if delta <= MAX_DECAY_FILE_DELTA_GEV:
        return True, delta
    message = (
        f"{source} mismatch at flavour={flavour} requested={requested_mass:.3f} GeV "
        f"selected={selected.mass_GeV:.3f} GeV delta={delta:.3f} GeV "
        f"(threshold={MAX_DECAY_FILE_DELTA_GEV:.3f}) path={selected.path}"
    )
    if _allow_large_mismatch_from_env():
        warnings.warn(f"{message} (allowed by HNL_ALLOW_DECAY_MASS_MISMATCH)", UserWarning)
        return True, delta
    return False, delta


def _nearest_for_source(flavour: str, source: str, mass_GeV: float) -> DecayFileEntry | None:
    entries = [e for e in list_decay_files(flavour) if e.source == source]
    if not entries:
        return None
    return min(
        entries,
        key=lambda e: (
            abs(float(e.mass_GeV) - float(mass_GeV)),
            e.path.name,
        ),
    )


def _mean_metrics(
    events: Sequence[Sequence[Tuple[float, float, float, float, float, int]]]
) -> Tuple[int, float, float]:
    if not events:
        return 0, 0.0, 0.0
    mean_daughters = 0.0
    mean_visible_frac = 0.0
    n_events = len(events)
    for event in events:
        mean_daughters += float(len(event))
        total_e = 0.0
        visible_e = 0.0
        for daughter in event:
            energy = float(daughter[0])
            pid = int(daughter[5])
            total_e += energy
            if abs(pid) not in INVISIBLE_NEUTRINO_PIDS:
                visible_e += energy
        event_visible_frac = 0.0 if total_e <= 0.0 else visible_e / total_e
        mean_visible_frac += event_visible_frac
    return n_events, mean_daughters / n_events, mean_visible_frac / n_events


def _relative_diff(a: float, b: float) -> float:
    denom = max(abs(float(b)), 1.0e-12)
    return abs(float(a) - float(b)) / denom


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Validate generated vs external decay-library overlap consistency."
    )
    parser.add_argument("--flavours", type=str, default="electron,muon,tau")
    parser.add_argument("--masses", type=str, default=None, help="Comma-separated mass list in GeV.")
    parser.add_argument(
        "--from-mass-grid",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Use MASS_GRID when --masses is not provided (default: true).",
    )
    parser.add_argument("--min-mass", type=float, default=4.0)
    parser.add_argument("--max-mass", type=float, default=5.0, help="Exclusive upper bound.")
    parser.add_argument("--max-visible-frac-diff", type=float, default=0.05)
    parser.add_argument("--max-mean-daughters-rel-diff", type=float, default=0.15)
    parser.add_argument("--out", type=str, default=None, help="Optional CSV report path.")
    args = parser.parse_args()

    if args.max_mass <= args.min_mass:
        raise ValueError("--max-mass must be greater than --min-mass.")
    if args.max_visible_frac_diff < 0.0:
        raise ValueError("--max-visible-frac-diff must be non-negative.")
    if args.max_mean_daughters_rel_diff < 0.0:
        raise ValueError("--max-mean-daughters-rel-diff must be non-negative.")

    flavours = _resolve_flavours(args.flavours)
    masses = _resolve_masses(args)

    rows: List[dict] = []
    n_fail = 0

    for flavour in flavours:
        threshold = float(FLAVOUR_CONFIG[flavour]["low_mass_threshold"])
        overlap_masses = [m for m in masses if m > threshold]
        for mass in overlap_masses:
            row = {
                "flavour": flavour,
                "requested_mass_GeV": f"{mass:.6g}",
                "external_mass_GeV": "",
                "external_delta_GeV": "",
                "external_path": "",
                "generated_mass_GeV": "",
                "generated_delta_GeV": "",
                "generated_path": "",
                "n_events_external": "",
                "n_events_generated": "",
                "mean_daughters_external": "",
                "mean_daughters_generated": "",
                "mean_visible_frac_external": "",
                "mean_visible_frac_generated": "",
                "visible_frac_abs_diff": "",
                "mean_daughters_rel_diff": "",
                "status": "",
                "error": "",
            }
            try:
                ext = _nearest_for_source(flavour, "external", mass)
                gen = _nearest_for_source(flavour, "generated", mass)
                if ext is None:
                    raise FileNotFoundError(
                        f"Missing external overlap candidate for flavour={flavour} mass={mass:.3f} GeV"
                    )
                if gen is None:
                    raise FileNotFoundError(
                        f"Missing generated overlap candidate for flavour={flavour} mass={mass:.3f} GeV"
                    )

                ext_ok, ext_delta = _check_delta("external", flavour, mass, ext)
                gen_ok, gen_delta = _check_delta("generated", flavour, mass, gen)
                if not ext_ok:
                    raise ValueError(
                        f"external strict mismatch for flavour={flavour} mass={mass:.3f} GeV: "
                        f"delta={ext_delta:.3f} GeV > {MAX_DECAY_FILE_DELTA_GEV:.3f} GeV"
                    )
                if not gen_ok:
                    raise ValueError(
                        f"generated strict mismatch for flavour={flavour} mass={mass:.3f} GeV: "
                        f"delta={gen_delta:.3f} GeV > {MAX_DECAY_FILE_DELTA_GEV:.3f} GeV"
                    )

                ext_events = load_decay_events(ext.path)
                gen_events = load_decay_events(gen.path)

                n_ext, mean_dau_ext, mean_vis_ext = _mean_metrics(ext_events)
                n_gen, mean_dau_gen, mean_vis_gen = _mean_metrics(gen_events)

                visible_frac_abs_diff = abs(mean_vis_gen - mean_vis_ext)
                mean_daughters_rel_diff = _relative_diff(mean_dau_gen, mean_dau_ext)

                row["external_mass_GeV"] = f"{ext.mass_GeV:.6g}"
                row["external_delta_GeV"] = f"{ext_delta:.6g}"
                row["external_path"] = str(ext.path)
                row["generated_mass_GeV"] = f"{gen.mass_GeV:.6g}"
                row["generated_delta_GeV"] = f"{gen_delta:.6g}"
                row["generated_path"] = str(gen.path)
                row["n_events_external"] = str(n_ext)
                row["n_events_generated"] = str(n_gen)
                row["mean_daughters_external"] = f"{mean_dau_ext:.6g}"
                row["mean_daughters_generated"] = f"{mean_dau_gen:.6g}"
                row["mean_visible_frac_external"] = f"{mean_vis_ext:.6g}"
                row["mean_visible_frac_generated"] = f"{mean_vis_gen:.6g}"
                row["visible_frac_abs_diff"] = f"{visible_frac_abs_diff:.6g}"
                row["mean_daughters_rel_diff"] = f"{mean_daughters_rel_diff:.6g}"

                visible_ok = visible_frac_abs_diff <= args.max_visible_frac_diff
                daughters_ok = mean_daughters_rel_diff <= args.max_mean_daughters_rel_diff
                if visible_ok and daughters_ok:
                    row["status"] = "ok"
                else:
                    status_parts = []
                    if not visible_ok:
                        status_parts.append("visible_frac_fail")
                    if not daughters_ok:
                        status_parts.append("mean_daughters_fail")
                    row["status"] = "+".join(status_parts)
                    row["error"] = (
                        f"visible_abs_diff={visible_frac_abs_diff:.6g} "
                        f"(max={args.max_visible_frac_diff:.6g}), "
                        f"mean_daughters_rel_diff={mean_daughters_rel_diff:.6g} "
                        f"(max={args.max_mean_daughters_rel_diff:.6g})"
                    )
                    n_fail += 1
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
                    "external_mass_GeV",
                    "external_delta_GeV",
                    "external_path",
                    "generated_mass_GeV",
                    "generated_delta_GeV",
                    "generated_path",
                    "n_events_external",
                    "n_events_generated",
                    "mean_daughters_external",
                    "mean_daughters_generated",
                    "mean_visible_frac_external",
                    "mean_visible_frac_generated",
                    "visible_frac_abs_diff",
                    "mean_daughters_rel_diff",
                    "status",
                    "error",
                ],
            )
            writer.writeheader()
            writer.writerows(rows)
        print(f"Wrote overlap validation report: {out_path}")

    n_ok = sum(r["status"] == "ok" for r in rows)
    print(
        "Overlap validation complete: "
        f"points={len(rows)} ok={n_ok} fail={n_fail} "
        f"window=[{args.min_mass:.3f},{args.max_mass:.3f}) GeV"
    )
    if n_fail > 0:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
