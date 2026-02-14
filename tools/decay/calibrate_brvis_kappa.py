#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import List

import numpy as np

THIS_FILE = Path(__file__).resolve()
REPO_ROOT = THIS_FILE.parents[2]
ANALYSIS_ROOT = REPO_ROOT / "analysis_pbc"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(ANALYSIS_ROOT) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_ROOT))

from hnl_models.hnl_model_hnlcalc import HNLModel
from limits.expected_signal import couplings_from_eps2, expected_signal_events
from tools.decay.brvis_kappa_workflow import (
    FLAVOUR_TO_BENCHMARK,
    available_mass_points,
    load_geom_for_sim_files,
    parse_csv_list,
    resolve_flavours,
    resolve_masses,
    threshold_for_mass,
)

EPS2_REF = 1e-6
DEFAULT_OUT = REPO_ROOT / "output" / "csv" / "analysis" / "decay_kappa_table.csv"
DEFAULT_REPORT = REPO_ROOT / "output" / "csv" / "analysis" / "decay_kappa_validation.csv"


def _parse_seed_list(raw: str) -> List[int]:
    vals = [int(x) for x in parse_csv_list(raw)]
    if not vals:
        raise ValueError("--calib-seeds must contain at least one integer seed.")
    return vals


def _write_csv(path: Path, rows: List[dict], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Calibrate BR_vis*kappa decay surrogate against library-mode acceptance."
    )
    parser.add_argument("--flavours", type=str, default="electron,muon,tau")
    parser.add_argument("--masses", type=str, default=None, help="Comma-separated masses in GeV.")
    parser.add_argument(
        "--from-mass-grid",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Use MASS_GRID when --masses is not provided (default: true).",
    )
    parser.add_argument("--min-mass", type=float, default=0.0)
    parser.add_argument("--max-mass", type=float, default=17.0)
    parser.add_argument("--switch-mass", type=float, default=5.0)
    parser.add_argument("--p-min-gev", type=float, default=0.6)
    parser.add_argument("--separation-mm", type=float, default=1.0)
    parser.add_argument(
        "--max-signal-events",
        type=int,
        default=None,
        help="Optional cap on signal events per input file during calibration (default: no cap).",
    )
    parser.add_argument(
        "--allow-variant-drop",
        action="store_true",
        help="Allow dropping lower-priority pTHat/QCD variants instead of erroring.",
    )
    parser.add_argument("--calib-seeds", type=str, default="12345,22345,32345")
    parser.add_argument("--val-seed", type=int, default=42345)
    parser.add_argument("--out", type=str, default=str(DEFAULT_OUT))
    parser.add_argument("--report", type=str, default=str(DEFAULT_REPORT))
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    if args.switch_mass <= 0.0:
        raise ValueError("--switch-mass must be positive.")
    if args.p_min_gev <= 0.0:
        raise ValueError("--p-min-gev must be positive.")
    if args.separation_mm <= 0.0:
        raise ValueError("--separation-mm must be positive.")
    if args.max_signal_events is not None and args.max_signal_events <= 0:
        raise ValueError("--max-signal-events must be positive when provided.")

    flavours = resolve_flavours(args.flavours)
    masses = resolve_masses(args.masses, args.from_mass_grid, args.min_mass, args.max_mass)
    calib_seeds = _parse_seed_list(args.calib_seeds)

    source_policy = f"hybrid_external_lt{args.switch_mass:g}_generated_ge{args.switch_mass:g}"

    report_rows: List[dict] = []
    table_rows: List[dict] = []

    n_total = 0
    n_ok = 0

    for flavour in flavours:
        benchmark = FLAVOUR_TO_BENCHMARK[flavour]
        points = available_mass_points(
            flavour,
            selected_masses=masses,
            max_mass=args.max_mass,
            allow_variant_drop=args.allow_variant_drop,
        )
        if not args.quiet:
            print(f"[{flavour}] points={len(points)}")

        for mass_val, mass_str, sim_files in points:
            n_total += 1
            row = {
                "flavour": flavour,
                "mass_GeV": f"{mass_val:.6g}",
                "kappa": "",
                "p_min_GeV": f"{args.p_min_gev:.6g}",
                "separation_mm": f"{args.separation_mm:.6g}",
                "source_policy": source_policy,
                "status": "",
                "N_lib_mean": "",
                "N_brvis_unit": "",
                "N_lib_val": "",
                "N_pred": "",
                "rel_diff": "",
                "threshold": "",
                "error": "",
            }
            try:
                geom_df = load_geom_for_sim_files(
                    sim_files,
                    mass_str=mass_str,
                    flavour=flavour,
                    show_progress=False,
                    max_signal_events=args.max_signal_events,
                )
                if len(geom_df) == 0:
                    raise RuntimeError("No geometry rows loaded for point.")

                Ue2, Umu2, Utau2 = couplings_from_eps2(EPS2_REF, benchmark)
                model = HNLModel(mass_GeV=mass_val, Ue2=Ue2, Umu2=Umu2, Utau2=Utau2)
                ctau0_m = model.ctau0_m
                br_per_parent = model.production_brs()
                br_vis = model.visible_branching_ratio()

                n_brvis_unit = expected_signal_events(
                    geom_df=geom_df,
                    mass_GeV=mass_val,
                    eps2=EPS2_REF,
                    benchmark=benchmark,
                    lumi_fb=3000.0,
                    separation_m=args.separation_mm * 1e-3,
                    p_min_GeV=args.p_min_gev,
                    decay_mode="brvis_kappa",
                    br_vis=br_vis,
                    kappa_eff=1.0,
                    ctau0_m=ctau0_m,
                    br_per_parent=br_per_parent,
                    br_scale=1.0,
                )

                if n_brvis_unit <= 0.0:
                    raise RuntimeError(f"N_brvis(kappa=1) <= 0 ({n_brvis_unit:.6g}).")

                n_lib_samples: List[float] = []
                for seed in calib_seeds:
                    n_lib = expected_signal_events(
                        geom_df=geom_df,
                        mass_GeV=mass_val,
                        eps2=EPS2_REF,
                        benchmark=benchmark,
                        lumi_fb=3000.0,
                        separation_m=args.separation_mm * 1e-3,
                        decay_seed=seed,
                        p_min_GeV=args.p_min_gev,
                        decay_mode="library",
                        ctau0_m=ctau0_m,
                        br_per_parent=br_per_parent,
                        br_scale=1.0,
                    )
                    n_lib_samples.append(float(n_lib))

                n_lib_mean = float(np.mean(np.array(n_lib_samples, dtype=float)))
                kappa = n_lib_mean / float(n_brvis_unit)
                if not np.isfinite(kappa) or kappa <= 0.0:
                    raise RuntimeError(f"Invalid calibrated kappa={kappa}.")

                n_lib_val = expected_signal_events(
                    geom_df=geom_df,
                    mass_GeV=mass_val,
                    eps2=EPS2_REF,
                    benchmark=benchmark,
                    lumi_fb=3000.0,
                    separation_m=args.separation_mm * 1e-3,
                    decay_seed=args.val_seed,
                    p_min_GeV=args.p_min_gev,
                    decay_mode="library",
                    ctau0_m=ctau0_m,
                    br_per_parent=br_per_parent,
                    br_scale=1.0,
                )
                n_pred = float(n_brvis_unit) * float(kappa)

                threshold = threshold_for_mass(mass_val, args.switch_mass)
                rel_diff = abs(float(n_pred) - float(n_lib_val)) / max(abs(float(n_lib_val)), 1.0e-12)
                status = "ok" if rel_diff <= threshold else "fail"

                row.update(
                    {
                        "kappa": f"{kappa:.8g}",
                        "status": status,
                        "N_lib_mean": f"{n_lib_mean:.8g}",
                        "N_brvis_unit": f"{float(n_brvis_unit):.8g}",
                        "N_lib_val": f"{float(n_lib_val):.8g}",
                        "N_pred": f"{float(n_pred):.8g}",
                        "rel_diff": f"{rel_diff:.8g}",
                        "threshold": f"{threshold:.6g}",
                    }
                )
                if status == "ok":
                    table_rows.append(
                        {
                            "flavour": flavour,
                            "mass_GeV": f"{mass_val:.8g}",
                            "kappa": f"{kappa:.10g}",
                            "p_min_GeV": f"{args.p_min_gev:.8g}",
                            "separation_mm": f"{args.separation_mm:.8g}",
                            "source_policy": source_policy,
                            "status": "ok",
                        }
                    )
                    n_ok += 1
            except Exception as exc:
                row["status"] = "fail"
                row["error"] = str(exc)

            report_rows.append(row)
            if not args.quiet:
                print(
                    f"  m={mass_val:.3f} GeV status={row['status']} "
                    f"kappa={row['kappa'] or 'NA'} rel_diff={row['rel_diff'] or 'NA'}"
                )

    out_path = Path(args.out).resolve()
    report_path = Path(args.report).resolve()

    table_rows_sorted = sorted(table_rows, key=lambda r: (r["flavour"], float(r["mass_GeV"])))
    report_rows_sorted = sorted(report_rows, key=lambda r: (r["flavour"], float(r["mass_GeV"])))

    _write_csv(
        out_path,
        table_rows_sorted,
        [
            "flavour",
            "mass_GeV",
            "kappa",
            "p_min_GeV",
            "separation_mm",
            "source_policy",
            "status",
        ],
    )
    _write_csv(
        report_path,
        report_rows_sorted,
        [
            "flavour",
            "mass_GeV",
            "kappa",
            "p_min_GeV",
            "separation_mm",
            "source_policy",
            "status",
            "N_lib_mean",
            "N_brvis_unit",
            "N_lib_val",
            "N_pred",
            "rel_diff",
            "threshold",
            "error",
        ],
    )

    n_fail = n_total - n_ok
    print(
        f"Calibration complete: total={n_total} ok={n_ok} fail={n_fail} "
        f"table={out_path} report={report_path}"
    )

    if n_fail > 0:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
