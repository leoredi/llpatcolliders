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

from decay.brvis_kappa import lookup_kappa, resolve_kappa_table_path
from hnl_models.hnl_model_hnlcalc import HNLModel
from limits.expected_signal import couplings_from_eps2, expected_signal_events
from tools.decay.brvis_kappa_workflow import (
    FLAVOUR_TO_BENCHMARK,
    available_mass_points,
    resolve_flavours,
    resolve_masses,
    threshold_for_mass,
)

EPS2_REF = 1e-6
DEFAULT_OUT = REPO_ROOT / "output" / "csv" / "analysis" / "decay_kappa_validation_check.csv"


def _write_csv(path: Path, rows: List[dict], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Validate calibrated BR_vis*kappa mode against library mode on selected points."
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
    parser.add_argument("--val-seed", type=int, default=42345)
    parser.add_argument("--kappa-table", type=str, default=None)
    parser.add_argument("--out", type=str, default=str(DEFAULT_OUT))
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    if args.switch_mass <= 0.0:
        raise ValueError("--switch-mass must be positive.")
    if args.p_min_gev <= 0.0:
        raise ValueError("--p-min-gev must be positive.")
    if args.separation_mm <= 0.0:
        raise ValueError("--separation-mm must be positive.")

    flavours = resolve_flavours(args.flavours)
    masses = resolve_masses(args.masses, args.from_mass_grid, args.min_mass, args.max_mass)
    kappa_table_path = resolve_kappa_table_path(args.kappa_table)

    rows: List[dict] = []
    n_total = 0
    n_fail = 0

    for flavour in flavours:
        benchmark = FLAVOUR_TO_BENCHMARK[flavour]
        points = available_mass_points(flavour, selected_masses=masses, max_mass=args.max_mass)
        if not args.quiet:
            print(f"[{flavour}] points={len(points)}")

        for mass_val, mass_str, sim_files in points:
            n_total += 1
            row = {
                "flavour": flavour,
                "mass_GeV": f"{mass_val:.6g}",
                "N_lib": "",
                "N_pred": "",
                "rel_diff": "",
                "threshold": "",
                "status": "",
                "error": "",
            }
            try:
                from tools.decay.brvis_kappa_workflow import load_geom_for_sim_files

                geom_df = load_geom_for_sim_files(sim_files, mass_str=mass_str, flavour=flavour, show_progress=False)
                if len(geom_df) == 0:
                    raise RuntimeError("No geometry rows loaded for point.")

                Ue2, Umu2, Utau2 = couplings_from_eps2(EPS2_REF, benchmark)
                model = HNLModel(mass_GeV=mass_val, Ue2=Ue2, Umu2=Umu2, Utau2=Utau2)
                ctau0_m = model.ctau0_m
                br_per_parent = model.production_brs()
                br_vis = model.visible_branching_ratio()

                kappa_eff = lookup_kappa(
                    flavour=flavour,
                    mass_GeV=mass_val,
                    p_min_GeV=args.p_min_gev,
                    separation_mm=args.separation_mm,
                    table_path=kappa_table_path,
                )

                n_lib = expected_signal_events(
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
                n_pred = expected_signal_events(
                    geom_df=geom_df,
                    mass_GeV=mass_val,
                    eps2=EPS2_REF,
                    benchmark=benchmark,
                    lumi_fb=3000.0,
                    separation_m=args.separation_mm * 1e-3,
                    p_min_GeV=args.p_min_gev,
                    decay_mode="brvis_kappa",
                    br_vis=br_vis,
                    kappa_eff=kappa_eff,
                    ctau0_m=ctau0_m,
                    br_per_parent=br_per_parent,
                    br_scale=1.0,
                )

                threshold = threshold_for_mass(mass_val, args.switch_mass)
                rel_diff = abs(float(n_pred) - float(n_lib)) / max(abs(float(n_lib)), 1.0e-12)
                status = "ok" if rel_diff <= threshold else "fail"

                row.update(
                    {
                        "N_lib": f"{float(n_lib):.8g}",
                        "N_pred": f"{float(n_pred):.8g}",
                        "rel_diff": f"{rel_diff:.8g}",
                        "threshold": f"{threshold:.6g}",
                        "status": status,
                    }
                )
                if status != "ok":
                    n_fail += 1
            except Exception as exc:
                row["status"] = "fail"
                row["error"] = str(exc)
                n_fail += 1

            rows.append(row)
            if not args.quiet:
                print(
                    f"  m={mass_val:.3f} GeV status={row['status']} rel_diff={row['rel_diff'] or 'NA'}"
                )

    out_path = Path(args.out).resolve()
    rows_sorted = sorted(rows, key=lambda r: (r["flavour"], float(r["mass_GeV"])))
    _write_csv(
        out_path,
        rows_sorted,
        ["flavour", "mass_GeV", "N_lib", "N_pred", "rel_diff", "threshold", "status", "error"],
    )

    print(
        f"Validation complete: points={n_total} fail={n_fail} "
        f"kappa_table={kappa_table_path} out={out_path}"
    )

    if n_fail > 0:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
