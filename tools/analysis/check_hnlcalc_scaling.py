#!/usr/bin/env python3

from __future__ import annotations

import argparse
import math
import random
from typing import Dict, Iterable, Tuple

import sys
from pathlib import Path

THIS_FILE = Path(__file__).resolve()
REPO_ROOT = THIS_FILE.parents[2]
ANALYSIS_ROOT = REPO_ROOT / "analysis_pbc"
if str(ANALYSIS_ROOT) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_ROOT))

from hnl_models.hnl_model_hnlcalc import HNLModel


def _parse_list(text: str, cast):
    items = [t.strip() for t in text.split(",") if t.strip()]
    return [cast(x) for x in items]


def _couplings_from_eps2(eps2: float, flavour: str) -> Tuple[float, float, float]:
    if flavour == "electron":
        return eps2, 0.0, 0.0
    if flavour == "muon":
        return 0.0, eps2, 0.0
    if flavour == "tau":
        return 0.0, 0.0, eps2
    raise ValueError(f"Unsupported flavour: {flavour}")


def _relative_error(a: float, b: float) -> float:
    if b == 0.0:
        return math.inf if a != 0.0 else 0.0
    return abs(a - b) / abs(b)


def _br_stats(
    br_ref: Dict[int, float],
    br_new: Dict[int, float],
    scale: float,
    br_min: float,
) -> Tuple[float, float, int]:
    rel_errs = []
    keys = set(br_ref.keys()) & set(br_new.keys())
    for pid in keys:
        ref_val = br_ref.get(pid, 0.0)
        if abs(ref_val) < br_min:
            continue
        pred = ref_val * scale
        err = _relative_error(br_new.get(pid, 0.0), pred)
        rel_errs.append(err)
    if not rel_errs:
        return 0.0, 0.0, 0
    rel_errs.sort()
    mid = rel_errs[len(rel_errs) // 2]
    return max(rel_errs), mid, len(rel_errs)


def check_scaling(
    masses: Iterable[float],
    flavours: Iterable[str],
    eps2_list: Iterable[float],
    eps2_ref: float,
    tol: float,
    br_min: float,
    seed: int | None,
) -> int:
    failures = 0
    eps2_list = list(eps2_list)
    if eps2_ref not in eps2_list:
        eps2_list = [eps2_ref] + eps2_list

    for mass in masses:
        for flavour in flavours:
            print(f"\nMass {mass:.4f} GeV, flavour {flavour}")
            Ue2, Umu2, Utau2 = _couplings_from_eps2(eps2_ref, flavour)
            if seed is not None:
                random.seed(seed)
            model_ref = HNLModel(mass_GeV=mass, Ue2=Ue2, Umu2=Umu2, Utau2=Utau2)
            ctau_ref = float(model_ref.ctau0_m)
            br_ref = model_ref.production_brs()
            print(f"  eps2_ref={eps2_ref:.3e} ctau_ref={ctau_ref:.6e} m, parents={len(br_ref)}")

            for eps2 in eps2_list:
                if eps2 == eps2_ref:
                    continue
                Ue2, Umu2, Utau2 = _couplings_from_eps2(eps2, flavour)
                if seed is not None:
                    random.seed(seed)
                model = HNLModel(mass_GeV=mass, Ue2=Ue2, Umu2=Umu2, Utau2=Utau2)
                ctau = float(model.ctau0_m)
                brs = model.production_brs()

                ctau_pred = ctau_ref * (eps2_ref / eps2)
                ctau_err = _relative_error(ctau, ctau_pred)

                scale = eps2 / eps2_ref
                br_max_err, br_med_err, n_br = _br_stats(br_ref, brs, scale, br_min)

                status = "OK"
                if ctau_err > tol or br_max_err > tol:
                    status = "FAIL"
                    failures += 1

                print(
                    f"  eps2={eps2:.3e} ctau_err={ctau_err:.2e} "
                    f"br_max_err={br_max_err:.2e} br_med_err={br_med_err:.2e} "
                    f"(n_br={n_br}) [{status}]"
                )

    return failures


def main() -> int:
    parser = argparse.ArgumentParser(description="Check HNLCalc scaling with eps2.")
    parser.add_argument(
        "--masses",
        type=str,
        default="2.6",
        help="Comma-separated masses in GeV (default: 2.6)",
    )
    parser.add_argument(
        "--flavours",
        type=str,
        default="muon",
        help="Comma-separated flavours: electron,muon,tau (default: muon)",
    )
    parser.add_argument(
        "--eps2",
        type=str,
        default="1e-8,1e-6,1e-4",
        help="Comma-separated eps2 values to test (default: 1e-8,1e-6,1e-4)",
    )
    parser.add_argument(
        "--eps2-ref",
        type=float,
        default=1e-6,
        help="Reference eps2 used for scaling (default: 1e-6)",
    )
    parser.add_argument(
        "--tol",
        type=float,
        default=5e-4,
        help="Relative error tolerance for ctau and BR scaling (default: 5e-4)",
    )
    parser.add_argument(
        "--br-min",
        type=float,
        default=1e-15,
        help="Ignore BRs smaller than this when comparing (default: 1e-15)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Seed Python's random generator before each HNLCalc call for reproducible BRs.",
    )
    args = parser.parse_args()

    masses = _parse_list(args.masses, float)
    flavours = _parse_list(args.flavours, str)
    eps2_list = _parse_list(args.eps2, float)

    failures = check_scaling(
        masses=masses,
        flavours=flavours,
        eps2_list=eps2_list,
        eps2_ref=args.eps2_ref,
        tol=args.tol,
        br_min=args.br_min,
        seed=args.seed,
    )

    if failures:
        print(f"\nScaling check failed for {failures} case(s).")
        return 1
    print("\nScaling check passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
