#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "output"
sys.path.insert(0, str(ROOT / "scripts"))
from generate_meson_grids import fit_charm_d0_feeddown_weight, fonll_env  # noqa: E402

# Public FONLL web v1.3.2, CTEQ6.6, pp 14 TeV, central prediction,
# dsigma/dpT/dy at pT=5 GeV, y=0, fragmentation fraction 1.
# Queried from https://www.lpthe.jussieu.fr/~cacciari/fonll/fonllform.html
# on 2026-05-29.
PUBLIC_POINTS = {
    "bottom": {
        "filename": "fonll_pp14tev_cteq66_fonll_meson_dsdpTdy_pt0-50_y-3to3_central_bottom.dat",
        "pt": 5.0,
        "y": 0.0,
        "public_pb_per_GeV": 7.5329e6,
        "public_fragmentation": "B hadron default, N=5 Kartvelishvili alpha=24.2",
    },
    "charm": {
        "pt": 5.0,
        "y": 0.0,
        "public_pb_per_GeV": 3.5221e7,
        "public_fragmentation": "D0 default, BCFY r=0.1 including public D* feeddown convention",
    },
}


def find_point(path: Path, pt: float, y: float) -> float:
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        cols = line.split()
        if len(cols) < 3:
            continue
        row_pt = float(cols[0].replace("D", "E"))
        row_y = float(cols[1].replace("D", "E"))
        if math.isclose(row_pt, pt, abs_tol=1e-12) and math.isclose(row_y, y, abs_tol=1e-12):
            return float(cols[2].replace("D", "E"))
    raise ValueError(f"point pT={pt}, y={y} not found in {path}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tolerance", type=float, default=0.05, help="maximum allowed relative difference")
    args = parser.parse_args()

    failed = False
    for quark, cfg in PUBLIC_POINTS.items():
        if quark == "charm":
            continue
        path = OUTDIR / cfg["filename"]
        if not path.exists():
            raise FileNotFoundError(path)
        local = find_point(path, cfg["pt"], cfg["y"])
        public = cfg["public_pb_per_GeV"]
        rel = (local - public) / public
        status = "OK" if abs(rel) <= args.tolerance else "FAIL"
        print(
            f"{status} {quark}: local={local:.6e} public={public:.6e} "
            f"rel_diff={rel:+.3%} ({cfg['public_fragmentation']})"
        )
        failed = failed or status == "FAIL"

    charm = fit_charm_d0_feeddown_weight(fonll_env())
    max_abs_rel = float(charm["max_abs_relative_difference"])
    status = "OK" if max_abs_rel <= args.tolerance else "FAIL"
    print(
        f"{status} charm: calibrated public-D0 feeddown max_abs_rel_diff="
        f"{max_abs_rel:.3%} over {len(charm['points'])} CTEQ6.6 y=0 public points"
    )
    for point in charm["points"]:
        if math.isclose(point["pt"], PUBLIC_POINTS["charm"]["pt"], abs_tol=1e-12):
            print(
                f"  pT=5,y=0 local={point['local']:.6e} public={point['public']:.6e} "
                f"rel_diff={point['relative_difference']:+.3%} "
                f"({PUBLIC_POINTS['charm']['public_fragmentation']})"
            )
            break
    failed = failed or status == "FAIL"

    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
