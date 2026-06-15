#!/usr/bin/env python3
"""Hard-abort validation gate for the variation campaign.

Generates the central (member 0, scale 1,1, central mass) grid with the
refactored generator and asserts its dsigma column reproduces the committed
reference grid within tolerance. This proves the PDF-member/scale/mass
plumbing did not perturb the central result before the full campaign commits
days of compute. Also reports a grid-coverage tail bound (fraction of the
cross section in the outermost pT and rapidity bins).

Exit status is nonzero if any quark's reproduction exceeds --rtol.
"""
from __future__ import annotations

import argparse
import importlib.util
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
REFERENCE_DIR_DEFAULT = ROOT.parent / "fonll-nnpdf40" / "grids"


def _load_generator_module():
    path = Path(__file__).resolve().parent / "generate_meson_grids.py"
    spec = importlib.util.spec_from_file_location("generate_meson_grids", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


GMG = _load_generator_module()


def data_column(path: Path) -> np.ndarray:
    data = np.loadtxt(path, comments="#")
    return data[:, 0], data[:, 1], data[:, 2]


def coverage_tail_bound(pt: np.ndarray, y: np.ndarray, ds: np.ndarray) -> dict:
    n_pt, n_y = len(GMG.PT_VALUES), len(GMG.Y_VALUES)
    grid = ds.reshape(n_pt, n_y)
    total = float(np.trapz(np.trapz(grid, GMG.Y_VALUES, axis=1), GMG.PT_VALUES))
    # Fraction of |sigma| in the outermost pT row (near 50 GeV) and y columns (|y|~3).
    pt_tail = float(np.trapz(grid[-1, :], GMG.Y_VALUES)) * (GMG.PT_VALUES[-1] - GMG.PT_VALUES[-2])
    y_tail = (
        float(np.trapz(grid[:, 0], GMG.PT_VALUES)) + float(np.trapz(grid[:, -1], GMG.PT_VALUES))
    ) * (GMG.Y_VALUES[1] - GMG.Y_VALUES[0])
    return {
        "total_pb": total,
        "pt_edge_fraction": pt_tail / total if total else float("nan"),
        "y_edge_fraction": y_tail / total if total else float("nan"),
    }


def gate_quark(quark: str, out_root: Path, reference_dir: Path, rtol: float, grid_workers: int = 1) -> bool:
    feeddown = None
    if quark == "charm":
        feeddown = GMG.load_or_fit_charm_feeddown_weight(
            GMG.fonll_env(), out_root / "charm_feeddown_calibration.json"
        )
    central = GMG.Variation("central", "central", 1.0, 1.0, 0, float(GMG.QUARKS[quark]["mass"]))
    summary = GMG.generate(
        "nlo", quark, central, reuse_existing_grids=False, grid_workers=grid_workers,
        feeddown_calibration=feeddown, out_root=out_root, run_root=out_root / "run",
    )
    fresh_path = Path(summary["path"])
    ref_path = reference_dir / (
        f"fonll_pp14tev_nnpdf40_nlo_as_01180_fonll_meson_dsdpTdy_"
        f"pt0-50_y-3to3_central_{quark}.dat"
    )
    fpt, fy, fds = data_column(fresh_path)
    rpt, ry, rds = data_column(ref_path)
    if fds.shape != rds.shape:
        print(f"[{quark}] FAIL: shape {fds.shape} vs reference {rds.shape}")
        return False
    nz = rds != 0
    rel = np.zeros_like(rds)
    rel[nz] = np.abs(fds[nz] - rds[nz]) / np.abs(rds[nz])
    max_rel = float(rel.max())
    cov = coverage_tail_bound(fpt, fy, fds)
    status = "PASS" if max_rel <= rtol else "FAIL"
    print(f"[{quark}] {status}: max relative diff vs committed central = {max_rel:.3e} (rtol={rtol:.1e})")
    print(f"[{quark}] coverage: total={cov['total_pb']:.6e} pb; "
          f"pT-edge fraction={cov['pt_edge_fraction']:.3e}; "
          f"y-edge fraction={cov['y_edge_fraction']:.3e}")
    return max_rel <= rtol


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quark", action="append", choices=sorted(GMG.QUARKS),
                        help="quark to check; may be repeated. Default: bottom and charm.")
    parser.add_argument("--reference-dir", default=str(REFERENCE_DIR_DEFAULT),
                        help="directory holding the committed central reference grids")
    parser.add_argument("--rtol", type=float, default=1e-3,
                        help="max allowed relative difference vs the committed central grid")
    parser.add_argument("--out-root", default="",
                        help="output dir for fresh grids (default: a temp dir)")
    parser.add_argument("--grid-workers", type=int, default=1,
                        help="parallel rapidity chunks for the gate grid (speeds up the one-off check)")
    args = parser.parse_args()
    quarks = args.quark or ["bottom", "charm"]
    out_root = Path(args.out_root) if args.out_root else Path(tempfile.mkdtemp(prefix="vargate_"))
    out_root.mkdir(parents=True, exist_ok=True)
    print(f"gate output dir: {out_root}")

    ok = True
    for quark in quarks:
        ok = gate_quark(quark, out_root, Path(args.reference_dir), args.rtol, args.grid_workers) and ok
    if not ok:
        print("GATE FAILED -- do not launch the full campaign until resolved")
        sys.exit(1)
    print("GATE PASSED -- central reproduction confirmed")


if __name__ == "__main__":
    main()
