"""Export the SensCalc BC1 baseline exclusion envelope as a portable NPZ.

The source JSON stores mass in GeV and kinetic mixing epsilon.  Conversion to
epsilon squared is deliberately left to the plotter so the vendored data stay
faithful to their source convention.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


_ROOT = Path(__file__).resolve().parents[1]
_OUTPUT = _ROOT / "data" / "constraints" / "existing_exclusion_senscalc.npz"
SOURCE_REVISION = "0bca050633aae16e148d47f21840fa07ff4b8724"


def _default_source() -> Path:
    rel = Path("shared/vendored/SensCalc-v1.3.3/contours/DP/Constraints-DP.json")
    for parent in Path(__file__).resolve().parents:
        candidate = parent / rel
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"workspace SensCalc source not found: {rel}")


def export(source: Path, output: Path = _OUTPUT) -> Path:
    baseline = json.loads(Path(source).read_text())["Baseline"]
    polygons = [np.asarray(points, dtype=float) for points in baseline]
    arrays = {f"polygon_{i}": points
              for i, points in enumerate(polygons)}
    for i, points in enumerate(polygons):
        arrays[f"closed_{i}"] = np.bool_(np.allclose(points[0], points[-1]))
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output,
        **arrays,
        n_polygons=np.int64(len(baseline)),
        source=np.str_("SensCalc v1.3.3 contours/DP/Constraints-DP.json:Baseline"),
        source_revision=np.str_(SOURCE_REVISION),
        mass_units=np.str_("GeV"),
        ordinate=np.str_("epsilon"),
    )
    return output


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--source", type=Path, default=None)
    parser.add_argument("--output", type=Path, default=_OUTPUT)
    args = parser.parse_args(argv)
    output = export(args.source or _default_source(), args.output)
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
