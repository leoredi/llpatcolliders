#!/usr/bin/env python3
"""Validate and combine the completed exact BC10 uncertainty campaign."""

from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path

import numpy as np
import pandas as pd

from alp_fermion.uncertainty_campaign import (
    all_variations,
    atomic_json,
    combine_band,
    sha256_file,
)


def _atomic_csv(frame: pd.DataFrame, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(temporary, index=False)
    temporary.replace(path)


def load_completed(scratch: Path, variations) -> tuple[pd.DataFrame, list[dict]]:
    rows = []
    registry = []
    for variation in variations:
        marker_path = scratch / "runs" / variation["name"] / "variation.complete.json"
        if not marker_path.exists():
            raise FileNotFoundError(
                f"variation {variation['name']} is incomplete: {marker_path}"
            )
        marker = json.loads(marker_path.read_text())
        if marker["variation"] != variation:
            raise ValueError(f"variation definition mismatch in {marker_path}")
        sensitivity = Path(marker["sensitivity_csv"])
        actual_sha = sha256_file(sensitivity)
        if actual_sha != marker["sensitivity_csv_sha256"]:
            raise ValueError(f"sensitivity checksum mismatch: {sensitivity}")
        frame = pd.read_csv(sensitivity)
        if len(frame) != 99:
            raise ValueError(f"{sensitivity} has {len(frame)} rows, expected 99")
        frame.insert(0, "axis", variation["axis"])
        frame.insert(0, "variation", variation["name"])
        rows.append(frame)
        registry.append({
            "variation": variation,
            "completion_marker": str(marker_path),
            "completion_marker_sha256": sha256_file(marker_path),
            "sensitivity_csv": str(sensitivity),
            "sensitivity_csv_sha256": actual_sha,
            "production_marker": marker["production_marker"],
            "production_marker_sha256": marker["production_marker_sha256"],
            "template": marker["template"],
            "code": marker["code"],
        })
    return pd.concat(rows, ignore_index=True), registry


def _headline(band: pd.DataFrame) -> dict:
    sensitive = band[band["has_sensitivity"]].copy()
    output = {
        "n_mass_points": len(band),
        "n_central_sensitive": int(band["has_sensitivity"].sum()),
        "n_any_variation_sensitive": int(band["any_variation_sensitive"].sum()),
    }
    for boundary in ("invf_min", "invf_max"):
        for direction in ("up", "dn"):
            column = f"{boundary}_total_{direction}_dex"
            values = sensitive[column].replace([np.inf, -np.inf], np.nan)
            if values.notna().any():
                index = values.idxmax()
                output[f"max_{column}"] = float(values.loc[index])
                output[f"max_{column}_mass_GeV"] = float(band.loc[index, "mass_GeV"])
    return output


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--scratch-root",
        default=os.environ.get("ALP_UNCERTAINTY_SCRATCH"),
        required=os.environ.get("ALP_UNCERTAINTY_SCRATCH") is None,
    )
    parser.add_argument("--grid-dir", type=Path, required=True)
    parser.add_argument(
        "--out-dir", type=Path,
        default=Path(__file__).resolve().parent / "data" / "published",
    )
    args = parser.parse_args(argv)
    scratch = Path(args.scratch_root).expanduser().resolve()
    variations = all_variations(args.grid_dir)
    raw, registry = load_completed(scratch, variations)
    band = combine_band(raw)
    out_dir = Path(args.out_dir)
    raw_path = out_dir / "bc10_uncertainty_variations.csv"
    band_path = out_dir / "bc10_uncertainty_band.csv"
    manifest_path = out_dir / "UNCERTAINTY_MANIFEST.json"
    _atomic_csv(raw, raw_path)
    _atomic_csv(band, band_path)
    code_states = {json.dumps(item["code"], sort_keys=True) for item in registry}
    if len(code_states) != 1:
        raise ValueError("campaign variations were produced from different code states")
    atomic_json(manifest_path, {
        "artifact": "GRENDEL BC10 exact theory-uncertainty band",
        "generated_unix": time.time(),
        "method": {
            "production": (
                "independent 600000-parent production, geometry, reconstruction, "
                "and sensitivity run for central plus 6 scale, 100 NNPDF replica, "
                "and 2 bottom-mass FONLL grids; no importance reweighting"
            ),
            "decay_gg": (
                "independently simulated full-branching Pythia u/d/s surrogate "
                "templates, each propagated through its own geometry/reconstruction run"
            ),
            "cbs": (
                "+/-20% C_bs amplitude, implemented as fresh production rates "
                "scaled by 0.8^2 and 1.2^2 and full downstream runs"
            ),
            "combination": (
                "log10(1/f): asymmetric scale/gg/C_bs envelopes, NNPDF replica "
                "sample standard deviation, max absolute mb displacement, sources "
                "combined in quadrature separately up/down"
            ),
        },
        "variation_counts": {
            "central": 1, "scale": 6, "pdf": 100, "mb": 2,
            "decay_gg": 3, "cbs": 2, "total": 114,
        },
        "outputs": {
            raw_path.name: sha256_file(raw_path),
            band_path.name: sha256_file(band_path),
        },
        "headline": _headline(band),
        "code": json.loads(next(iter(code_states))),
        "registry": registry,
    })
    print(f"wrote {raw_path} ({len(raw)} rows)")
    print(f"wrote {band_path} ({len(band)} rows)")
    print(f"wrote {manifest_path} ({len(registry)} exact variations)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
