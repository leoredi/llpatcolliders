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
    structural_alternative_table,
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
            "storage_state": marker.get("storage_state", "full"),
            "pre_compaction_hashes": {
                "vector_tree_sha256": json.loads(
                    Path(marker["production_marker"]).read_text()
                ).get("vector_tree_sha256"),
                "geometry_tree_sha256": marker["geometry_tree_sha256"],
            },
        })
    return pd.concat(rows, ignore_index=True), registry


def _headline(band: pd.DataFrame) -> dict:
    sensitive = band[band["has_sensitivity"]].copy()
    output = {
        "n_mass_points": len(band),
        "n_central_sensitive": int(band["has_sensitivity"].sum()),
        "n_any_variation_sensitive": int(band["any_variation_sensitive"].sum()),
        "n_decay_structure_sensitive": int(
            band["decay_structure_has_sensitivity"].sum()
        ),
        "n_decay_structure_restored": int(
            band["decay_structure_restores_sensitivity"].sum()
        ),
        "n_decay_structure_removed": int(
            band["decay_structure_removes_sensitivity"].sum()
        ),
        "n_decay_structure_topology_differences": int(
            band["decay_structure_topology_differs"].sum()
        ),
    }
    for boundary in ("invf_min", "invf_max"):
        central = sensitive[f"{boundary}_central"]
        for direction, edge in (("up", "hi"), ("dn", "lo")):
            envelope = sensitive[f"{boundary}_envelope_{edge}"]
            if direction == "up":
                values = np.log10(envelope / central)
            else:
                values = np.log10(central / envelope)
            values = values.replace([np.inf, -np.inf], np.nan)
            if values.notna().any():
                index = values.idxmax()
                key = f"max_{boundary}_envelope_{direction}_dex"
                output[key] = float(values.loc[index])
                output[f"{key}_mass_GeV"] = float(band.loc[index, "mass_GeV"])
        repeat = sensitive[f"{boundary}_repeat_max_abs_dex"].replace(
            [np.inf, -np.inf], np.nan
        )
        flagged = sensitive[sensitive[f"{boundary}_repeat_not_subdominant"]]
        output[f"{boundary}_numerical_control"] = {
            "median_repeat_max_abs_dex": (
                float(repeat.median()) if repeat.notna().any() else None
            ),
            "max_repeat_abs_dex": (
                float(repeat.max()) if repeat.notna().any() else None
            ),
            "max_repeat_abs_dex_mass_GeV": (
                float(sensitive.loc[repeat.idxmax(), "mass_GeV"])
                if repeat.notna().any() else None
            ),
            "n_not_subdominant": len(flagged),
            "not_subdominant_masses_GeV": [
                float(value) for value in flagged["mass_GeV"]
            ],
        }
    return output


def _campaign_code_states(registry: list[dict]) -> tuple[dict, dict]:
    """Require consistency within the halo campaign and structural add-on."""
    halo = {
        json.dumps(item["code"], sort_keys=True)
        for item in registry
        if item["variation"]["axis"] != "decay_structure"
    }
    structural = {
        json.dumps(item["code"], sort_keys=True)
        for item in registry
        if item["variation"]["axis"] == "decay_structure"
    }
    if len(halo) != 1:
        raise ValueError(
            "pointwise-halo/control variations were produced from different code states"
        )
    if len(structural) != 1:
        raise ValueError(
            "decay-structure variations were produced from different code states"
        )
    return json.loads(next(iter(halo))), json.loads(next(iter(structural)))


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
        default=(
            Path(__file__).resolve().parent / "data" / "published" / "bundle"
        ),
    )
    args = parser.parse_args(argv)
    scratch = Path(args.scratch_root).expanduser().resolve()
    variations = all_variations(args.grid_dir)
    raw, registry = load_completed(scratch, variations)
    band = combine_band(raw)
    structural = structural_alternative_table(raw)
    out_dir = Path(args.out_dir)
    raw_path = out_dir / "bc10_uncertainty_variations.csv"
    band_path = out_dir / "bc10_single_source_variation_envelope.csv"
    structural_path = out_dir / "bc10_decay_2310_structural_alternative.csv"
    manifest_path = out_dir / "UNCERTAINTY_MANIFEST.json"
    _atomic_csv(raw, raw_path)
    _atomic_csv(band, band_path)
    _atomic_csv(structural, structural_path)
    halo_code, structural_code = _campaign_code_states(registry)
    atomic_json(manifest_path, {
        "artifact": "GRENDEL BC10 exact single-source variation envelope",
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
            "decay_structure": (
                "one exact SensCalc arXiv:2310.03524 width, exclusive-BR, and "
                "three-body-matrix-element model with independent 20000-event "
                "Pythia templates and full geometry/reconstruction; exact central "
                "production-vector reuse; published as a dashed one-sided structural "
                "comparison and excluded from the pointwise halo"
            ),
            "numerical_control": (
                "two same-physics central repeats with fresh independent 600000-"
                "parent production pools and distinct reconstruction RNG offsets; "
                "excluded from the physical envelope and compared pointwise with it"
            ),
            "combination": (
                "pointwise one-source-at-a-time intervals in log10(1/f): named "
                "scale/mb/gg/C_bs extrema and NNPDF replica 16th/84th "
                "percentiles; display envelope is their outermost boundary; "
                "the 2310 structural contour and numerical repeats are excluded; "
                "no quadrature combination and no confidence-interval claim"
            ),
            "label": "single_source_variation_envelope",
        },
        "variation_counts": {
            "central": 1, "scale": 6, "pdf": 100, "mb": 2,
            "decay_gg": 3, "cbs": 2, "decay_structure": 1,
            "numerical_control": 2,
            "pointwise_halo_total": 114,
            "physics_plus_structural_total": 115,
            "total_with_controls": 117,
        },
        "outputs": {
            raw_path.name: sha256_file(raw_path),
            band_path.name: sha256_file(band_path),
            structural_path.name: sha256_file(structural_path),
        },
        "headline": _headline(band),
        "code": halo_code,
        "decay_structure_code": structural_code,
        "registry": registry,
    })
    print(f"wrote {raw_path} ({len(raw)} rows)")
    print(f"wrote {band_path} ({len(band)} rows)")
    print(f"wrote {structural_path} ({len(structural)} rows)")
    print(f"wrote {manifest_path} ({len(registry)} exact variations)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
