"""Definitions and provenance helpers for the exact BC10 uncertainty campaign."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


GRID_STEM = (
    "fonll_pp14tev_nnpdf40_nlo_as_01180_fonll_meson_"
    "dsdpTdy_pt0-50_y-3to3"
)
CBS_SCHEME_RELATIVE_AMPLITUDE = 0.20
BOUNDARIES = ("invf_min", "invf_max")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_tree(paths) -> str:
    """Hash file names and contents in stable order."""
    digest = hashlib.sha256()
    for path in sorted((Path(item) for item in paths), key=lambda item: str(item)):
        digest.update(path.name.encode())
        digest.update(b"\0")
        digest.update(sha256_file(path).encode())
        digest.update(b"\n")
    return digest.hexdigest()


def atomic_json(path: Path, payload) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def stable_seed(name: str) -> int:
    """Stable independent NumPy seed for one non-central production run."""
    digest = hashlib.sha256(f"GRENDEL-BC10-FONLL:{name}".encode()).digest()
    return int.from_bytes(digest[:4], "big") or 1


def discover_fonll_variations(grid_dir: Path) -> list[dict]:
    """Select the complete bottom-grid ensemble from its pinned manifest."""
    grid_dir = Path(grid_dir).expanduser().resolve()
    manifest_path = grid_dir / "variation_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    variations = []
    for entry in manifest["grids"]:
        if entry["quark"] != "bottom":
            continue
        tag = entry["variation_tag"]
        kind = entry["variation_kind"]
        if kind == "pdf" and int(entry["lhapdf_member"]) == 0:
            continue
        if kind not in {"central", "scale", "pdf", "mass"}:
            continue
        path = grid_dir / f"{GRID_STEM}_{tag}_bottom.dat"
        if not path.exists():
            raise FileNotFoundError(path)
        actual_sha = sha256_file(path)
        if actual_sha != entry["sha256"]:
            raise ValueError(f"FONLL checksum mismatch for {path}")
        axis = "mb" if kind == "mass" else kind
        variations.append({
            "name": tag,
            "axis": axis,
            "campaign_axis": "fonll",
            "grid_path": str(path),
            "grid_sha256": actual_sha,
            "muR": entry.get("muR"),
            "muF": entry.get("muF"),
            "pdf_member": entry.get("lhapdf_member"),
            "mb_GeV": entry.get("heavy_quark_mass_GeV"),
            "production_seed": 42 if kind == "central" else stable_seed(tag),
            "reco_seed_offset": 0,
            "cbs_amplitude_scale": 1.0,
            "template_variant": "central",
            "production_mode": "fresh_600k",
        })
    order = {"central": 0, "scale": 1, "pdf": 2, "mb": 3}
    variations.sort(key=lambda item: (order[item["axis"]], item["name"]))
    counts = {axis: sum(v["axis"] == axis for v in variations) for axis in order}
    expected = {"central": 1, "scale": 6, "pdf": 100, "mb": 2}
    if counts != expected:
        raise ValueError(f"incomplete FONLL bottom ensemble: {counts} != {expected}")
    return variations


def auxiliary_variations(central_grid: Path) -> list[dict]:
    grid_path = str(Path(central_grid).resolve())
    grid_sha = sha256_file(Path(central_grid))
    output = []
    for flavor in ("u", "d", "s"):
        output.append({
            "name": f"gg_{flavor}",
            "axis": "decay_gg",
            "campaign_axis": "decay_gg",
            "grid_path": grid_path,
            "grid_sha256": grid_sha,
            "production_seed": 42,
            "reco_seed_offset": 0,
            "cbs_amplitude_scale": 1.0,
            "template_variant": f"gg_{flavor}",
            "production_mode": "central_vectors_exact_reuse",
        })
    for direction, scale in (
        ("down", 1.0 - CBS_SCHEME_RELATIVE_AMPLITUDE),
        ("up", 1.0 + CBS_SCHEME_RELATIVE_AMPLITUDE),
    ):
        output.append({
            "name": f"cbs_{direction}",
            "axis": "cbs",
            "campaign_axis": "cbs",
            "grid_path": grid_path,
            "grid_sha256": grid_sha,
            "production_seed": 42,
            "reco_seed_offset": 0,
            "cbs_amplitude_scale": scale,
            "template_variant": "central",
            "production_mode": "fresh_600k",
        })
    return output


def numerical_control_variations(central_grid: Path) -> list[dict]:
    """Same-physics central repeats used only to measure numerical spread."""
    grid_path = str(Path(central_grid).resolve())
    grid_sha = sha256_file(Path(central_grid))
    output = []
    for index in (1, 2):
        name = f"central_repeat_{index}"
        output.append({
            "name": name,
            "axis": "numerical_control",
            "campaign_axis": "numerical_control",
            "grid_path": grid_path,
            "grid_sha256": grid_sha,
            "production_seed": stable_seed(f"{name}:production"),
            "reco_seed_offset": stable_seed(f"{name}:reconstruction"),
            "cbs_amplitude_scale": 1.0,
            "template_variant": "central",
            "production_mode": "fresh_600k",
        })
    return output


def all_variations(grid_dir: Path) -> list[dict]:
    fonll = discover_fonll_variations(grid_dir)
    central_grid = Path(fonll[0]["grid_path"])
    return [
        *fonll,
        *auxiliary_variations(central_grid),
        *numerical_control_variations(central_grid),
    ]


def _log_boundary(raw, mass, variation, boundary):
    selected = raw[
        (raw["mass_GeV"] == mass) & (raw["variation"] == variation)
    ]
    if len(selected) != 1:
        return None, False, False
    row = selected.iloc[0]
    sensitive = bool(row["has_sensitivity"])
    opened = bool(row.get(f"{boundary}_open", False))
    value = float(row.get(boundary, np.nan))
    if not sensitive or not np.isfinite(value) or value <= 0.0:
        return None, opened, sensitive
    return float(np.log10(value)), opened, sensitive


def combine_band(raw: pd.DataFrame) -> pd.DataFrame:
    """Build a pointwise one-source-at-a-time variation envelope.

    Named alternatives use their extrema. The NNPDF ensemble uses its 16th and
    84th percentiles; its raw extrema are deliberately not used as a headline
    interval. No sources are combined in quadrature and this is not a
    confidence interval. Same-physics central repeats are excluded from the
    envelope and reported separately as numerical controls.
    """
    axes = {
        axis: list(raw.loc[raw["axis"] == axis, "variation"].unique())
        for axis in (
            "central", "scale", "pdf", "mb", "decay_gg", "cbs",
            "numerical_control",
        )
    }
    expected_counts = {
        "central": 1, "scale": 6, "pdf": 100, "mb": 2,
        "decay_gg": 3, "cbs": 2, "numerical_control": 2,
    }
    counts = {axis: len(names) for axis, names in axes.items()}
    if counts != expected_counts:
        raise ValueError(f"incomplete contour ensemble: {counts} != {expected_counts}")
    rows = []
    for mass in sorted(raw["mass_GeV"].unique()):
        central = raw[
            (raw["mass_GeV"] == mass) & (raw["variation"] == "central")
        ].iloc[0]
        record = {
            "mass_GeV": mass,
            "has_sensitivity": bool(central["has_sensitivity"]),
            "any_variation_sensitive": bool(
                raw.loc[
                    (raw["mass_GeV"] == mass)
                    & (raw["axis"] != "numerical_control"),
                    "has_sensitivity",
                ].any()
            ),
            "envelope_definition": "single_source_variation_envelope",
        }
        for boundary in BOUNDARIES:
            xc, central_open, _ = _log_boundary(
                raw, mass, "central", boundary
            )
            record[f"{boundary}_central"] = float(central.get(boundary, np.nan))
            record[f"{boundary}_open"] = central_open
            if xc is None:
                record[f"{boundary}_envelope_lo"] = np.nan
                record[f"{boundary}_envelope_hi"] = np.nan
                record[f"{boundary}_any_variation_open"] = central_open
                record[f"{boundary}_variation_missing"] = False
                record[f"{boundary}_envelope_lo_source"] = ""
                record[f"{boundary}_envelope_hi_source"] = ""
                record[f"{boundary}_repeat_missing"] = False
                record[f"{boundary}_repeat_not_subdominant"] = False
                continue

            missing = False
            any_open = central_open

            def values(names):
                nonlocal missing, any_open
                result = []
                for name in names:
                    value, opened, sensitive = _log_boundary(
                        raw, mass, name, boundary
                    )
                    any_open = any_open or opened
                    if value is None:
                        missing = True
                    else:
                        result.append((name, value))
                return result

            source_intervals = {}
            for axis, label in (
                ("scale", "scale"),
                ("mb", "mb"),
                ("decay_gg", "gg"),
                ("cbs", "cbs"),
            ):
                candidates = [("central", xc), *values(axes[axis])]
                lo_name, lo = min(candidates, key=lambda item: item[1])
                hi_name, hi = max(candidates, key=lambda item: item[1])
                source_intervals[label] = (lo, hi)
                record.update({
                    f"{boundary}_{label}_lo": 10.0 ** lo,
                    f"{boundary}_{label}_hi": 10.0 ** hi,
                    f"{boundary}_{label}_lo_variation": lo_name,
                    f"{boundary}_{label}_hi_variation": hi_name,
                })

            pdf_pairs = values(axes["pdf"])
            pdf = np.asarray([value for _, value in pdf_pairs], dtype=float)
            if len(pdf):
                pdf_p16, pdf_p84 = np.quantile(pdf, [0.16, 0.84])
                pdf_std = float(np.std(pdf, ddof=1)) if len(pdf) >= 2 else np.nan
            else:
                pdf_p16 = pdf_p84 = xc
                pdf_std = np.nan
            pdf_lo = min(xc, float(pdf_p16))
            pdf_hi = max(xc, float(pdf_p84))
            source_intervals["pdf"] = (pdf_lo, pdf_hi)
            record.update({
                f"{boundary}_pdf_p16": 10.0 ** float(pdf_p16),
                f"{boundary}_pdf_p84": 10.0 ** float(pdf_p84),
                f"{boundary}_pdf_log10_std": pdf_std,
                f"{boundary}_pdf_n_finite": len(pdf),
            })

            lo_source, lo = min(
                ((source, interval[0]) for source, interval in source_intervals.items()),
                key=lambda item: item[1],
            )
            hi_source, hi = max(
                ((source, interval[1]) for source, interval in source_intervals.items()),
                key=lambda item: item[1],
            )
            record.update({
                f"{boundary}_envelope_lo": 10.0 ** lo,
                f"{boundary}_envelope_hi": 10.0 ** hi,
                f"{boundary}_any_variation_open": any_open,
                f"{boundary}_variation_missing": missing,
                f"{boundary}_envelope_lo_source": lo_source,
                f"{boundary}_envelope_hi_source": hi_source,
            })

            repeat_values = []
            repeat_missing = False
            for name in axes["numerical_control"]:
                value, _, sensitive = _log_boundary(
                    raw, mass, name, boundary
                )
                if value is None:
                    repeat_missing = True
                    record[f"{boundary}_{name}"] = np.nan
                else:
                    repeat_values.append(value)
                    record[f"{boundary}_{name}"] = 10.0 ** value
            repeat_abs_dex = np.abs(np.asarray(repeat_values) - xc)
            repeat_abs_fraction = np.abs(10.0 ** (
                np.asarray(repeat_values) - xc
            ) - 1.0)
            physical_max_abs_dex = max(xc - lo, hi - xc)
            repeat_max_abs_dex = (
                float(np.max(repeat_abs_dex)) if len(repeat_abs_dex) else np.nan
            )
            record.update({
                f"{boundary}_physical_envelope_max_abs_dex": physical_max_abs_dex,
                f"{boundary}_repeat_median_abs_dex": (
                    float(np.median(repeat_abs_dex))
                    if len(repeat_abs_dex) else np.nan
                ),
                f"{boundary}_repeat_max_abs_dex": repeat_max_abs_dex,
                f"{boundary}_repeat_median_abs_fraction": (
                    float(np.median(repeat_abs_fraction))
                    if len(repeat_abs_fraction) else np.nan
                ),
                f"{boundary}_repeat_max_abs_fraction": (
                    float(np.max(repeat_abs_fraction))
                    if len(repeat_abs_fraction) else np.nan
                ),
                f"{boundary}_repeat_missing": repeat_missing,
                f"{boundary}_repeat_not_subdominant": bool(
                    repeat_missing
                    or (
                        np.isfinite(repeat_max_abs_dex)
                        and repeat_max_abs_dex > 0.0
                        and repeat_max_abs_dex >= physical_max_abs_dex
                    )
                ),
            })
        rows.append(record)
    return pd.DataFrame(rows).sort_values("mass_GeV")
