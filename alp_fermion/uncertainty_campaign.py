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
            "cbs_amplitude_scale": scale,
            "template_variant": "central",
            "production_mode": "fresh_600k",
        })
    return output


def all_variations(grid_dir: Path) -> list[dict]:
    fonll = discover_fonll_variations(grid_dir)
    return [*fonll, *auxiliary_variations(Path(fonll[0]["grid_path"]))]


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
    """Combine exact full-run contours into the stable published band schema."""
    axes = {
        axis: list(raw.loc[raw["axis"] == axis, "variation"].unique())
        for axis in ("central", "scale", "pdf", "mb", "decay_gg", "cbs")
    }
    expected_counts = {
        "central": 1, "scale": 6, "pdf": 100, "mb": 2,
        "decay_gg": 3, "cbs": 2,
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
                raw.loc[raw["mass_GeV"] == mass, "has_sensitivity"].any()
            ),
        }
        for boundary in BOUNDARIES:
            xc, central_open, _ = _log_boundary(
                raw, mass, "central", boundary
            )
            record[f"{boundary}_central"] = float(central.get(boundary, np.nan))
            record[f"{boundary}_open"] = central_open
            if xc is None:
                record[f"{boundary}_band_lo"] = np.nan
                record[f"{boundary}_band_hi"] = np.nan
                record[f"{boundary}_variation_missing"] = False
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
                        missing = missing or not sensitive
                    else:
                        result.append(value)
                return result

            scale = [xc, *values(axes["scale"])]
            scale_up = max(scale) - xc
            scale_dn = xc - min(scale)
            pdf = values(axes["pdf"])
            pdf_sigma = float(np.std(pdf, ddof=1)) if len(pdf) >= 2 else 0.0
            mb_values = values(axes["mb"])
            mb_dev = max((abs(value - xc) for value in mb_values), default=0.0)

            def envelope(axis):
                candidates = [xc, *values(axes[axis])]
                return max(candidates) - xc, xc - min(candidates)

            gg_up, gg_dn = envelope("decay_gg")
            cbs_up, cbs_dn = envelope("cbs")
            total_up = float(np.sqrt(
                scale_up**2 + pdf_sigma**2 + mb_dev**2
                + gg_up**2 + cbs_up**2
            ))
            total_dn = float(np.sqrt(
                scale_dn**2 + pdf_sigma**2 + mb_dev**2
                + gg_dn**2 + cbs_dn**2
            ))
            record.update({
                f"{boundary}_band_lo": 10.0 ** (xc - total_dn),
                f"{boundary}_band_hi": 10.0 ** (xc + total_up),
                f"{boundary}_open": any_open,
                f"{boundary}_variation_missing": missing,
                f"{boundary}_scale_up_dex": scale_up,
                f"{boundary}_scale_dn_dex": scale_dn,
                f"{boundary}_pdf_sigma_dex": pdf_sigma,
                f"{boundary}_mb_dev_dex": mb_dev,
                f"{boundary}_gg_up_dex": gg_up,
                f"{boundary}_gg_dn_dex": gg_dn,
                f"{boundary}_cbs_up_dex": cbs_up,
                f"{boundary}_cbs_dn_dex": cbs_dn,
                f"{boundary}_total_up_dex": total_up,
                f"{boundary}_total_dn_dex": total_dn,
            })
        rows.append(record)
    return pd.DataFrame(rows).sort_values("mass_GeV")
