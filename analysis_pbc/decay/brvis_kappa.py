from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

from geometry.per_parent_efficiency import GeometryConfig, geometry_tag, is_default_geometry_config, normalize_geometry_config

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_KAPPA_TABLE = REPO_ROOT / "output" / "csv" / "analysis" / "decay_kappa_table.csv"
REQUIRED_P_MIN_GEV = 0.6
REQUIRED_SEPARATION_MM = 1.0
REQUIRED_COLUMNS = (
    "flavour",
    "mass_GeV",
    "kappa",
    "p_min_GeV",
    "separation_mm",
    "source_policy",
    "status",
)
MASS_MATCH_TOL = 1.0e-6
CUT_MATCH_TOL = 1.0e-6
OPTIONAL_GEOMETRY_COLUMNS = (
    "geometry_tag",
    "geometry_model",
    "tube_radius_m",
    "detector_thickness_m",
    "profile_inset_floor",
)


class KappaTableError(ValueError):
    """Raised when the calibrated kappa table is missing or inconsistent."""


def resolve_kappa_table_path(table_path: str | Path | None = None) -> Path:
    path = DEFAULT_KAPPA_TABLE if table_path is None else Path(table_path)
    return path.expanduser().resolve()


@lru_cache(maxsize=16)
def _load_kappa_table_cached(path_str: str) -> pd.DataFrame:
    path = Path(path_str)
    if not path.exists():
        raise FileNotFoundError(
            f"Calibrated kappa table not found: {path}. "
            "Run tools/decay/calibrate_brvis_kappa.py first."
        )

    df = pd.read_csv(path)
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise KappaTableError(
            f"Kappa table {path} is missing required columns: {missing}. "
            f"Required: {list(REQUIRED_COLUMNS)}"
        )

    if len(df) == 0:
        raise KappaTableError(f"Kappa table {path} is empty.")

    df = df.copy()
    df["flavour"] = df["flavour"].astype(str).str.strip().str.lower()
    df["status"] = df["status"].astype(str).str.strip().str.lower()
    df["source_policy"] = df["source_policy"].astype(str).str.strip()

    for col in ("mass_GeV", "kappa", "p_min_GeV", "separation_mm"):
        df[col] = pd.to_numeric(df[col], errors="coerce")

    for col in ("tube_radius_m", "detector_thickness_m"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    if "profile_inset_floor" in df.columns:
        val = df["profile_inset_floor"].astype(str).str.strip().str.lower()
        mapper = {
            "true": True,
            "1": True,
            "yes": True,
            "false": False,
            "0": False,
            "no": False,
        }
        parsed = val.map(mapper)
        bad_mask = parsed.isna()
        if bad_mask.any():
            bad_idx = int(df.index[bad_mask][0])
            raise KappaTableError(
                f"Kappa table {path} has non-boolean profile_inset_floor at row {bad_idx}."
            )
        df["profile_inset_floor"] = parsed.astype(bool)

    bad = df[
        (~np.isfinite(df["mass_GeV"]))
        | (~np.isfinite(df["kappa"]))
        | (~np.isfinite(df["p_min_GeV"]))
        | (~np.isfinite(df["separation_mm"]))
    ]
    if len(bad) > 0:
        raise KappaTableError(
            f"Kappa table {path} has non-finite numeric rows; first bad index={int(bad.index[0])}."
        )

    ok_df = df[df["status"] == "ok"].copy()
    if len(ok_df) == 0:
        raise KappaTableError(f"Kappa table {path} has no rows with status=ok.")

    if (ok_df["kappa"] <= 0.0).any():
        bad_idx = int(ok_df.index[(ok_df["kappa"] <= 0.0)][0])
        raise KappaTableError(f"Kappa table {path} has non-positive kappa at row {bad_idx}.")

    ok_df = ok_df.sort_values(["flavour", "mass_GeV"]).reset_index(drop=True)
    return ok_df


def _require_cut_consistency(df: pd.DataFrame, p_min_GeV: float, separation_mm: float, path: Path) -> None:
    p_vals = np.unique(df["p_min_GeV"].to_numpy(dtype=float))
    s_vals = np.unique(df["separation_mm"].to_numpy(dtype=float))

    if len(p_vals) != 1:
        raise KappaTableError(
            f"Kappa table {path} mixes multiple p_min_GeV values: {p_vals.tolist()}"
        )
    if len(s_vals) != 1:
        raise KappaTableError(
            f"Kappa table {path} mixes multiple separation_mm values: {s_vals.tolist()}"
        )

    p_table = float(p_vals[0])
    s_table = float(s_vals[0])
    if abs(p_table - float(p_min_GeV)) > CUT_MATCH_TOL:
        raise KappaTableError(
            f"Kappa table p_min mismatch: table={p_table:.6g} GeV, runtime={float(p_min_GeV):.6g} GeV"
        )
    if abs(s_table - float(separation_mm)) > CUT_MATCH_TOL:
        raise KappaTableError(
            f"Kappa table separation mismatch: table={s_table:.6g} mm, runtime={float(separation_mm):.6g} mm"
        )


def _require_geometry_consistency(
    df: pd.DataFrame,
    geometry_config: GeometryConfig | None,
    path: Path,
) -> None:
    cfg = normalize_geometry_config(geometry_config)
    runtime_tag = geometry_tag(cfg)

    has_any_geometry_col = any(col in df.columns for col in OPTIONAL_GEOMETRY_COLUMNS)
    if not has_any_geometry_col:
        if not is_default_geometry_config(cfg):
            raise KappaTableError(
                f"Kappa table {path} has no geometry metadata and is treated as legacy default-only. "
                f"Runtime geometry '{runtime_tag}' is non-default; recalibrate kappa for this geometry."
            )
        return

    if "geometry_tag" not in df.columns:
        raise KappaTableError(
            f"Kappa table {path} has partial geometry metadata but is missing 'geometry_tag'. "
            "Recalibrate to produce a consistent table."
        )

    tags = np.unique(df["geometry_tag"].astype(str).str.strip().to_numpy())
    if len(tags) != 1:
        raise KappaTableError(
            f"Kappa table {path} mixes multiple geometry_tag values: {tags.tolist()}"
        )
    table_tag = str(tags[0])
    if table_tag != runtime_tag:
        raise KappaTableError(
            f"Kappa table geometry_tag mismatch: table={table_tag}, runtime={runtime_tag}. "
            "Recalibrate kappa table for the requested geometry."
        )

    if "geometry_model" in df.columns:
        vals = np.unique(df["geometry_model"].astype(str).str.strip().str.lower().to_numpy())
        if len(vals) != 1:
            raise KappaTableError(
                f"Kappa table {path} mixes multiple geometry_model values: {vals.tolist()}"
            )
        if str(vals[0]) != cfg.model:
            raise KappaTableError(
                f"Kappa table geometry_model mismatch: table={vals[0]}, runtime={cfg.model}."
            )

    if "tube_radius_m" in df.columns:
        vals = np.unique(df["tube_radius_m"].to_numpy(dtype=float))
        if len(vals) != 1:
            raise KappaTableError(
                f"Kappa table {path} mixes multiple tube_radius_m values: {vals.tolist()}"
            )
        if not np.isfinite(vals[0]) or abs(float(vals[0]) - float(cfg.tube_radius_m)) > CUT_MATCH_TOL:
            raise KappaTableError(
                f"Kappa table tube_radius_m mismatch: table={float(vals[0]):.6g} m, "
                f"runtime={float(cfg.tube_radius_m):.6g} m"
            )

    if "detector_thickness_m" in df.columns:
        vals = np.unique(df["detector_thickness_m"].to_numpy(dtype=float))
        if len(vals) != 1:
            raise KappaTableError(
                f"Kappa table {path} mixes multiple detector_thickness_m values: {vals.tolist()}"
            )
        if not np.isfinite(vals[0]) or abs(float(vals[0]) - float(cfg.detector_thickness_m)) > CUT_MATCH_TOL:
            raise KappaTableError(
                f"Kappa table detector_thickness_m mismatch: table={float(vals[0]):.6g} m, "
                f"runtime={float(cfg.detector_thickness_m):.6g} m"
            )

    if "profile_inset_floor" in df.columns:
        vals = np.unique(df["profile_inset_floor"].astype(bool).to_numpy())
        if len(vals) != 1:
            raise KappaTableError(
                f"Kappa table {path} mixes multiple profile_inset_floor values: {vals.tolist()}"
            )
        if bool(vals[0]) != bool(cfg.profile_inset_floor):
            raise KappaTableError(
                f"Kappa table profile_inset_floor mismatch: table={bool(vals[0])}, "
                f"runtime={bool(cfg.profile_inset_floor)}"
            )


def lookup_kappa(
    flavour: str,
    mass_GeV: float,
    p_min_GeV: float,
    separation_mm: float,
    table_path: str | Path | None = None,
    geometry_config: GeometryConfig | None = None,
) -> float:
    flavour = str(flavour).strip().lower()
    mass = float(mass_GeV)
    path = resolve_kappa_table_path(table_path)
    df = _load_kappa_table_cached(str(path))

    _require_cut_consistency(df, float(p_min_GeV), float(separation_mm), path)
    _require_geometry_consistency(df, geometry_config, path)

    df_flavour = df[df["flavour"] == flavour].copy()
    if len(df_flavour) == 0:
        raise KappaTableError(f"No kappa rows for flavour='{flavour}' in table {path}.")

    masses = df_flavour["mass_GeV"].to_numpy(dtype=float)
    kappas = df_flavour["kappa"].to_numpy(dtype=float)

    exact = np.where(np.abs(masses - mass) <= MASS_MATCH_TOL)[0]
    if len(exact) > 0:
        return float(kappas[int(exact[0])])

    order = np.argsort(masses)
    masses = masses[order]
    kappas = kappas[order]

    if mass < float(masses[0]) or mass > float(masses[-1]):
        raise KappaTableError(
            f"Mass {mass:.6g} GeV out of calibrated range for flavour='{flavour}': "
            f"[{masses[0]:.6g}, {masses[-1]:.6g}] GeV"
        )

    hi = int(np.searchsorted(masses, mass, side="right"))
    lo = hi - 1
    if lo < 0 or hi >= len(masses):
        raise KappaTableError(
            f"Could not bracket mass {mass:.6g} GeV for flavour='{flavour}' in table {path}."
        )

    x0, x1 = float(masses[lo]), float(masses[hi])
    y0, y1 = float(kappas[lo]), float(kappas[hi])
    if abs(x1 - x0) <= MASS_MATCH_TOL:
        return float(y0)

    t = (mass - x0) / (x1 - x0)
    interp = y0 + t * (y1 - y0)
    if not np.isfinite(interp) or interp <= 0.0:
        raise KappaTableError(
            f"Interpolated non-positive/invalid kappa for flavour='{flavour}' mass={mass:.6g}: {interp}"
        )
    return float(interp)
