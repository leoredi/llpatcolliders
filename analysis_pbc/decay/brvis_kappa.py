from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

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


def lookup_kappa(
    flavour: str,
    mass_GeV: float,
    p_min_GeV: float,
    separation_mm: float,
    table_path: str | Path | None = None,
) -> float:
    flavour = str(flavour).strip().lower()
    mass = float(mass_GeV)
    path = resolve_kappa_table_path(table_path)
    df = _load_kappa_table_cached(str(path))

    _require_cut_consistency(df, float(p_min_GeV), float(separation_mm), path)

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
