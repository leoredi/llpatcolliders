#!/usr/bin/env python3

from __future__ import annotations

import re
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

THIS_FILE = Path(__file__).resolve()
REPO_ROOT = THIS_FILE.parents[2]
ANALYSIS_ROOT = REPO_ROOT / "analysis_pbc"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(ANALYSIS_ROOT) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_ROOT))

from config_mass_grid import MASS_GRID
from geometry.per_parent_efficiency import (
    GeometryConfig,
    build_drainage_gallery_mesh,
    geometry_tag,
    is_default_geometry_config,
    normalize_geometry_config,
    preprocess_hnl_csv,
)
from limits import run as limits_run

MASS_FILTER_TOL = 5e-4
ALL_FLAVOURS = ("electron", "muon", "tau")


def threshold_for_mass(mass_GeV: float, switch_mass: float) -> float:
    """Validation threshold: tighter below the overlay switch mass."""
    return 0.10 if float(mass_GeV) < float(switch_mass) else 0.15
FLAVOUR_TO_BENCHMARK = {
    "electron": "100",
    "muon": "010",
    "tau": "001",
}


def parse_csv_list(raw: str) -> List[str]:
    return [x.strip() for x in raw.split(",") if x.strip()]


def parse_mass_list(raw: str) -> List[float]:
    return [float(x) for x in parse_csv_list(raw)]


def resolve_flavours(raw: str) -> List[str]:
    flavours = parse_csv_list(raw)
    if not flavours:
        return list(ALL_FLAVOURS)
    unknown = [f for f in flavours if f not in ALL_FLAVOURS]
    if unknown:
        raise ValueError(f"Unknown flavour(s): {unknown}. Allowed: {ALL_FLAVOURS}")
    return flavours


def resolve_masses(raw_masses: str | None, from_mass_grid: bool, min_mass: float, max_mass: float) -> List[float]:
    if raw_masses:
        masses = parse_mass_list(raw_masses)
    elif from_mass_grid:
        masses = list(MASS_GRID)
    else:
        raise ValueError("No masses selected. Provide --masses or keep --from-mass-grid enabled.")

    if max_mass < min_mass:
        raise ValueError("max_mass must be >= min_mass")

    masses = sorted({float(m) for m in masses if float(min_mass) <= float(m) <= float(max_mass)})
    return masses


def _mass_selected(mass_val: float, selected_masses: List[float]) -> bool:
    if not selected_masses:
        return True
    return any(abs(float(mass_val) - float(m)) <= MASS_FILTER_TOL for m in selected_masses)


def discover_sim_files_for_flavour(
    flavour: str,
    selected_masses: List[float] | None = None,
    max_mass: float | None = None,
    allow_variant_drop: bool = False,
) -> Dict[Tuple[float, str], List[Tuple[Path, str]]]:
    selected_masses = selected_masses or []

    pattern = re.compile(
        rf"^HNL_([0-9]+p[0-9]{{1,2}})GeV_{flavour}_"
        r"((?:kaon|charm|beauty|Bc|ew|all|combined)(?:_ff)?)"
        r"(?:_(direct|fromTau))?"
        r"(?:_(hardBc|hardccbar|hardbbbar)(?:_pTHat([0-9]+(?:p[0-9]+)?))?)?"
        r"\.csv$"
    )

    files = []
    for path in limits_run.SIM_DIR.glob(f"*{flavour}*.csv"):
        match = pattern.search(path.name)
        if not match:
            continue
        if path.stat().st_size < 1000:
            continue

        mass_str = match.group(1)
        mass_val = float(mass_str.replace("p", "."))
        if max_mass is not None and mass_val > float(max_mass):
            continue
        if not _mass_selected(mass_val, selected_masses):
            continue

        regime_token = match.group(2)
        mode = match.group(3)
        qcd_mode = match.group(4) or "auto"
        pthat_token = match.group(5)
        pthat_min = float(pthat_token.replace("p", ".")) if pthat_token else None

        is_ff = regime_token.endswith("_ff")
        base_regime = regime_token.replace("_ff", "")

        files.append((mass_val, mass_str, base_regime, mode, is_ff, qcd_mode, pthat_min, path))

    files_by_mass: Dict[Tuple[float, str], List[Tuple[str, str | None, bool, str, float | None, Path]]] = {}
    for mass_val, mass_str, base_regime, mode, is_ff, qcd_mode, pthat_min, path in files:
        key = (mass_val, mass_str)
        files_by_mass.setdefault(key, []).append((base_regime, mode, is_ff, qcd_mode, pthat_min, path))

    def _label(
        base_regime: str,
        mode: str | None,
        is_ff: bool,
        qcd_mode: str,
        pthat_min: float | None,
    ) -> str:
        label = base_regime
        if is_ff:
            label += "_ff"
        if mode:
            label += f"_{mode}"
        if qcd_mode != "auto":
            label += f"_{qcd_mode}"
            if pthat_min is not None:
                pthat_text = f"{pthat_min:g}".replace(".", "p")
                label += f"_pTHat{pthat_text}"
        return label

    def _variant_priority(
        base_regime: str,
        is_ff: bool,
        qcd_mode: str,
        pthat_min: float | None,
    ) -> tuple[int, int, float]:
        if base_regime == "charm" and qcd_mode == "hardccbar":
            qcd_priority = 3
        elif base_regime in {"beauty", "Bc"} and qcd_mode in {"hardbbbar", "hardBc"}:
            qcd_priority = 3
        elif qcd_mode != "auto":
            qcd_priority = 2
        else:
            qcd_priority = 1
        ff_priority = 1 if is_ff else 0
        pthat_priority = float(pthat_min) if pthat_min is not None else -1.0
        return (qcd_priority, ff_priority, pthat_priority)

    def _sort_key(item):
        base_regime, mode, is_ff, qcd_mode, pthat_min, _ = item
        regime_order = {"kaon": 0, "charm": 1, "beauty": 2, "Bc": 3, "ew": 4, "all": 5, "combined": 5}
        mode_order = {None: 0, "direct": 1, "fromTau": 2}
        variant_order = _variant_priority(base_regime, is_ff, qcd_mode, pthat_min)
        return (regime_order.get(base_regime, 99), mode_order.get(mode, 99), base_regime, mode or "", variant_order)

    all_regimes = {"all", "combined"}
    selected_by_mass: Dict[Tuple[float, str], List[Tuple[Path, str]]] = {}
    for key, items in files_by_mass.items():
        all_chan = [it for it in items if it[0] in all_regimes]
        if all_chan:
            chosen = max(all_chan, key=lambda it: _variant_priority(it[0], it[2], it[3], it[4]))
            selected_by_mass[key] = [
                (chosen[5], _label(chosen[0], chosen[1], chosen[2], chosen[3], chosen[4]))
            ]
            continue

        chosen = {}
        for base_regime, mode, is_ff, qcd_mode, pthat_min, path in items:
            k2 = (base_regime, mode)
            if k2 not in chosen:
                chosen[k2] = (base_regime, mode, is_ff, qcd_mode, pthat_min, path)
                continue
            current = chosen[k2]
            new_priority = _variant_priority(base_regime, is_ff, qcd_mode, pthat_min)
            old_priority = _variant_priority(current[0], current[2], current[3], current[4])
            if new_priority > old_priority:
                chosen[k2] = (base_regime, mode, is_ff, qcd_mode, pthat_min, path)

        selected = [v for _, v in sorted(chosen.items(), key=lambda kv: _sort_key(kv[1]))]
        for k2 in {(base_regime, mode) for base_regime, mode, *_ in items}:
            candidates = [it for it in items if (it[0], it[1]) == k2]
            if len(candidates) <= 1:
                continue
            kept = chosen[k2]
            kept_path = str(kept[5])
            dropped = [c for c in candidates if str(c[5]) != kept_path]
            if dropped and not allow_variant_drop:
                raise ValueError(
                    f"Multiple variants for {k2} at m={key[0]:.2f}: "
                    f"keeping {_label(*kept[:5])}, would drop "
                    f"{[_label(*d[:5]) for d in dropped]}. "
                    f"Pass --allow-variant-drop to override."
                )
        selected_by_mass[key] = [
            (path, _label(base_regime, mode, is_ff, qcd_mode, pthat_min))
            for base_regime, mode, is_ff, qcd_mode, pthat_min, path in selected
        ]

    return selected_by_mass


def load_geom_for_sim_files(
    sim_files: List[Tuple[Path, str]],
    mass_str: str,
    flavour: str,
    show_progress: bool | None = False,
    max_signal_events: int | None = None,
    geometry_config: GeometryConfig | None = None,
) -> pd.DataFrame:
    geom_dfs: List[pd.DataFrame] = []
    mesh = None
    geometry_cfg = normalize_geometry_config(geometry_config)
    geom_tag = geometry_tag(geometry_cfg)
    default_geometry = is_default_geometry_config(geometry_cfg)

    for sim_csv, regime in sim_files:
        geom_cache_name = f"{sim_csv.stem}_geom_{geom_tag}.csv"
        geom_csv = limits_run.GEOM_CACHE_DIR / geom_cache_name

        if default_geometry and not geom_csv.exists():
            for legacy_stem in (
                sim_csv.stem + "_geom",
                f"HNL_{mass_str}GeV_{flavour}_combined_geom",
                f"HNL_{mass_str}GeV_{flavour}_geom",
            ):
                legacy_path = limits_run.GEOM_CACHE_DIR / f"{legacy_stem}.csv"
                if legacy_path.exists():
                    geom_csv = legacy_path
                    break

        if geom_csv.exists() and sim_csv.stat().st_mtime > geom_csv.stat().st_mtime:
            geom_csv.unlink()

        if geom_csv.exists():
            geom_df = limits_run._load_geom_cached(sim_csv, geom_csv)
        else:
            if mesh is None:
                mesh = build_drainage_gallery_mesh(geometry_cfg)
            geom_df = preprocess_hnl_csv(sim_csv, mesh, show_progress=show_progress)
            limits_run._save_geom_cache(geom_df, geom_csv)

        geom_df = limits_run._attach_sim_metadata(geom_df, sim_csv)
        if max_signal_events is not None and max_signal_events > 0 and len(geom_df) > max_signal_events:
            geom_df = geom_df.sample(n=max_signal_events, random_state=42)
        geom_dfs.append(geom_df)

    if not geom_dfs:
        return pd.DataFrame()
    return pd.concat(geom_dfs, ignore_index=True)


def available_mass_points(
    flavour: str,
    selected_masses: List[float] | None = None,
    max_mass: float | None = None,
    allow_variant_drop: bool = False,
) -> List[Tuple[float, str, List[Tuple[Path, str]]]]:
    selected = discover_sim_files_for_flavour(
        flavour=flavour,
        selected_masses=selected_masses,
        max_mass=max_mass,
        allow_variant_drop=allow_variant_drop,
    )
    return [
        (mass_val, mass_str, selected[(mass_val, mass_str)])
        for mass_val, mass_str in sorted(selected.keys(), key=lambda x: x[0])
    ]
