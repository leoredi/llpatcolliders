#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd
import pytest

HERE = Path(__file__).resolve()
REPO_ROOT = HERE.parents[3]
ANALYSIS_ROOT = REPO_ROOT / "analysis_pbc"
if str(ANALYSIS_ROOT) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_ROOT))

from decay import brvis_kappa
from geometry.per_parent_efficiency import GeometryConfig, geometry_tag, normalize_geometry_config


@pytest.fixture(autouse=True)
def _clear_cache():
    brvis_kappa._load_kappa_table_cached.cache_clear()
    yield
    brvis_kappa._load_kappa_table_cached.cache_clear()


def _write_table(path: Path, rows: list[dict]) -> None:
    df = pd.DataFrame(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def test_lookup_exact_match(tmp_path: Path):
    path = tmp_path / "kappa.csv"
    _write_table(
        path,
        [
            {
                "flavour": "electron",
                "mass_GeV": 2.0,
                "kappa": 0.91,
                "p_min_GeV": 0.6,
                "separation_mm": 1.0,
                "source_policy": "hybrid",
                "status": "ok",
            }
        ],
    )

    value = brvis_kappa.lookup_kappa("electron", 2.0, 0.6, 1.0, path)
    assert value == pytest.approx(0.91)


def test_lookup_interpolates_between_bracketing_points(tmp_path: Path):
    path = tmp_path / "kappa.csv"
    _write_table(
        path,
        [
            {
                "flavour": "electron",
                "mass_GeV": 2.0,
                "kappa": 0.8,
                "p_min_GeV": 0.6,
                "separation_mm": 1.0,
                "source_policy": "hybrid",
                "status": "ok",
            },
            {
                "flavour": "electron",
                "mass_GeV": 4.0,
                "kappa": 1.2,
                "p_min_GeV": 0.6,
                "separation_mm": 1.0,
                "source_policy": "hybrid",
                "status": "ok",
            },
        ],
    )

    value = brvis_kappa.lookup_kappa("electron", 3.0, 0.6, 1.0, path)
    assert value == pytest.approx(1.0)


def test_lookup_rejects_cut_metadata_mismatch(tmp_path: Path):
    path = tmp_path / "kappa.csv"
    _write_table(
        path,
        [
            {
                "flavour": "muon",
                "mass_GeV": 3.0,
                "kappa": 1.0,
                "p_min_GeV": 0.6,
                "separation_mm": 1.0,
                "source_policy": "hybrid",
                "status": "ok",
            }
        ],
    )

    with pytest.raises(brvis_kappa.KappaTableError, match="p_min mismatch"):
        brvis_kappa.lookup_kappa("muon", 3.0, 0.5, 1.0, path)


def test_lookup_rejects_mass_outside_calibrated_range(tmp_path: Path):
    path = tmp_path / "kappa.csv"
    _write_table(
        path,
        [
            {
                "flavour": "tau",
                "mass_GeV": 4.0,
                "kappa": 0.9,
                "p_min_GeV": 0.6,
                "separation_mm": 1.0,
                "source_policy": "hybrid",
                "status": "ok",
            },
            {
                "flavour": "tau",
                "mass_GeV": 6.0,
                "kappa": 1.1,
                "p_min_GeV": 0.6,
                "separation_mm": 1.0,
                "source_policy": "hybrid",
                "status": "ok",
            },
        ],
    )

    with pytest.raises(brvis_kappa.KappaTableError, match="out of calibrated range"):
        brvis_kappa.lookup_kappa("tau", 7.0, 0.6, 1.0, path)


def test_lookup_ignores_non_ok_rows_and_fails_when_none_ok(tmp_path: Path):
    path = tmp_path / "kappa.csv"
    _write_table(
        path,
        [
            {
                "flavour": "electron",
                "mass_GeV": 4.0,
                "kappa": 1.0,
                "p_min_GeV": 0.6,
                "separation_mm": 1.0,
                "source_policy": "hybrid",
                "status": "fail",
            }
        ],
    )

    with pytest.raises(brvis_kappa.KappaTableError, match="no rows with status=ok"):
        brvis_kappa.lookup_kappa("electron", 4.0, 0.6, 1.0, path)


def test_lookup_rejects_geometry_tag_mismatch(tmp_path: Path):
    path = tmp_path / "kappa.csv"
    cfg_default = normalize_geometry_config()
    _write_table(
        path,
        [
            {
                "flavour": "electron",
                "mass_GeV": 3.0,
                "kappa": 1.0,
                "p_min_GeV": 0.6,
                "separation_mm": 1.0,
                "source_policy": "hybrid",
                "status": "ok",
                "geometry_tag": geometry_tag(cfg_default),
                "geometry_model": cfg_default.model,
                "tube_radius_m": cfg_default.tube_radius_m,
                "detector_thickness_m": cfg_default.detector_thickness_m,
                "profile_inset_floor": str(cfg_default.profile_inset_floor).lower(),
            }
        ],
    )

    cfg_profile = normalize_geometry_config(
        GeometryConfig(model="profile", detector_thickness_m=0.24, profile_inset_floor=False)
    )
    with pytest.raises(brvis_kappa.KappaTableError, match="geometry_tag mismatch"):
        brvis_kappa.lookup_kappa("electron", 3.0, 0.6, 1.0, path, geometry_config=cfg_profile)


def test_lookup_legacy_table_allows_only_default_geometry(tmp_path: Path):
    path = tmp_path / "kappa.csv"
    _write_table(
        path,
        [
            {
                "flavour": "muon",
                "mass_GeV": 4.0,
                "kappa": 0.95,
                "p_min_GeV": 0.6,
                "separation_mm": 1.0,
                "source_policy": "hybrid",
                "status": "ok",
            }
        ],
    )

    # Legacy table: default geometry remains valid.
    val = brvis_kappa.lookup_kappa("muon", 4.0, 0.6, 1.0, path)
    assert val == pytest.approx(0.95)

    non_default_tube = normalize_geometry_config(GeometryConfig(model="tube", tube_radius_m=1.60))
    with pytest.raises(brvis_kappa.KappaTableError, match="legacy default-only"):
        brvis_kappa.lookup_kappa(
            "muon",
            4.0,
            0.6,
            1.0,
            path,
            geometry_config=non_default_tube,
        )
