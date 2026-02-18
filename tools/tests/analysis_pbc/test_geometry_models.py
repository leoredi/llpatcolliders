#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path
import sys

import pytest

HERE = Path(__file__).resolve()
REPO_ROOT = HERE.parents[3]
ANALYSIS_ROOT = REPO_ROOT / "analysis_pbc"
if str(ANALYSIS_ROOT) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_ROOT))

from geometry.per_parent_efficiency import (  # noqa: E402
    GeometryConfig,
    build_drainage_gallery_mesh,
    geometry_tag,
    normalize_geometry_config,
)


def test_geometry_tag_is_deterministic():
    cfg = normalize_geometry_config(GeometryConfig(model="tube", tube_radius_m=1.54))
    tag_a = geometry_tag(cfg)
    tag_b = geometry_tag(cfg)
    assert tag_a == tag_b
    assert len(tag_a.split("_")[-1]) == 8


def test_geometry_tag_changes_with_config():
    cfg_tube = normalize_geometry_config(GeometryConfig(model="tube", tube_radius_m=1.54))
    cfg_profile = normalize_geometry_config(
        GeometryConfig(model="profile", detector_thickness_m=0.24, profile_inset_floor=False)
    )
    assert geometry_tag(cfg_tube) != geometry_tag(cfg_profile)


@pytest.mark.parametrize(
    "cfg",
    [
        GeometryConfig(model="tube", tube_radius_m=1.54),
        GeometryConfig(model="profile", detector_thickness_m=0.24, profile_inset_floor=False),
    ],
)
def test_geometry_mesh_is_watertight_and_positive_volume(cfg: GeometryConfig):
    mesh = build_drainage_gallery_mesh(cfg)
    assert mesh.volume > 0.0
    assert mesh.is_watertight
