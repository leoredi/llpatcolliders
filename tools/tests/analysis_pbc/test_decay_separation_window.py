#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

HERE = Path(__file__).resolve()
REPO_ROOT = HERE.parents[3]
ANALYSIS_ROOT = REPO_ROOT / "analysis_pbc"
if str(ANALYSIS_ROOT) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_ROOT))

from decay.decay_detector import (  # noqa: E402
    DecayCache,
    DecaySelection,
    compute_decay_acceptance,
    compute_separation_pass_static,
    pairwise_separation_pass,
)


def test_all_pairs_min_baseline_behavior():
    pts_pass = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.010],
            [0.0, 0.010, 0.0],
        ],
        dtype=float,
    )
    pts_fail = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0005],
            [0.0, 0.010, 0.0],
        ],
        dtype=float,
    )

    assert pairwise_separation_pass(pts_pass, min_separation_m=1e-3)
    assert not pairwise_separation_pass(pts_fail, min_separation_m=1e-3)


def test_any_pair_window_behavior():
    pts = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.001],
            [0.0, 0.0, 0.050],
        ],
        dtype=float,
    )
    assert pairwise_separation_pass(
        pts,
        min_separation_m=1e-3,
        max_separation_m=2e-2,
        separation_policy="any-pair-window",
    )
    assert not pairwise_separation_pass(
        pts,
        min_separation_m=6e-2,
        max_separation_m=8e-2,
        separation_policy="any-pair-window",
    )


def test_max_separation_edges():
    pts = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.010],
            [0.0, 0.030, 0.0],
        ],
        dtype=float,
    )
    assert not pairwise_separation_pass(
        pts,
        min_separation_m=1e-3,
        max_separation_m=2e-2,
        separation_policy="all-pairs-min",
    )
    assert pairwise_separation_pass(
        pts,
        min_separation_m=1e-3,
        max_separation_m=5e-2,
        separation_policy="all-pairs-min",
    )


def test_invalid_policy_rejected():
    pts = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.01]], dtype=float)
    with pytest.raises(ValueError, match="Unsupported separation policy"):
        pairwise_separation_pass(pts, min_separation_m=1e-3, separation_policy="bad-policy")


def _one_hit_geom_df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "hits_tube": [True],
            "entry_distance": [10.0],
            "path_length": [2.0],
            "beta_gamma": [4.0],
            "eta": [0.2],
            "phi": [0.1],
        }
    )


def _one_hit_decay_cache() -> DecayCache:
    return DecayCache(
        charged_directions=[[np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0])]],
        decay_u=np.array([0.5], dtype=float),
        hit_indices=np.array([0], dtype=int),
    )


def test_dynamic_static_structural_consistency(monkeypatch: pytest.MonkeyPatch):
    def _fake_intersections(mesh, origin, directions):
        return [np.array([0.0, 0.0, 0.0]), np.array([0.0, 0.0, 0.015])]

    import decay.decay_detector as dd  # noqa: WPS433

    monkeypatch.setattr(dd, "_batch_first_intersections", _fake_intersections)

    geom_df = _one_hit_geom_df()
    cache = _one_hit_decay_cache()
    selection = DecaySelection(
        separation_m=1e-3,
        max_separation_m=2e-2,
        separation_policy="all-pairs-min",
        seed=123,
        p_min_GeV=0.6,
    )

    dyn = compute_decay_acceptance(
        geom_df=geom_df,
        mass_GeV=4.0,
        flavour="electron",
        ctau0_m=1.0,
        mesh=object(),
        selection=selection,
        decay_cache=cache,
    )
    sta = compute_separation_pass_static(
        geom_df=geom_df,
        decay_cache=cache,
        mesh=object(),
        separation_m=1e-3,
        max_separation_m=2e-2,
        separation_policy="all-pairs-min",
    )

    assert bool(dyn[0]) is True
    assert bool(sta[0]) is True
