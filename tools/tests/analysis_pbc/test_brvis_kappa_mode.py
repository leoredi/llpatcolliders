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

from config.production_xsecs import get_parent_sigma_pb
from limits import expected_signal as es


def _geom_df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "parent_id": [511, 511],
            "weight": [1.0, 1.0],
            "beta_gamma": [2.0, 2.0],
            "hits_tube": [True, False],
            "entry_distance": [10.0, 0.0],
            "path_length": [2.0, 0.0],
            "eta": [0.1, 0.2],
            "phi": [0.0, 0.1],
        }
    )


def test_brvis_kappa_mode_skips_decay_library_calls(monkeypatch: pytest.MonkeyPatch):
    def _boom(*args, **kwargs):
        raise AssertionError("compute_decay_acceptance should not be called in brvis_kappa mode")

    monkeypatch.setattr(es, "compute_decay_acceptance", _boom)

    n_sig = es.expected_signal_events(
        geom_df=_geom_df(),
        mass_GeV=4.0,
        eps2=1e-6,
        benchmark="100",
        lumi_fb=1.0,
        separation_m=1e-3,
        decay_mode="brvis_kappa",
        br_vis=0.8,
        kappa_eff=0.5,
        ctau0_m=1.0,
        br_per_parent={511: 0.1},
        br_scale=1.0,
    )

    assert np.isfinite(n_sig)
    assert n_sig > 0.0


def test_brvis_kappa_scaling_matches_manual_expectation():
    geom_df = _geom_df()

    br_vis = 0.8
    kappa = 0.5
    scale = br_vis * kappa
    ctau0_m = 1.0

    n_sig = es.expected_signal_events(
        geom_df=geom_df,
        mass_GeV=4.0,
        eps2=1e-6,
        benchmark="100",
        lumi_fb=1.0,
        separation_m=1e-3,
        decay_mode="brvis_kappa",
        br_vis=br_vis,
        kappa_eff=kappa,
        ctau0_m=ctau0_m,
        br_per_parent={511: 0.1},
        br_scale=1.0,
    )

    sigma_pb = get_parent_sigma_pb(511)
    lam = 2.0 * ctau0_m
    p_decay_hit = np.exp(-10.0 / lam) * (1.0 - np.exp(-2.0 / lam))
    eff_parent = (p_decay_hit * scale + 0.0) / 2.0
    expected = (sigma_pb * 1e3) * 0.1 * eff_parent

    assert n_sig == pytest.approx(expected, rel=1e-12, abs=0.0)


def test_brvis_kappa_requires_br_vis_and_kappa():
    geom_df = _geom_df()

    with pytest.raises(ValueError, match="requires finite br_vis"):
        es.expected_signal_events(
            geom_df=geom_df,
            mass_GeV=4.0,
            eps2=1e-6,
            benchmark="100",
            lumi_fb=1.0,
            separation_m=1e-3,
            decay_mode="brvis_kappa",
            br_vis=None,
            kappa_eff=1.0,
            ctau0_m=1.0,
            br_per_parent={511: 0.1},
            br_scale=1.0,
        )

    with pytest.raises(ValueError, match="requires finite kappa_eff"):
        es.expected_signal_events(
            geom_df=geom_df,
            mass_GeV=4.0,
            eps2=1e-6,
            benchmark="100",
            lumi_fb=1.0,
            separation_m=1e-3,
            decay_mode="brvis_kappa",
            br_vis=1.0,
            kappa_eff=None,
            ctau0_m=1.0,
            br_per_parent={511: 0.1},
            br_scale=1.0,
        )


def test_brvis_kappa_ignores_decay_cache_and_separation_pass_inputs():
    geom_df = _geom_df()

    kwargs = dict(
        geom_df=geom_df,
        mass_GeV=4.0,
        eps2=1e-6,
        benchmark="100",
        lumi_fb=1.0,
        separation_m=1e-3,
        decay_mode="brvis_kappa",
        br_vis=0.8,
        kappa_eff=0.5,
        ctau0_m=1.0,
        br_per_parent={511: 0.1},
        br_scale=1.0,
    )

    n_base = es.expected_signal_events(**kwargs)
    n_with_library_artifacts = es.expected_signal_events(
        **kwargs,
        decay_cache=object(),
        separation_pass=np.array([False, False], dtype=bool),
    )

    assert n_with_library_artifacts == pytest.approx(n_base, rel=0.0, abs=0.0)
