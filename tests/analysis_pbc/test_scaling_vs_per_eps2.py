#!/usr/bin/env python3

from __future__ import annotations

import math
import random
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
import sys
REPO_ROOT = HERE.parents[2]
ANALYSIS_ROOT = REPO_ROOT / "analysis_pbc"
if str(ANALYSIS_ROOT) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_ROOT))

from limits.expected_signal import expected_signal_events, couplings_from_eps2
from hnl_models.hnl_model_hnlcalc import HNLModel


def _make_geom_df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "parent_id": [511, 521, 511],
            "weight": [1.0, 1.0, 2.0],
            "beta_gamma": [2.0, 5.0, 3.0],
            "hits_tube": [True, True, False],
            "entry_distance": [10.0, 5.0, 0.0],
            "path_length": [2.0, 1.0, 0.0],
            "eta": [0.1, -0.2, 0.3],
            "phi": [0.0, 1.0, 2.0],
        }
    )


def test_scaling_matches_per_eps2():
    geom_df = _make_geom_df()
    separation_pass = np.ones(len(geom_df), dtype=bool)

    mass = 2.6
    benchmark = "010"
    eps2 = 1e-8
    eps2_ref = 1e-6
    lumi_fb = 3000.0
    separation_m = 1e-3

    random.seed(123)
    n_direct = expected_signal_events(
        geom_df=geom_df,
        mass_GeV=mass,
        eps2=eps2,
        benchmark=benchmark,
        lumi_fb=lumi_fb,
        separation_m=separation_m,
        separation_pass=separation_pass,
    )

    random.seed(123)
    Ue2, Umu2, Utau2 = couplings_from_eps2(eps2_ref, benchmark)
    model = HNLModel(mass_GeV=mass, Ue2=Ue2, Umu2=Umu2, Utau2=Utau2)
    ctau0_ref = model.ctau0_m
    br_ref = model.production_brs()

    n_scaled = expected_signal_events(
        geom_df=geom_df,
        mass_GeV=mass,
        eps2=eps2,
        benchmark=benchmark,
        lumi_fb=lumi_fb,
        separation_m=separation_m,
        separation_pass=separation_pass,
        ctau0_m=ctau0_ref * (eps2_ref / eps2),
        br_per_parent=br_ref,
        br_scale=eps2 / eps2_ref,
    )

    if n_direct == 0.0:
        assert n_scaled == 0.0
        return

    rel_err = abs(n_scaled - n_direct) / abs(n_direct)
    assert rel_err < 1e-10, f"Scaling mismatch: n_direct={n_direct}, n_scaled={n_scaled}, rel_err={rel_err}"


if __name__ == "__main__":
    test_scaling_matches_per_eps2()
    print("OK: scaling matches per-eps2 for single test point.")
