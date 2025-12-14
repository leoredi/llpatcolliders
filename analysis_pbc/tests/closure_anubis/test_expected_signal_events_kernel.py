"""
Algorithmic closure tests for limits.expected_signal_events.

These tests do NOT use the real HNL model or cross-sections.
They patch HNLModel and get_parent_sigma_pb inside
limits.expected_signal so that expected_signal_events
should reproduce a known analytic probability.

Run with:
  cd analysis_pbc
  python tests/closure_anubis/test_expected_signal_events_kernel.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd

# ----------------------------------------------------------------------
# 0. Ensure analysis_pbc is on sys.path
# ----------------------------------------------------------------------

THIS_FILE = Path(__file__).resolve()
# .../analysis_pbc/tests/closure_anubis/test_expected_signal_events_kernel.py
ANALYSIS_DIR = THIS_FILE.parents[2]  # .../analysis_pbc
if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))

from limits.expected_signal import expected_signal_events


# ----------------------------------------------------------------------
# 1. Single-HNL shell test
# ----------------------------------------------------------------------

def test_single_hnl_segment() -> None:
    """
    Check that for a single HNL in a segment [d, d+L] with known λ,
    expected_signal_events reproduces:

        P_analytic = exp(-d / λ) * (1 - exp(-L / λ))

    when we normalise such that L_int * σ * BR = 1.
    """
    # Choose convenient parameters
    beta_gamma = 10.0           # dimensionless
    lam = 20.0                  # mean decay length in metres
    ctau0_m = lam / beta_gamma  # so that λ = βγ * cτ₀

    d = 5.0                     # entry distance [m]
    L = 10.0                    # path length [m] inside detection region

    # Analytic probability
    P_analytic = math.exp(-d / lam) * (1.0 - math.exp(-L / lam))

    # Build geometry DataFrame: single HNL from parent 521
    geom_df = pd.DataFrame(
        {
            "parent_id": [521.0],
            "weight": [1.0],
            "beta_gamma": [beta_gamma],
            "hits_tube": [True],
            "entry_distance": [d],
            "path_length": [L],
        }
    )

    # Choose luminosity so that L_int * σ * BR = 1
    lumi_fb = 1.0 / 1000.0      # fb^-1
    mass_GeV = 1.0
    eps2 = 1.0                  # arbitrary (only used by the real HNLModel)
    benchmark = "010"

    with patch("limits.expected_signal.HNLModel") as MockModel, \
         patch("limits.expected_signal.get_parent_sigma_pb") as mock_sigma:

        # Configure the fake HNL model
        mock_instance = MockModel.return_value
        mock_instance.ctau0_m = ctau0_m
        mock_instance.production_brs.return_value = {521: 1.0}

        # Fake production cross-section
        mock_sigma.return_value = 1.0  # pb

        N_sig = expected_signal_events(
            geom_df=geom_df,
            mass_GeV=mass_GeV,
            eps2=eps2,
            benchmark=benchmark,
            lumi_fb=lumi_fb,
        )

    if not np.isfinite(N_sig):
        raise RuntimeError(f"N_sig is not finite: {N_sig}")

    rel_diff = abs(N_sig - P_analytic) / P_analytic
    print(f"[single] N_sig = {N_sig:.6e}, P_analytic = {P_analytic:.6e}")
    print(f"[single] Relative difference = {rel_diff:.3e}")

    if rel_diff > 1e-3:
        raise AssertionError(
            f"Algorithmic closure failed: rel diff = {rel_diff:.3e} > 1e-3"
        )


# ----------------------------------------------------------------------
# 2. Two-HNL weighted test
# ----------------------------------------------------------------------

def test_two_hnls_weighted_average() -> None:
    """
    Check that with two HNLs and weights w1, w2, the result matches
    the weighted average of the two analytic probabilities.
    """
    beta_gamma = 5.0
    lam = 10.0
    ctau0_m = lam / beta_gamma

    d1, L1, w1 = 3.0, 5.0, 2.0
    d2, L2, w2 = 15.0, 8.0, 1.0

    P1 = math.exp(-d1 / lam) * (1.0 - math.exp(-L1 / lam))
    P2 = math.exp(-d2 / lam) * (1.0 - math.exp(-L2 / lam))

    P_avg = (w1 * P1 + w2 * P2) / (w1 + w2)

    geom_df = pd.DataFrame(
        {
            "parent_id": [521.0, 521.0],
            "weight": [w1, w2],
            "beta_gamma": [beta_gamma, beta_gamma],
            "hits_tube": [True, True],
            "entry_distance": [d1, d2],
            "path_length": [L1, L2],
        }
    )

    lumi_fb = 1.0 / 1000.0
    mass_GeV = 1.0
    eps2 = 1.0
    benchmark = "010"

    with patch("limits.expected_signal.HNLModel") as MockModel, \
         patch("limits.expected_signal.get_parent_sigma_pb") as mock_sigma:

        mock_instance = MockModel.return_value
        mock_instance.ctau0_m = ctau0_m
        mock_instance.production_brs.return_value = {521: 1.0}
        mock_sigma.return_value = 1.0

        N_sig = expected_signal_events(
            geom_df=geom_df,
            mass_GeV=mass_GeV,
            eps2=eps2,
            benchmark=benchmark,
            lumi_fb=lumi_fb,
        )

    rel_diff = abs(N_sig - P_avg) / P_avg
    print(f"[weighted] N_sig = {N_sig:.6e}, P_avg = {P_avg:.6e}")
    print(f"[weighted] Relative difference = {rel_diff:.3e}")

    if rel_diff > 1e-3:
        raise AssertionError(
            f"Weighted closure failed: rel diff = {rel_diff:.3e} > 1e-3"
        )


# ----------------------------------------------------------------------
# 3. Entry point for manual execution
# ----------------------------------------------------------------------

def main() -> None:
    print("=" * 70)
    print("ALGORITHMIC CLOSURE TESTS FOR expected_signal_events")
    print("=" * 70)
    test_single_hnl_segment()
    test_two_hnls_weighted_average()
    print("✓ All algorithmic closure tests passed.\n")


if __name__ == "__main__":
    main()
