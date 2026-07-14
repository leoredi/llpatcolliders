import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from scalar.production import MASS_GRID


PUBLISHED = Path(__file__).resolve().parents[1] / "data" / "published"


def _sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_published_curve_and_manifest_are_self_consistent():
    curve_path = PUBLISHED / "bc4_island.csv"
    manifest = json.loads((PUBLISHED / "MANIFEST.json").read_text())
    frame = pd.read_csv(curve_path)

    assert _sha256(curve_path) == manifest["csv_sha256"]
    assert not frame["mass_GeV"].duplicated().any()
    assert frame["mass_GeV"].is_monotonic_increasing
    assert np.array_equal(
        frame["mass_GeV"].to_numpy(float), np.asarray(MASS_GRID, dtype=float)
    )

    headline = manifest["headline_reach"]
    sensitive = frame[frame["has_sensitivity"].astype(bool)]
    best = sensitive.loc[sensitive["u2_min"].idxmin()]
    assert len(frame) == headline["n_grid_points"]
    assert len(sensitive) == headline["n_sensitive_points"]
    assert best["mass_GeV"] == pytest.approx(
        headline["best_sin2theta_min_at_mass_GeV"]
    )
    assert best["u2_min"] == pytest.approx(
        headline["best_sin2theta_min"], rel=1e-12
    )
    assert sensitive["mass_GeV"].min() == pytest.approx(
        headline["mass_lo_grid_GeV"]
    )
    assert sensitive["mass_GeV"].max() == pytest.approx(
        headline["mass_hi_grid_GeV"]
    )

    last = sensitive.index[-1]
    following = last + 1
    fraction = (
        np.log(3.0) - np.log(frame.loc[last, "peak_N"])
    ) / (
        np.log(frame.loc[following, "peak_N"])
        - np.log(frame.loc[last, "peak_N"])
    )
    closure = frame.loc[last, "mass_GeV"] + fraction * (
        frame.loc[following, "mass_GeV"] - frame.loc[last, "mass_GeV"]
    )
    assert closure == pytest.approx(headline["mass_hi_closure_GeV"], rel=1e-12)

    upper_open = frame[frame["u2_max_open"].astype(bool)]
    assert [upper_open["mass_GeV"].min(), upper_open["mass_GeV"].max()] == (
        manifest["topology"]["upper_open_grid_run_GeV"]
    )
