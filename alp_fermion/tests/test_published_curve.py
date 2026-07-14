import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from alp_fermion.mass_grid import ALP_MASS_GRID


PUBLISHED = Path(__file__).resolve().parents[1] / "data" / "published"


def _sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _sensitive_runs(frame):
    runs = []
    start = None
    masses = frame["mass_GeV"].to_numpy(float)
    flags = frame["has_sensitivity"].astype(bool).to_numpy()
    for index, sensitive in enumerate(flags):
        if sensitive and start is None:
            start = index
        if start is not None and (not sensitive or index == len(flags) - 1):
            end = index if sensitive and index == len(flags) - 1 else index - 1
            runs.append([masses[start], masses[end]])
            start = None
    return runs


def test_published_curve_and_manifest_are_self_consistent():
    curve_path = PUBLISHED / "bc10_sensitivity.csv"
    manifest = json.loads((PUBLISHED / "MANIFEST.json").read_text())
    frame = pd.read_csv(curve_path)

    assert _sha256(curve_path) == manifest["csv_sha256"]
    assert not frame["mass_GeV"].duplicated().any()
    assert frame["mass_GeV"].is_monotonic_increasing
    assert np.array_equal(
        frame["mass_GeV"].to_numpy(float), np.asarray(ALP_MASS_GRID, dtype=float)
    )

    headline = manifest["headline_reach"]
    sensitive = frame[frame["has_sensitivity"].astype(bool)]
    best = sensitive.loc[sensitive["invf_min"].idxmin()]
    assert len(frame) == headline["n_grid_points"]
    assert len(sensitive) == headline["n_sensitive_points"]
    assert best["mass_GeV"] == pytest.approx(
        headline["best_invf_min_at_mass_GeV"]
    )
    assert best["invf_min"] == pytest.approx(
        headline["best_invf_min_GeV_inv"], rel=1e-12
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
    assert _sensitive_runs(frame) == manifest["topology"]["sensitive_grid_runs_GeV"]
