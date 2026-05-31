"""Sanity checks on the FONLL table parser using the vendored FONLL tables."""

import json
import os
import subprocess
import sys

import numpy as np
import pytest

from production.fonll.fonll_parser import (
    FONLL_FILES, FONLL_FILE_SETS, parse_fonll_file, get_sigma_total,
)


def test_table_files_exist():
    for label, files in FONLL_FILE_SETS.items():
        for q, path in files.items():
            assert path.exists(), f"Missing vendored {label} FONLL table for {q}: {path}"


def test_parse_shape_is_rectangular_and_covers_target_range():
    for files in FONLL_FILE_SETS.values():
        for q in ("charm", "bottom"):
            pt, y, ds = parse_fonll_file(files[q])
            assert pt.ndim == 1
            assert y.ndim == 1
            assert ds.shape == (100, 100)
            assert pt.min() <= 0.0 and pt.max() >= 50.0
            assert y.min() <= -3.0 and y.max() >= 3.0


def test_nnpdf40_grid_matches_dense_production_contract():
    for q in ("charm", "bottom"):
        pt, y, ds = parse_fonll_file(FONLL_FILE_SETS["nnpdf40_nlo"][q])
        assert pt.ndim == 1
        assert y.ndim == 1
        assert ds.shape == (100, 100)
        assert np.allclose(pt, np.linspace(0.0, 50.0, 100))
        assert np.allclose(y, np.linspace(-3.0, 3.0, 100))
        assert np.isfinite(ds).all()


def test_parser_rejects_rows_outside_pt_outer_y_inner_order(tmp_path):
    path = tmp_path / "shuffled.dat"
    np.savetxt(path, np.array([
        [0.0, -1.0, 1.0],
        [1.0, -1.0, 2.0],
        [0.0, 1.0, 3.0],
        [1.0, 1.0, 4.0],
    ]))
    with pytest.raises(ValueError, match="pT outermost and y innermost"):
        parse_fonll_file(path)


def test_parser_accepts_consistent_rounding(tmp_path):
    """A rectangular table with ~1e-12 float noise must pass the new allclose check.

    Also confirms the previous np.array_equal check would have rejected it,
    so the regression value of relaxing to np.allclose is explicit.
    """
    path = tmp_path / "noisy.dat"
    pt_nodes = np.array([0.0, 1.0, 2.0])
    y_nodes = np.array([-1.0, 0.0, 1.0])
    rows = []
    for i, pt in enumerate(pt_nodes):
        for j, y in enumerate(y_nodes):
            # Add subtle 1e-12 noise to pT; y stays clean.
            noisy_pt = pt + 1e-12 * ((i + j) % 2)
            rows.append([noisy_pt, y, 1.0 + 0.1 * i + 0.01 * j])
    data = np.array(rows)
    np.savetxt(path, data)

    # Sanity: np.array_equal would NOT have accepted this — the noise is real.
    pt_all = data[:, 0]
    pt_unique = np.unique(pt_all)
    expected_pt = np.repeat(pt_unique, len(y_nodes))
    assert np.array_equal(pt_all, expected_pt) is False

    # The new allclose-based parser must accept the noisy table.
    pt_arr, y_arr, ds = parse_fonll_file(path)
    assert pt_arr.shape == (3,)
    assert y_arr.shape == (3,)
    assert ds.shape == (3, 3)


def test_parser_rejects_descending_y(tmp_path):
    """A file with descending y rows must raise a clear ascending-order error."""
    path = tmp_path / "descending_y.dat"
    pt_nodes = np.array([0.0, 1.0])
    y_nodes_desc = np.array([1.0, 0.0, -1.0])  # descending
    rows = []
    for pt in pt_nodes:
        for y in y_nodes_desc:
            rows.append([pt, y, 1.0])
    np.savetxt(path, np.array(rows))
    with pytest.raises(ValueError, match="y column must be in ascending order"):
        parse_fonll_file(path)


def _run_dispatch(env_value=None):
    env = os.environ.copy()
    if env_value is None:
        env.pop("HNL_FONLL_SET", None)
    else:
        env["HNL_FONLL_SET"] = env_value
    code = """
import json
from production.fonll.fonll_parser import FONLL_DEFAULT_SET, FONLL_FILES, get_sigma_total
from production.fonll.meson_sampler import FONLL_FILES as SAMPLER_FILES
print(json.dumps({
    "set": FONLL_DEFAULT_SET,
    "bottom": FONLL_FILES["bottom"].name,
    "sampler_bottom": SAMPLER_FILES["bottom"].name,
    "sigma_bottom": get_sigma_total("bottom"),
}))
"""
    return subprocess.run(
        [sys.executable, "-c", code], env=env, check=False,
        capture_output=True, text=True,
    )


def test_backend_dispatch_default_and_legacy():
    default = _run_dispatch()
    assert default.returncode == 0, default.stderr
    default_data = json.loads(default.stdout)
    assert default_data["set"] == "nnpdf40_nlo"
    assert "nnpdf40_nlo_as_01180" in default_data["bottom"]
    assert default_data["sampler_bottom"] == default_data["bottom"]
    assert default_data["sigma_bottom"] > 0

    legacy = _run_dispatch("cteq66_legacy")
    assert legacy.returncode == 0, legacy.stderr
    legacy_data = json.loads(legacy.stdout)
    assert legacy_data["set"] == "cteq66_legacy"
    assert "cteq66" in legacy_data["bottom"]
    assert legacy_data["sampler_bottom"] == legacy_data["bottom"]
    assert legacy_data["sigma_bottom"] > 0


def test_invalid_backend_fails_fast():
    result = _run_dispatch("typo")
    assert result.returncode != 0
    assert "unknown HNL_FONLL_SET='typo'" in result.stderr
    assert "cteq66_legacy, nnpdf40_nlo" in result.stderr


def test_dsigma_nonnegative_majority():
    """At least 95% of bins should be non-negative (small numerical negatives allowed near edges)."""
    for q in ("charm", "bottom"):
        _, _, ds = parse_fonll_file(FONLL_FILES[q])
        frac_nn = np.mean(ds >= 0)
        assert frac_nn > 0.95, f"{q}: only {frac_nn:.1%} non-negative bins"


def test_total_sigma_order_of_magnitude():
    """Total integrated cross section should match expected LHC 14 TeV ballpark."""
    # FONLL at 14 TeV: sigma_cc ~ 1e10 pb, sigma_bb ~ 1e8 pb (very approximate)
    sigma_c = get_sigma_total("charm")
    sigma_b = get_sigma_total("bottom")
    assert 1e9 < sigma_c < 1e11, f"sigma_charm = {sigma_c:.3e} pb out of expected range"
    assert 1e7 < sigma_b < 1e9, f"sigma_bottom = {sigma_b:.3e} pb out of expected range"
    # charm > bottom at LHC
    assert sigma_c > sigma_b
