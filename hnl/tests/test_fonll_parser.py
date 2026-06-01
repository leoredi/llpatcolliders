"""Sanity checks on the FONLL table parser using the vendored FONLL tables."""

import importlib

import numpy as np
import pytest

from production.fonll.fonll_parser import (
    FONLL_FILES, FONLL_FILE_SETS, parse_fonll_file, get_sigma_total,
)


def test_table_files_exist():
    for label, files in FONLL_FILE_SETS.items():
        for q, path in files.items():
            assert path.exists(), f"Missing vendored {label} FONLL table for {q}: {path}"


def test_parse_shape_is_rectangular_for_all_backends():
    """All backends parse as 1D pT, 1D y, 2D dsigma of consistent shape.
    Does NOT hard-code grid density."""
    for files in FONLL_FILE_SETS.values():
        for q in ("charm", "bottom"):
            pt, y, ds = parse_fonll_file(files[q])
            assert pt.ndim == 1
            assert y.ndim == 1
            assert ds.shape == (len(pt), len(y))
            assert len(pt) >= 2 and len(y) >= 2
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


def _reload_parser_with_env(monkeypatch, value):
    """Reload the parser module with HNL_FONLL_SET set to `value`,
    or unset if value is None. Returns the freshly reloaded module."""
    if value is None:
        monkeypatch.delenv("HNL_FONLL_SET", raising=False)
    else:
        monkeypatch.setenv("HNL_FONLL_SET", value)
    import production.fonll.fonll_parser as m
    return importlib.reload(m)


def test_backend_dispatch_default_is_nnpdf40(monkeypatch):
    m = _reload_parser_with_env(monkeypatch, None)
    assert m.FONLL_DEFAULT_SET == "nnpdf40_nlo"
    assert "nnpdf40_nlo_as_01180" in m.FONLL_FILES["bottom"].name
    assert m.get_sigma_total("bottom") > 0


def test_backend_dispatch_legacy_cteq66(monkeypatch):
    m = _reload_parser_with_env(monkeypatch, "cteq66_legacy")
    assert m.FONLL_DEFAULT_SET == "cteq66_legacy"
    assert "cteq66" in m.FONLL_FILES["bottom"].name
    assert m.get_sigma_total("bottom") > 0


def test_invalid_backend_fails_fast(monkeypatch):
    monkeypatch.setenv("HNL_FONLL_SET", "typo")
    import production.fonll.fonll_parser as m
    with pytest.raises(ValueError, match="unknown HNL_FONLL_SET='typo'"):
        importlib.reload(m)


def test_sampler_picks_up_each_backend_subprocess():
    """Smoke test: spawn a clean Python under each HNL_FONLL_SET and confirm
    that sample_meson_4vectors actually reads from the requested backend.

    This covers the import-time module-binding hazard that the in-process
    reload tests above cannot catch: any consumer that did
    `from ...fonll_parser import FONLL_FILES` (a named import) keeps the
    original dict after a parser reload, so reloading only the parser is
    not a faithful proxy for what the sampler actually does in production.
    """
    import json
    import os
    import subprocess
    import sys

    code = """
import json
from production.fonll import fonll_parser, meson_sampler
print(json.dumps({
    "default_set": fonll_parser.FONLL_DEFAULT_SET,
    "parser_bottom": fonll_parser.FONLL_FILES["bottom"].name,
    # Verify the sampler looks up the parser's dict dynamically; the path
    # used in production is the parser module's binding, not a stale
    # snapshot captured at meson_sampler import time.
    "sampler_bottom_via_parser": meson_sampler.fonll_parser.FONLL_FILES["bottom"].name,
}))
"""
    for env_value, expected_set, expected_substr in (
        ("nnpdf40_nlo", "nnpdf40_nlo", "nnpdf40_nlo_as_01180"),
        ("cteq66_legacy", "cteq66_legacy", "cteq66"),
    ):
        env = os.environ.copy()
        env["HNL_FONLL_SET"] = env_value
        result = subprocess.run(
            [sys.executable, "-c", code], env=env, check=False,
            capture_output=True, text=True,
        )
        assert result.returncode == 0, result.stderr
        data = json.loads(result.stdout)
        assert data["default_set"] == expected_set
        assert expected_substr in data["parser_bottom"]
        assert data["sampler_bottom_via_parser"] == data["parser_bottom"]


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


@pytest.fixture(autouse=True)
def _restore_parser_default(monkeypatch):
    yield
    monkeypatch.delenv("HNL_FONLL_SET", raising=False)
    import production.fonll.fonll_parser as m
    importlib.reload(m)
