"""Sanity checks on the FONLL table parser using the vendored CTEQ6.6 tables."""

import numpy as np

from production.fonll.fonll_parser import (
    FONLL_FILES, parse_fonll_file, get_sigma_total,
)


def test_table_files_exist():
    for q, path in FONLL_FILES.items():
        assert path.exists(), f"Missing vendored FONLL table for {q}: {path}"


def test_parse_shape_is_100_by_100():
    for q in ("charm", "bottom"):
        pt, y, ds = parse_fonll_file(FONLL_FILES[q])
        assert pt.shape == (100,)
        assert y.shape == (100,)
        assert ds.shape == (100, 100)
        # pt grid covers (0, 50] GeV, y grid covers [-3, 3]
        assert pt.min() >= 0 and pt.max() <= 50.5
        assert -3.05 <= y.min() and y.max() <= 3.05


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
