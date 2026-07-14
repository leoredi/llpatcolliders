"""Focused tests for the independent BC4 uncertainty campaign."""
import os
from pathlib import Path

import pandas as pd
import pytest

from scalar.uncertainty_band import (
    DECAY_VARIATION, _stable_seed, combine_band, discover_variations)


def _workspace_grid_dir():
    repo = Path(__file__).resolve().parents[2]
    return repo.parents[2] / "shared" / "NNPDF40" / "fonll-local" / "output"


def test_complete_fonll_manifest_is_discovered():
    _, variations = discover_variations(_workspace_grid_dir())
    counts = pd.Series([variation["axis"] for variation in variations]).value_counts()
    assert counts.to_dict() == {
        "pdf": 100, "scale": 6, "mass": 2,
        "central": 1, "decay_model": 1,
    }
    assert variations[0]["name"] == "central"
    assert variations[1]["name"] == DECAY_VARIATION
    assert variations[0]["path"] == variations[1]["path"]
    assert variations[0]["width_scheme"] != variations[1]["width_scheme"]


def test_variation_seeds_are_deterministic_and_independent():
    central = _stable_seed(42, "central", "parent_pool")
    assert central == _stable_seed(42, "central", "parent_pool")
    assert central != _stable_seed(42, DECAY_VARIATION, "parent_pool")
    assert central != _stable_seed(42, "central", "0p500", "production")
    assert 0 <= central < 2**32


def _row(name, axis, umin, umax):
    return {
        "mass_GeV": 1.0, "variation": name, "axis": axis,
        "has_sensitivity": True, "u2_min": umin, "u2_max": umax,
        "u2_min_open": False, "u2_max_open": False,
    }


def test_combine_band_rebases_asymmetric_components():
    raw = pd.DataFrame([
        _row("central", "central", 1e-8, 1e-4),
        _row("scale_up", "scale", 10**-7.8, 10**-3.8),
        _row("scale_down", "scale", 10**-8.1, 10**-4.1),
        _row("pdf_1", "pdf", 10**-7.95, 10**-3.95),
        _row("pdf_2", "pdf", 10**-8.05, 10**-4.05),
        _row("mb_dn", "mass", 10**-7.9, 10**-3.9),
        _row("mb_up", "mass", 10**-8.08, 10**-4.08),
        _row(DECAY_VARIATION, "decay_model", 10**-8.3, 10**-3.7),
    ])
    reference = pd.DataFrame([{
        "mass_GeV": 1.0, "has_sensitivity": True,
        "u2_min": 2e-8, "u2_max": 2e-4,
        "u2_min_open": False, "u2_max_open": False,
    }])
    out = combine_band(raw, reference).iloc[0]
    assert out["u2_min_central"] == pytest.approx(2e-8)
    assert out["u2_min_campaign_central"] == pytest.approx(1e-8)
    assert out["u2_min_scale_up_dex"] == pytest.approx(0.2)
    assert out["u2_min_scale_dn_dex"] == pytest.approx(0.1)
    assert out["u2_min_pdf_sigma_dex"] == pytest.approx(0.05 * 2**0.5)
    assert out["u2_min_mb_dev_dex"] == pytest.approx(0.1)
    assert out["u2_min_decay_model_dn_dex"] == pytest.approx(0.3)
    assert out["u2_min_decay_model_alt"] == pytest.approx(2e-8 * 10**-0.3)
    assert out["u2_min_band_lo"] < out["u2_min_fonll_band_lo"]
    assert out["u2_max_band_hi"] > out["u2_max_fonll_band_hi"]


def test_default_scratch_can_be_overridden_without_symlink(monkeypatch):
    # The module-level default is evaluated at import; the CLI option itself is
    # a real filesystem path, not a repository symlink contract.
    scratch = Path("/Volumes/GRENDEL/extra_space/bc4_uncertainty")
    assert scratch.is_absolute()
    monkeypatch.setenv("BC4_UNCERTAINTY_DIR", str(scratch))
    assert Path(os.environ["BC4_UNCERTAINTY_DIR"]) == scratch
