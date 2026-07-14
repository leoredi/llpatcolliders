"""Focused tests for the independent BC4 uncertainty campaign."""
import hashlib
import json
import os
from pathlib import Path

import pandas as pd
import pytest

from scalar.uncertainty_band import (
    DECAY_VARIATION, _finalize_and_compact, _stable_seed,
    _validate_retained_completion, combine_band, discover_variations)


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


def test_combine_band_rebases_single_source_intervals():
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
    assert out["envelope_definition"] == "single_source_variation_envelope"
    assert bool(out["any_variation_sensitive"])
    assert out["u2_min_central"] == pytest.approx(2e-8)
    assert out["u2_min_campaign_central"] == pytest.approx(1e-8)
    assert out["u2_min_scale_envelope_lo_dex"] == pytest.approx(-0.1)
    assert out["u2_min_scale_envelope_hi_dex"] == pytest.approx(0.2)
    assert out["u2_min_pdf_p16_shift_dex"] == pytest.approx(-0.034)
    assert out["u2_min_pdf_p84_shift_dex"] == pytest.approx(0.034)
    assert out["u2_min_pdf_std_dex"] == pytest.approx(0.05 * 2**0.5)
    assert out["u2_min_bottom_mass_envelope_lo_dex"] == pytest.approx(-0.08)
    assert out["u2_min_bottom_mass_envelope_hi_dex"] == pytest.approx(0.1)
    assert out["u2_min_decay_model_shift_dex"] == pytest.approx(-0.3)
    assert out["u2_min_decay_model_alt"] == pytest.approx(2e-8 * 10**-0.3)
    assert out["u2_min_envelope_lo"] == pytest.approx(2e-8 * 10**-0.3)
    assert out["u2_min_envelope_hi"] == pytest.approx(2e-8 * 10**0.2)
    assert out["u2_min_envelope_lo_source"] == (
        f"decay_model:{DECAY_VARIATION}")
    assert out["u2_min_envelope_hi_source"] == "scale:scale_up"
    assert out["u2_max_envelope_hi"] == pytest.approx(2e-4 * 10**0.3)
    assert out["u2_max_envelope_hi_source"] == (
        f"decay_model:{DECAY_VARIATION}")


def test_pdf_replica_outlier_does_not_define_display_envelope():
    rows = [_row("central", "central", 1e-8, 1e-4)]
    rows.extend(_row(f"pdf_{index}", "pdf", 1e-8, 1e-4)
                for index in range(99))
    rows.append(_row("pdf_99", "pdf", 1e-5, 1e-1))
    rows.append(_row(DECAY_VARIATION, "decay_model", 1e-8, 1e-4))
    reference = pd.DataFrame([{
        "mass_GeV": 1.0, "has_sensitivity": True,
        "u2_min": 2e-8, "u2_max": 2e-4,
        "u2_min_open": False, "u2_max_open": False,
    }])
    out = combine_band(pd.DataFrame(rows), reference).iloc[0]
    assert out["u2_min_pdf_p84"] == pytest.approx(2e-8)
    assert out["u2_min_envelope_hi"] == pytest.approx(2e-8)
    assert out["u2_min_pdf_std_dex"] > 0


def test_completed_variation_is_checksummed_then_compacted(tmp_path):
    run_dir = tmp_path / "runs" / "central"
    vector_dir = run_dir / "llp_4vectors"
    geometry_dir = run_dir / "geometry_cache"
    result_dir = run_dir / "results"
    vector_dir.mkdir(parents=True)
    geometry_dir.mkdir()
    result_dir.mkdir()

    mass = 0.5
    config = {
        "config_sha256": "config-sha", "base_seed": 42,
        "parent_pool_seed": 123, "code_sha256": {"code.py": "abc"},
        "grid_sha256": "grid-sha", "ray_backend": {"implementation": "test"},
    }
    vector = vector_dir / "mS_0p500.csv"
    vector.write_text("1,2,3,4,5\n")
    vector_sha = hashlib.sha256(vector.read_bytes()).hexdigest()
    (vector_dir / "mS_0p500.meta.json").write_text(json.dumps({
        "config_sha256": config["config_sha256"], "mass_GeV": mass,
        "bytes": vector.stat().st_size, "sha256": vector_sha,
    }))
    (geometry_dir / "geom_mS_0p500.npz").write_bytes(b"geometry")
    (result_dir / "mS_0p500.json").write_text(json.dumps({
        "config_sha256": config["config_sha256"], "mass_GeV": mass,
        "has_sensitivity": True, "u2_min": 1e-8, "u2_max": 1e-4,
    }))

    _finalize_and_compact(
        run_dir, {"name": "central"}, config, [mass])
    assert not vector_dir.exists()
    assert not geometry_dir.exists()
    assert (run_dir / "results" / "mS_0p500.json").exists()
    assert json.loads((run_dir / "complete.json").read_text())["state"] == "compacted"
    assert _validate_retained_completion(run_dir, config, [mass]) is not None


def test_default_scratch_can_be_overridden_without_symlink(monkeypatch):
    # The module-level default is evaluated at import; the CLI option itself is
    # a real filesystem path, not a repository symlink contract.
    scratch = Path("/Volumes/GRENDEL/extra_space/bc4_uncertainty")
    assert scratch.is_absolute()
    monkeypatch.setenv("BC4_UNCERTAINTY_DIR", str(scratch))
    assert Path(os.environ["BC4_UNCERTAINTY_DIR"]) == scratch
