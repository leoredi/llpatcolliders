import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from uncertainty_campaign import combine_band, stable_seed  # noqa: E402
from run_uncertainty_campaign import _compact_completed_run  # noqa: E402


def test_stable_seed_is_reproducible_and_member_specific():
    assert stable_seed("pdf_0001") == stable_seed("pdf_0001")
    assert stable_seed("pdf_0001") != stable_seed("pdf_0002")
    assert 0 < stable_seed("scale_muR2_muF2") < 2**32


def test_exact_campaign_builds_single_source_variation_envelope():
    xc = -8.0
    variations = [("central", "central", 0.0)]
    variations += [
        (f"scale_{index}", "scale", shift)
        for index, shift in enumerate((-0.10, -0.05, -0.01, 0.01, 0.05, 0.10))
    ]
    pdf_shifts = np.r_[np.linspace(-0.02, 0.02, 99), 1.0]
    variations += [
        (f"pdf_{index:04d}", "pdf", shift)
        for index, shift in enumerate(pdf_shifts, start=1)
    ]
    variations += [("mb_dn", "mb", -0.03), ("mb_up", "mb", 0.03)]
    variations += [
        ("gg_u", "decay_gg", -0.04),
        ("gg_d", "decay_gg", 0.0),
        ("gg_s", "decay_gg", 0.04),
        ("cbs_down", "cbs", -0.05),
        ("cbs_up", "cbs", 0.05),
    ]
    rows = []
    for name, axis, shift in variations:
        rows.append({
            "mass_GeV": 1.0,
            "variation": name,
            "axis": axis,
            "has_sensitivity": True,
            "invf_min": 10.0 ** (xc + shift),
            "invf_max": 10.0 ** (xc + 1.0 + shift),
            "invf_min_open": False,
            "invf_max_open": False,
        })
    band = combine_band(pd.DataFrame(rows)).iloc[0]
    p16, p84 = np.quantile(pdf_shifts, [0.16, 0.84])
    assert band["envelope_definition"] == "single_source_variation_envelope"
    assert np.isclose(band["invf_min_pdf_p16"], 10.0 ** (xc + p16))
    assert np.isclose(band["invf_min_pdf_p84"], 10.0 ** (xc + p84))
    assert np.isclose(
        band["invf_min_pdf_log10_std"], np.std(pdf_shifts, ddof=1)
    )
    assert np.isclose(band["invf_min_envelope_lo"], 10.0 ** (xc - 0.10))
    assert np.isclose(band["invf_min_envelope_hi"], 10.0 ** (xc + 0.10))
    assert band["invf_min_envelope_lo_source"] == "scale"
    assert band["invf_min_envelope_hi_source"] == "scale"
    assert band["invf_min_pdf_n_finite"] == 100
    assert not bool(band["invf_min_variation_missing"])


def test_missing_variation_boundary_is_explicit():
    rows = []
    definitions = [("central", "central")]
    definitions += [(f"scale_{i}", "scale") for i in range(6)]
    definitions += [(f"pdf_{i}", "pdf") for i in range(100)]
    definitions += [("mb_dn", "mb"), ("mb_up", "mb")]
    definitions += [(f"gg_{q}", "decay_gg") for q in "uds"]
    definitions += [("cbs_down", "cbs"), ("cbs_up", "cbs")]
    for name, axis in definitions:
        sensitive = name != "scale_0"
        rows.append({
            "mass_GeV": 1.0,
            "variation": name,
            "axis": axis,
            "has_sensitivity": sensitive,
            "invf_min": 1e-8 if sensitive else np.nan,
            "invf_max": 1e-7 if sensitive else np.nan,
            "invf_min_open": False,
            "invf_max_open": False,
        })
    band = combine_band(pd.DataFrame(rows)).iloc[0]
    assert bool(band["invf_min_variation_missing"])
    assert bool(band["invf_max_variation_missing"])
    assert np.isfinite(band["invf_min_envelope_lo"])


def test_central_insensitive_row_remains_an_explicit_gap():
    rows = []
    definitions = [("central", "central")]
    definitions += [(f"scale_{i}", "scale") for i in range(6)]
    definitions += [(f"pdf_{i}", "pdf") for i in range(100)]
    definitions += [("mb_dn", "mb"), ("mb_up", "mb")]
    definitions += [(f"gg_{q}", "decay_gg") for q in "uds"]
    definitions += [("cbs_down", "cbs"), ("cbs_up", "cbs")]
    for name, axis in definitions:
        sensitive = name == "scale_0"
        rows.append({
            "mass_GeV": 0.54,
            "variation": name,
            "axis": axis,
            "has_sensitivity": sensitive,
            "invf_min": 1e-8 if sensitive else np.nan,
            "invf_max": 1e-7 if sensitive else np.nan,
            "invf_min_open": False,
            "invf_max_open": False,
        })
    band = combine_band(pd.DataFrame(rows)).iloc[0]
    assert not bool(band["has_sensitivity"])
    assert bool(band["any_variation_sensitive"])
    assert np.isnan(band["invf_min_envelope_lo"])
    assert np.isnan(band["invf_max_envelope_hi"])


def test_post_validation_compaction_is_atomic_and_idempotent(tmp_path):
    run_dir = tmp_path / "runs" / "pdf_0001"
    vectors = run_dir / "llp_4vectors"
    geometry = run_dir / "analysis" / "geometry_cache"
    vectors.mkdir(parents=True)
    geometry.mkdir(parents=True)
    (vectors / "mA_1.csv").write_bytes(b"vectors")
    (geometry / "geom_1.npz").write_bytes(b"geometry")
    marker = {
        "variation": {
            "name": "pdf_0001",
            "production_mode": "fresh_600k",
        },
        "sensitivity_csv": str(run_dir / "analysis" / "bc10_sensitivity.csv"),
        "production_marker": str(run_dir / "production.complete.json"),
        "storage_state": "full",
    }
    completion = run_dir / "variation.complete.json"
    completion.write_text(json.dumps(marker))
    result = _compact_completed_run(completion, marker, run_dir, False)
    assert result["storage_state"] == "compacted"
    assert result["compaction"]["reclaimed_bytes"] == len(b"vectorsgeometry")
    assert not vectors.exists()
    assert not geometry.exists()
    assert json.loads(completion.read_text())["storage_state"] == "compacted"
    assert _compact_completed_run(completion, result, run_dir, False) == result
