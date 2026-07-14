import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from uncertainty_campaign import combine_band, stable_seed  # noqa: E402


def test_stable_seed_is_reproducible_and_member_specific():
    assert stable_seed("pdf_0001") == stable_seed("pdf_0001")
    assert stable_seed("pdf_0001") != stable_seed("pdf_0002")
    assert 0 < stable_seed("scale_muR2_muF2") < 2**32


def test_exact_campaign_combination_tracks_all_source_components():
    xc = -8.0
    variations = [("central", "central", 0.0)]
    variations += [
        (f"scale_{index}", "scale", shift)
        for index, shift in enumerate((-0.10, -0.05, -0.01, 0.01, 0.05, 0.10))
    ]
    pdf_shifts = np.linspace(-0.02, 0.02, 100)
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
    pdf_sigma = np.std(pdf_shifts, ddof=1)
    total = np.sqrt(0.10**2 + pdf_sigma**2 + 0.03**2 + 0.04**2 + 0.05**2)
    assert np.isclose(band["invf_min_total_up_dex"], total)
    assert np.isclose(band["invf_min_total_dn_dex"], total)
    assert np.isclose(band["invf_min_band_lo"], 10.0 ** (xc - total))
    assert np.isclose(band["invf_min_band_hi"], 10.0 ** (xc + total))
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
