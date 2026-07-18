"""Integrity and topology checks for the tracked HNL diagnostic bundle."""

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


PUBLISHED = Path(__file__).resolve().parents[1] / "data" / "published"
BUNDLE = PUBLISHED / "bundle"


def _sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def test_published_bundle_internal_hashes_match():
    manifest = json.loads((BUNDLE / "MANIFEST.json").read_text())
    for name, metadata in manifest["files"].items():
        path = BUNDLE / name
        assert path.exists(), name
        assert metadata["sha256"] == _sha256(path), name


def test_bundle_is_beta_fix_baseline_superseded_by_kaon_rerun():
    # The bundle is the 2026-07-07 beta=p/E baseline. The 2026-07-18 kaon rerun
    # superseded the published central at low-mass Ue/Umu WITHOUT regenerating the
    # bundle, so bundle central_curve.sha256 records the baseline it was built against
    # (echoed by the top MANIFEST) and no longer equals the live central curve.
    bundle = json.loads((BUNDLE / "MANIFEST.json").read_text())
    top = json.loads((PUBLISHED / "MANIFEST.json").read_text())
    baseline = bundle["central_curve"]["sha256"]
    assert top["diagnostic_bundle"]["baseline_central_sha256"] == baseline
    live = _sha256(PUBLISHED / "grendel_hnl_sensitivity.csv")
    assert top["csv_sha256"] == live      # top manifest tracks the live central
    assert live != baseline               # which the kaon rerun moved at low mass
    assert "kaon" in top["diagnostic_bundle"]["kaon_rerun_supersession"].lower()


def test_fonll_bundle_has_every_member_and_stable_topology():
    band = pd.read_csv(BUNDLE / "hnl_band_fonll.csv")
    raw = pd.read_csv(BUNDLE / "hnl_band_fonll_raw.csv")

    assert len(band) == 54
    assert len(raw) == 111 * 54
    assert raw["variation_name"].nunique() == 111
    assert set(raw.groupby("variation_name").size()) == {54}
    for boundary in ("u2_min", "u2_max"):
        assert not band[f"{boundary}_topology_changed"].any()
        assert set(band[f"{boundary}_n_members_expected"]) == {110}
        assert set(band[f"{boundary}_n_members_sensitive"]) == {110}


def test_channel_fractions_and_endpoint_controls_are_consistent():
    channels = pd.read_csv(BUNDLE / "channel_breakdown_u2min.csv")
    grouped = channels.groupby(["flavor", "mass_GeV", "boundary"])
    assert set(grouped["channel"].nunique()) == {8}
    assert np.allclose(grouped["frac"].sum(), 1.0, rtol=0.0, atol=1e-12)

    controls = json.loads((BUNDLE / "NUMERICAL_CONTROLS.json").read_text())
    assert controls["worker_count_independence"][
        "single_vs_two_worker_bitwise_equal"]
    assert controls["bc_endpoint_exact400"]["all_topologies_agree"]


def test_canonical_manifest_provenance():
    manifest = json.loads((PUBLISHED / "MANIFEST.json").read_text())
    # current publication is the 2026-07-18 kaon rerun
    assert manifest["producer_git_sha"] == "80668e6"
    # the earlier beta=p/E fix is preserved as prior_correction
    assert "fad3a8b" in manifest["prior_correction"]
    # the FONLL diagnostic campaign is unchanged
    assert manifest["diagnostic_bundle"]["campaign_git_sha"].startswith(
        "6482484")
