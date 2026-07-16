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


def test_published_bundle_hashes_match_manifest():
    manifest = json.loads((BUNDLE / "MANIFEST.json").read_text())
    assert manifest["central_curve"]["sha256"] == _sha256(
        PUBLISHED / "grendel_hnl_sensitivity.csv")
    for name, metadata in manifest["files"].items():
        path = BUNDLE / name
        assert path.exists(), name
        assert metadata["sha256"] == _sha256(path), name


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


def test_canonical_manifest_names_actual_beta_fix_commits():
    manifest = json.loads((PUBLISHED / "MANIFEST.json").read_text())
    assert manifest["producer_git_sha"] == "fad3a8b"
    assert manifest["diagnostic_bundle"]["campaign_git_sha"].startswith(
        "6482484")
