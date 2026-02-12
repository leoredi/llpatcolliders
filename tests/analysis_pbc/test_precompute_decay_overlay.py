#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path
import sys

HERE = Path(__file__).resolve()
REPO_ROOT = HERE.parents[2]
ANALYSIS_ROOT = REPO_ROOT / "analysis_pbc"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(ANALYSIS_ROOT) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_ROOT))

from decay.rhn_decay_library import FLAVOUR_CONFIG
from tools.decay.precompute_decay_library_overlay import (
    _hadronized_masses,
    _seed_for_point,
    couplings_for_flavour,
)


def test_coupling_routing_electron():
    assert couplings_for_flavour("electron", 1e-6) == (1e-6, 0.0, 0.0)


def test_coupling_routing_muon():
    assert couplings_for_flavour("muon", 1e-6) == (0.0, 1e-6, 0.0)


def test_coupling_routing_tau():
    assert couplings_for_flavour("tau", 1e-6) == (0.0, 0.0, 1e-6)


def test_hadronized_mass_filter_uses_flavour_threshold():
    masses = [0.3, 0.42, 0.43, 0.53, 0.54, 2.0]

    e_threshold = float(FLAVOUR_CONFIG["electron"]["low_mass_threshold"])
    mu_threshold = float(FLAVOUR_CONFIG["muon"]["low_mass_threshold"])

    e_expected = [m for m in masses if m > e_threshold]
    mu_expected = [m for m in masses if m > mu_threshold]

    assert _hadronized_masses(masses, "electron") == e_expected
    assert _hadronized_masses(masses, "muon") == mu_expected


def test_hadronized_mass_filter_respects_overlay_min_mass():
    masses = [0.3, 0.54, 2.0, 3.5, 4.0, 4.1, 6.0]
    assert _hadronized_masses(masses, "electron", overlay_min_mass_GeV=4.0) == [4.0, 4.1, 6.0]
    assert _hadronized_masses(masses, "muon", overlay_min_mass_GeV=4.1) == [4.1, 6.0]


def test_seed_for_point_deterministic():
    seed_a = _seed_for_point(12345, "electron", 6.0)
    seed_b = _seed_for_point(12345, "electron", 6.0)
    assert seed_a == seed_b


def test_seed_for_point_varies_with_mass():
    seed_a = _seed_for_point(12345, "electron", 6.0)
    seed_b = _seed_for_point(12345, "electron", 6.1)
    assert seed_a != seed_b


def test_seed_for_point_varies_with_flavour():
    seed_e = _seed_for_point(12345, "electron", 6.0)
    seed_mu = _seed_for_point(12345, "muon", 6.0)
    seed_tau = _seed_for_point(12345, "tau", 6.0)
    assert len({seed_e, seed_mu, seed_tau}) == 3


def test_seed_for_point_positive_and_bounded():
    seed = _seed_for_point(-99999999999, "tau", 17.0)
    assert 1 <= seed <= 900_000_000
