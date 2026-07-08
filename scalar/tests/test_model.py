"""Sanity checks for the BC4 scalar-portal model layer (scalar/model.py).

Validate against the published Winkler arXiv:1809.01876 / BC4 numbers, not
against our own re-derivation.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scalar import model


def test_higgs_vev():
    # v = (sqrt(2) G_F)^(-1/2) ~ 246 GeV.
    assert model.V_HIGGS == pytest.approx(246.2, abs=0.5)


def test_leptonic_width_threshold_and_scaling():
    # Below 2 m_mu the muon channel is closed.
    assert model.width_leptonic(0.2, model.M_MUON) == 0.0
    # Open above threshold and positive.
    assert model.width_leptonic(0.5, model.M_MUON) > 0.0
    # Linear in sin^2 theta is handled by the scan; the width here is at s2t=1.


def test_branching_ratios_normalize():
    for m_S in (0.25, 0.5, 1.0, 2.5, 4.0):
        br = model.branching_ratios(m_S)
        assert sum(br.values()) == pytest.approx(1.0, abs=1e-9)


def test_low_mass_is_dimuon():
    # Just above 2 m_mu and below 2 m_pi: essentially pure mu mu (ee is tiny).
    br = model.branching_ratios(0.25)
    assert br["mumu"] > 0.99


def test_pipi_dominates_mid_mass():
    # 0.5-0.9 GeV: pi pi is the largest channel (Winkler Fig. 4).
    br = model.branching_ratios(0.7)
    assert br["pipi"] > br["mumu"]
    assert br["pipi"] > 0.5


def test_spectator_region_opens_gg_cc_tautau():
    br = model.branching_ratios(4.0)
    # tau tau open above 2 m_tau ~ 3.55 GeV; c c above 2 m_D ~ 3.73 GeV.
    assert br["tautau"] > 0.0
    assert br["cc"] > 0.0
    assert br["pipi"] == 0.0   # dispersive region closed above 2 GeV


def test_production_br_matches_literature():
    # Exclusive BR(B+ -> K+ S)/sin^2 theta ~ 0.5 near threshold (Winkler eqs.
    # A3/A5/A6); inclusive BR(B -> X_s S) ~ 5-6 (eq. A7).
    assert model.br_B_to_K_S(0.5, "B+") == pytest.approx(0.48, abs=0.08)
    assert model.br_B_to_Xs_S(0.5, "B+") == pytest.approx(5.3, abs=0.8)


def test_production_linear_in_coupling():
    base = model.br_B_to_K_S(1.0, "B+", sin2theta=1.0)
    scaled = model.br_B_to_K_S(1.0, "B+", sin2theta=1e-8)
    assert scaled == pytest.approx(base * 1e-8, rel=1e-6)


def test_production_closes_above_kinematic_limit():
    assert float(model.br_B_to_K_S(4.9, "B+")) == 0.0
    assert model.M_S_MAX_BTOK == pytest.approx(4.785, abs=0.01)


def test_ctau_inverse_in_coupling():
    # c*tau ~ 1 / sin^2 theta.
    c1 = model.ctau(1.0, 1e-8)
    c2 = model.ctau(1.0, 1e-9)
    assert c2 == pytest.approx(10.0 * c1, rel=1e-6)
    # Reaches the ~metre GRENDEL scale in the sensitive coupling window.
    assert 1e-3 < model.ctau(1.0, 1e-8) < 1e2


def test_K_to_pi_S_only_below_threshold():
    assert model.br_K_to_pi_S(0.2) > 0.0
    assert float(model.br_K_to_pi_S(0.4)) == 0.0   # above m_K - m_pi ~ 0.354
