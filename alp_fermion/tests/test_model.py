"""Physics sanity checks for the BC10 ALP model layer."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import model as m


def test_lepton_width_threshold_and_scaling():
    # closed below 2 m_mu, open above
    assert m.alp_partial_widths(0.20, 1e-3)["mumu"] == 0.0
    assert m.alp_partial_widths(0.25, 1e-3)["mumu"] > 0.0
    # Gamma ~ (1/f)^2
    g1 = m.alp_total_width(1.0, 1e-3)
    g2 = m.alp_total_width(1.0, 2e-3)
    assert g2 == pytest.approx(4.0 * g1, rel=1e-6)


def test_ctau_scaling_inverse_square():
    c1 = m.alp_ctau(1.0, 1e-6)
    c2 = m.alp_ctau(1.0, 1e-7)
    assert c2 / c1 == pytest.approx(100.0, rel=1e-6)


def test_branchings_sum_to_one():
    for ma in (0.25, 0.5, 1.0, 3.0, 4.5):
        br = m.alp_branchings(ma, 1e-3)
        assert sum(br.values()) == pytest.approx(1.0, rel=1e-9)


def test_tautau_opens_above_threshold():
    assert m.alp_branchings(3.0, 1e-3)["tautau"] == 0.0
    assert m.alp_branchings(4.0, 1e-3)["tautau"] > 0.0


def test_BtoKa_threshold_and_coupling_scaling():
    # closed above the kinematic edge m_B - m_K
    assert m.br_B_to_K_a(5.0, 1e-3, "B+") == 0.0
    assert m.br_B_to_K_a(1.0, 1e-3, "B+") > 0.0
    # BR ~ (1/f)^2
    b1 = m.br_B_to_K_a(1.0, 1e-4, "B+")
    b2 = m.br_B_to_K_a(1.0, 2e-4, "B+")
    assert b2 == pytest.approx(4.0 * b1, rel=1e-6)


def test_BtoKa_order_of_magnitude():
    # leading-log estimate at 1/f = 1e-3 GeV^-1 (f = 1 TeV) is O(0.1).
    assert 0.05 < m.br_B_to_K_a(1.0, 1e-3, "B+") < 0.5
