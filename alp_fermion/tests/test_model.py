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


def test_mumu_br_anchor_at_1GeV():
    # GKOZ (arXiv:2310.03524): BR(a->mumu) "< 10% for m_a >~ 1 GeV".  With the
    # data-driven hadronic table in our conventions it lands at ~9%.
    br = m.alp_branchings(1.0, 1e-3)["mumu"]
    assert 0.05 < br < 0.10


def test_hadronic_table_seam_is_continuous():
    # the matched perturbative continuation equals the table at the ceiling
    mt = m.table_mass_max()
    w_tab = m.alp_partial_widths(mt - 1e-9, 1e-3)["hadronic"]
    w_cont = m.alp_partial_widths(mt + 1e-9, 1e-3)["hadronic"]
    assert w_cont == pytest.approx(w_tab, rel=0.02)


def test_visible_fraction_bounds():
    for ma in (0.25, 0.7, 1.0, 1.5, 2.5, 3.5):
        v = m.visible_fraction(ma)
        assert 0.0 < v <= 1.0
    # around 1--1.5 GeV the all-neutral modes (3pi0 etc.) are a real deficit
    assert m.visible_fraction(1.5) < 0.9


def test_visible_weights_subset_of_branchings():
    br = m.alp_branchings(1.0, 1e-3)
    vis = m.visible_channel_weights(1.0, 1e-3)
    assert vis["hadronic"] <= br["hadronic"] + 1e-12
    assert vis["mumu"] == pytest.approx(br["mumu"], rel=1e-9)
    assert "gammagamma" not in vis


def test_BtoKa_threshold_and_coupling_scaling():
    # closed above the kinematic edge m_B - m_K
    assert m.br_B_to_K_a(5.0, 1e-3, "B+") == 0.0
    assert m.br_B_to_K_a(1.0, 1e-3, "B+") > 0.0
    # BR ~ (1/f)^2
    b1 = m.br_B_to_K_a(1.0, 1e-4, "B+")
    b2 = m.br_B_to_K_a(1.0, 2e-4, "B+")
    assert b2 == pytest.approx(4.0 * b1, rel=1e-6)


def test_BtoKa_normalization_anchor():
    # g_bs = CBS_EFF/f with the Boiarska et al. form factors: BR(B+ -> K+ a)
    # ~ 2.4e-2 at 1/f = 1e-3 (checked against GKOZ Table 1 |C_bs| to ~20%).
    assert 0.01 < m.br_B_to_K_a(1.0, 1e-3, "B+") < 0.05


def test_kaon_tower_matches_GKOZ_factor():
    # GKOZ: full tower ~ 4x the K + K*(892) rate
    ma = 1.0
    k = m.br_B_to_Ka(ma, 1e-3, "B+", "K")
    kstar = m.br_B_to_Ka(ma, 1e-3, "B+", "Kstar_892")
    tower = m.br_B_tower_total(ma, 1e-3, "B+")
    assert 2.5 < tower / (k + kstar) < 5.5


def test_tower_channels_close_at_thresholds():
    # heavy resonances switch off as m_a approaches their kinematic edge
    assert m.br_B_to_Ka(4.0, 1e-3, "B+", "Kstar_1680") == 0.0
    assert m.br_B_to_Ka(1.0, 1e-3, "B+", "Kstar_1680") > 0.0
