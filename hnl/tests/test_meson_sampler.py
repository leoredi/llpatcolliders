"""Sampled meson 4-vectors should be on-shell and respect the FONLL grid bounds."""

import numpy as np

from production.constants import MESON_MASSES, FRAG_B, FRAG_C
from production.fonll.meson_sampler import sample_meson_4vectors


def test_sample_charm_on_shell_and_in_bounds():
    rng = np.random.default_rng(1)
    pool = sample_meson_4vectors(5000, "charm", rng=rng)

    E = pool['E']
    px, py, pz = pool['px'], pool['py'], pool['pz']
    m_recon = np.sqrt(np.maximum(E**2 - (px**2 + py**2 + pz**2), 0.0))
    m_exp = np.array([MESON_MASSES[int(p)] for p in pool['species_pdg']])
    # Allow tiny numerical slop from sub-bin smearing pushing pT < 0 then clipped.
    assert np.allclose(m_recon, m_exp, atol=1e-5)

    # Kinematic envelope: pt within [0, 50.5], |y| <= 3.05 after smearing.
    assert pool['pt'].min() >= 0
    assert pool['pt'].max() <= 50.5
    assert abs(pool['y']).max() <= 3.05


def test_charm_species_distribution_matches_FRAG_C():
    rng = np.random.default_rng(2)
    pool = sample_meson_4vectors(20000, "charm", rng=rng)
    species = pool['species_pdg']
    expected = np.array([FRAG_C[421], FRAG_C[411], FRAG_C[431]])
    expected = expected / expected.sum()
    counts = np.array([(species == pdg).mean() for pdg in (421, 411, 431)])
    assert np.allclose(counts, expected, atol=0.02)


def test_bottom_species_distribution_matches_FRAG_B():
    rng = np.random.default_rng(3)
    pool = sample_meson_4vectors(20000, "bottom", rng=rng)
    species = pool['species_pdg']
    expected = np.array([FRAG_B[521], FRAG_B[511], FRAG_B[531]])
    expected = expected / expected.sum()
    counts = np.array([(species == pdg).mean() for pdg in (521, 511, 531)])
    assert np.allclose(counts, expected, atol=0.02)
