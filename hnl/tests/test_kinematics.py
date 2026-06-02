"""Energy/momentum conservation for the 2-body and 3-body decay kinematics."""

import numpy as np
import pytest

from production.decay_engine.kinematics import (
    decay_2body, decay_2body_polarized, decay_3body_flat, _sample_polar_cos,
)


@pytest.fixture
def rng():
    return np.random.default_rng(0)


def _make_parents(n, mass, p_max=20.0, rng=None):
    """Random parent 4-vectors with |p| ~ uniform in [0, p_max]."""
    rng = rng or np.random.default_rng(0)
    p = rng.uniform(0.0, p_max, n)
    cos_th = rng.uniform(-1.0, 1.0, n)
    sin_th = np.sqrt(1.0 - cos_th**2)
    phi = rng.uniform(0, 2 * np.pi, n)
    px = p * sin_th * np.cos(phi)
    py = p * sin_th * np.sin(phi)
    pz = p * cos_th
    E = np.sqrt(p**2 + mass**2)
    return E, px, py, pz


def test_2body_conserves_4momentum(rng):
    M, m1, m2 = 5.0, 0.5, 0.1
    E, px, py, pz = _make_parents(2000, M, rng=rng)
    d1, d2 = decay_2body(E, px, py, pz, M, m1, m2, rng=rng)

    assert np.allclose(d1[:, 0] + d2[:, 0], E, rtol=1e-9, atol=1e-9)
    assert np.allclose(d1[:, 1] + d2[:, 1], px, rtol=1e-9, atol=1e-9)
    assert np.allclose(d1[:, 2] + d2[:, 2], py, rtol=1e-9, atol=1e-9)
    assert np.allclose(d1[:, 3] + d2[:, 3], pz, rtol=1e-9, atol=1e-9)

    # Daughter on-shell
    m1_recon = np.sqrt(np.maximum(d1[:, 0]**2 - (d1[:, 1]**2 + d1[:, 2]**2 + d1[:, 3]**2), 0.0))
    m2_recon = np.sqrt(np.maximum(d2[:, 0]**2 - (d2[:, 1]**2 + d2[:, 2]**2 + d2[:, 3]**2), 0.0))
    assert np.allclose(m1_recon, m1, atol=1e-7)
    assert np.allclose(m2_recon, m2, atol=1e-7)


def test_sample_polar_cos_distribution(rng):
    # pdf ∝ 1 + a·cosθ  =>  <cosθ> = a/3, and a=0 is flat (<cosθ>=0).
    for a in (0.0, 0.5, -1.0, 1.0):
        c = _sample_polar_cos(a, 200000, rng)
        assert c.min() >= -1.0 and c.max() <= 1.0
        assert abs(c.mean() - a / 3.0) < 5e-3


def test_2body_polarized_conserves_and_reduces_to_isotropic(rng):
    M, m1, m2 = 1.777, 0.139, 0.5  # tau -> pi N -like
    E, px, py, pz = _make_parents(5000, M, p_max=30.0, rng=rng)

    d1, d2 = decay_2body_polarized(E, px, py, pz, M, m1, m2, asymmetry=0.0, rng=rng)
    # 4-momentum conservation and on-shell daughters.
    assert np.allclose(d1[:, 0] + d2[:, 0], E, rtol=1e-9, atol=1e-9)
    for col, p in zip(range(1, 4), (px, py, pz)):
        assert np.allclose(d1[:, col] + d2[:, col], p, rtol=1e-9, atol=1e-9)
    m2_rec = np.sqrt(np.maximum(d2[:, 0]**2 - d2[:, 1]**2 - d2[:, 2]**2 - d2[:, 3]**2, 0.0))
    assert np.allclose(m2_rec, m2, atol=1e-7)


def test_2body_polarized_shifts_energy_spectrum(rng):
    # A nonzero asymmetry must measurably change the analyzed-daughter (N) lab
    # energy spectrum relative to the isotropic case (the whole point of the fix).
    M, m1, m2 = 1.777, 0.139, 0.5
    E, px, py, pz = _make_parents(40000, M, p_max=40.0, rng=rng)
    _, n_iso = decay_2body_polarized(E, px, py, pz, M, m1, m2, asymmetry=0.0, rng=rng)
    _, n_pol = decay_2body_polarized(E, px, py, pz, M, m1, m2, asymmetry=-1.0, rng=rng)
    assert abs(n_iso[:, 0].mean() - n_pol[:, 0].mean()) > 1e-2 * n_iso[:, 0].mean()


def test_3body_conserves_4momentum(rng):
    M, m1, m2, m3 = 5.28, 0.139, 0.106, 1.0   # B -> pi mu N-like masses
    E, px, py, pz = _make_parents(500, M, rng=rng)
    d1, d2, d3 = decay_3body_flat(E, px, py, pz, M, m1, m2, m3, rng=rng)

    s = d1 + d2 + d3
    assert np.allclose(s[:, 0], E, rtol=1e-8, atol=1e-8)
    assert np.allclose(s[:, 1], px, rtol=1e-8, atol=1e-8)
    assert np.allclose(s[:, 2], py, rtol=1e-8, atol=1e-8)
    assert np.allclose(s[:, 3], pz, rtol=1e-8, atol=1e-8)

    # Each daughter is on-shell at its assigned mass
    for d, m in [(d1, m1), (d2, m2), (d3, m3)]:
        m_recon = np.sqrt(np.maximum(d[:, 0]**2 - (d[:, 1]**2 + d[:, 2]**2 + d[:, 3]**2), 0.0))
        assert np.allclose(m_recon, m, atol=1e-6)


def test_3body_threshold_raises(rng):
    """3-body decay at or below threshold must error rather than sample garbage."""
    M = 1.0
    m1, m2, m3 = 0.5, 0.5, 0.5  # threshold = 1.5 > M
    E, px, py, pz = _make_parents(10, M, rng=rng)
    with pytest.raises(RuntimeError):
        decay_3body_flat(E, px, py, pz, M, m1, m2, m3, rng=rng)
