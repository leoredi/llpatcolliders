"""
Tests for core kinematics functions.

Validates Lorentz transformations, two-body decay kinematics,
and phase space calculations.
"""

import numpy as np
import pytest
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from llpdecay.core import (
    kallen,
    four_vector_mass,
    boost_to_lab,
    two_body_decay_momenta,
    sample_two_body_decay,
    invariant_mass_from_daughters,
    transverse_momentum,
)


class TestKallenFunction:
    """Test Källén triangle function."""

    def test_kallen_symmetric(self):
        """Källén function should be symmetric in its arguments."""
        a, b, c = 1.0, 2.0, 3.0
        assert np.isclose(kallen(a, b, c), kallen(b, c, a))
        assert np.isclose(kallen(a, b, c), kallen(c, a, b))

    def test_kallen_zero(self):
        """Källén function should be zero on threshold."""
        # λ(a, b, c) = 0 when a = (√b + √c)²
        a = (np.sqrt(1.0) + np.sqrt(4.0))**2  # (1 + 2)² = 9
        b, c = 1.0, 4.0
        assert np.isclose(kallen(a, b, c), 0.0, atol=1e-10)


class TestFourVectorMass:
    """Test invariant mass calculation."""

    def test_massive_particle_at_rest(self):
        """Particle at rest: m² = E²."""
        p = np.array([1.0, 0.0, 0.0, 0.0])
        assert np.isclose(four_vector_mass(p), 1.0)

    def test_massless_particle(self):
        """Photon: m² = 0."""
        p = np.array([1.0, 0.0, 0.0, 1.0])  # E = pz
        assert np.isclose(four_vector_mass(p), 0.0, atol=1e-10)

    def test_moving_particle(self):
        """Particle in motion: m² = E² - p²."""
        E = 5.0
        p_mag = 4.0
        p = np.array([E, p_mag, 0.0, 0.0])
        expected_mass = np.sqrt(E**2 - p_mag**2)
        assert np.isclose(four_vector_mass(p), expected_mass)

    def test_batch(self):
        """Test batch processing."""
        p = np.array([
            [1.0, 0.0, 0.0, 0.0],
            [5.0, 4.0, 0.0, 0.0],
        ])
        masses = four_vector_mass(p)
        assert masses.shape == (2,)
        assert np.isclose(masses[0], 1.0)
        assert np.isclose(masses[1], 3.0)


class TestTwoBodyDecayMomenta:
    """Test two-body decay kinematics."""

    def test_equal_mass_daughters(self):
        """Decay to two equal-mass daughters."""
        m_parent = 2.0
        m1 = m2 = 0.5

        p_mag, E1, E2 = two_body_decay_momenta(m_parent, m1, m2)

        # Equal masses → equal energies
        assert np.isclose(E1, E2)
        assert np.isclose(E1 + E2, m_parent)

        # Check on-shell: E² = p² + m²
        assert np.isclose(E1**2, p_mag**2 + m1**2)

    def test_pion_decay(self):
        """π⁺ → μ⁺ ν decay."""
        m_pi = 0.13957
        m_mu = 0.105658
        m_nu = 0.0

        p_mag, E_mu, E_nu = two_body_decay_momenta(m_pi, m_mu, m_nu)

        # Energy conservation
        assert np.isclose(E_mu + E_nu, m_pi)

        # Neutrino massless → E = p
        assert np.isclose(E_nu, p_mag)

        # Muon on-shell
        assert np.isclose(E_mu**2, p_mag**2 + m_mu**2)

    def test_forbidden_decay(self):
        """Kinematically forbidden decay should raise error."""
        m_parent = 1.0
        m1 = m2 = 1.0  # m1 + m2 > m_parent

        with pytest.raises(ValueError, match="Kinematically forbidden"):
            two_body_decay_momenta(m_parent, m1, m2)


class TestSampleTwoBodyDecay:
    """Test two-body decay sampling."""

    def test_energy_momentum_conservation(self):
        """Sampled decays should conserve energy-momentum."""
        m_parent = 1.0
        m1, m2 = 0.3, 0.2
        n_events = 100

        daughters = sample_two_body_decay(m_parent, m1, m2, n_events=n_events)

        assert daughters.shape == (n_events, 2, 4)

        for i in range(n_events):
            p1, p2 = daughters[i]

            # Energy conservation
            E_total = p1[0] + p2[0]
            assert np.isclose(E_total, m_parent, rtol=1e-6)

            # Momentum conservation (back-to-back in rest frame)
            p_total = p1[1:4] + p2[1:4]
            assert np.allclose(p_total, 0.0, atol=1e-10)

            # On-shell conditions
            m1_calc = four_vector_mass(p1)
            m2_calc = four_vector_mass(p2)
            assert np.isclose(m1_calc, m1, rtol=1e-6)
            assert np.isclose(m2_calc, m2, rtol=1e-6)

    def test_isotropic_distribution(self):
        """Unpolarized decay should be isotropic."""
        m_parent = 1.0
        m1, m2 = 0.1, 0.1
        n_events = 10000

        daughters = sample_two_body_decay(
            m_parent, m1, m2,
            n_events=n_events,
            polarization=0.0,  # Isotropic
            seed=42
        )

        # Extract cos(θ) from daughter 1
        p1 = daughters[:, 0, 1:4]
        p_mag = np.linalg.norm(p1, axis=1)
        cos_theta = p1[:, 2] / p_mag

        # Should be uniform in [-1, 1]
        # Check mean ≈ 0 and std ≈ 1/√3
        assert np.abs(np.mean(cos_theta)) < 0.05
        assert np.abs(np.std(cos_theta) - 1.0/np.sqrt(3)) < 0.05

    def test_polarized_distribution(self):
        """Polarized decay should have asymmetry."""
        m_parent = 1.0
        m1, m2 = 0.1, 0.1
        n_events = 10000

        daughters = sample_two_body_decay(
            m_parent, m1, m2,
            n_events=n_events,
            polarization=1.0,  # Fully polarized
            alpha=1.0,  # Maximal asymmetry
            seed=42
        )

        # Extract cos(θ)
        p1 = daughters[:, 0, 1:4]
        p_mag = np.linalg.norm(p1, axis=1)
        cos_theta = p1[:, 2] / p_mag

        # Should prefer forward direction (cos θ > 0)
        assert np.mean(cos_theta) > 0.2  # Significantly positive


class TestBoostToLab:
    """Test Lorentz boost transformations."""

    def test_boost_identity(self):
        """Boost by particle at rest is identity."""
        p_rest = np.array([[1.0, 0.5, 0.0, 0.866]])
        parent_at_rest = np.array([1.0, 0.0, 0.0, 0.0])

        p_lab = boost_to_lab(p_rest, parent_at_rest)

        assert np.allclose(p_lab, p_rest)

    def test_boost_along_z(self):
        """Boost along z-axis."""
        # Daughter at rest in parent frame
        m_daughter = 0.5
        p_rest = np.array([[m_daughter, 0.0, 0.0, 0.0]])

        # Parent moving along z
        m_parent = 1.0
        E_parent = 2.0
        pz_parent = np.sqrt(E_parent**2 - m_parent**2)
        parent = np.array([E_parent, 0.0, 0.0, pz_parent])

        p_lab = boost_to_lab(p_rest, parent)

        # Check that daughter now has pz
        assert p_lab[0, 3] > 0  # Should be moving in +z

        # Daughter should still be on-shell
        m_check = four_vector_mass(p_lab[0])
        assert np.isclose(m_check, m_daughter, rtol=1e-6)

    def test_decay_reconstruction(self):
        """Boosted decay products should reconstruct parent."""
        m_parent = 2.0
        m1, m2 = 0.5, 0.3

        # Sample decay in rest frame
        daughters_rest = sample_two_body_decay(m_parent, m1, m2, n_events=1, seed=42)

        # Parent moving in lab
        parent_4vec = np.array([10.0, 3.0, 2.0, 9.0])

        # Boost daughters to lab
        daughters_lab = boost_to_lab(daughters_rest[0], parent_4vec)

        # Sum should give parent 4-vector
        total_4vec = np.sum(daughters_lab, axis=0)
        assert np.allclose(total_4vec, parent_4vec, rtol=1e-6)


class TestInvariantMass:
    """Test invariant mass from daughters."""

    def test_two_body_reconstruction(self):
        """Reconstructed mass should match parent."""
        m_parent = 2.0
        m1, m2 = 0.5, 0.3

        daughters = sample_two_body_decay(m_parent, m1, m2, n_events=10)
        m_recon = invariant_mass_from_daughters(daughters)

        assert np.allclose(m_recon, m_parent, rtol=1e-6)


class TestTransverseMomentum:
    """Test pT calculation."""

    def test_pt_zero_for_beam_axis(self):
        """Particle along z has pT = 0."""
        p = np.array([5.0, 0.0, 0.0, 5.0])
        pt = transverse_momentum(p)
        assert np.isclose(pt, 0.0, atol=1e-10)

    def test_pt_calculation(self):
        """pT = sqrt(px² + py²)."""
        p = np.array([5.0, 3.0, 4.0, 1.0])
        pt = transverse_momentum(p)
        assert np.isclose(pt, 5.0)  # sqrt(9 + 16)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
