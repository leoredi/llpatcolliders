"""
Tests for three-body decay phase space.

Validates Dalitz plot sampling, energy-momentum conservation,
and matrix element weighting.
"""

import numpy as np
import pytest
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from llpdecay.decays.three_body import (
    ThreeBodyPhaseSpace,
    sample_three_body_decay,
    hnl_three_body_leptonic_me
)
from llpdecay.core import four_vector_mass, invariant_mass_from_daughters


class TestThreeBodyPhaseSpace:
    """Test 3-body phase space generator."""

    def test_initialization(self):
        """Valid initialization."""
        ps = ThreeBodyPhaseSpace(2.0, 0.1, 0.1, 0.0)
        assert ps.m_parent == 2.0
        assert ps.m1 == 0.1
        assert ps.m2 == 0.1
        assert ps.m3 == 0.0

    def test_kinematic_limits(self):
        """Dalitz plot limits should be correctly calculated."""
        ps = ThreeBodyPhaseSpace(2.0, 0.1, 0.1, 0.0)

        # s12 range
        assert ps.s_min['12'] == pytest.approx((0.1 + 0.1)**2)
        assert ps.s_max['12'] == pytest.approx((2.0 - 0.0)**2)

    def test_forbidden_decay_raises_error(self):
        """Kinematically forbidden should raise error."""
        with pytest.raises(ValueError, match="Kinematically forbidden"):
            ThreeBodyPhaseSpace(1.0, 0.5, 0.5, 0.5)

    def test_uniform_phase_space_sampling(self):
        """Sample uniform phase space."""
        ps = ThreeBodyPhaseSpace(2.0, 0.1, 0.1, 0.0)
        daughters, weights = ps.sample(n_events=50)

        assert daughters.shape == (50, 3, 4)
        assert np.all(weights == 1.0)  # Uniform → all weights = 1

    def test_energy_momentum_conservation(self):
        """3-body decays should conserve energy-momentum."""
        ps = ThreeBodyPhaseSpace(2.0, 0.1, 0.1, 0.0)
        daughters, _ = ps.sample(n_events=100, rng=np.random.default_rng(42))

        for i in range(len(daughters)):
            # Sum 4-vectors
            total = np.sum(daughters[i], axis=0)

            # Should be parent at rest
            expected = np.array([2.0, 0.0, 0.0, 0.0])
            assert np.allclose(total, expected, rtol=1e-5)

    def test_mass_reconstruction(self):
        """Invariant mass should match parent."""
        ps = ThreeBodyPhaseSpace(2.0, 0.1, 0.1, 0.0)
        daughters, _ = ps.sample(n_events=50, rng=np.random.default_rng(42))

        m_recon = invariant_mass_from_daughters(daughters)
        assert np.allclose(m_recon, 2.0, rtol=1e-5)

    def test_daughter_masses_on_shell(self):
        """Each daughter should have correct mass."""
        m_parent = 2.0
        m1, m2, m3 = 0.105, 0.105, 0.0  # μ⁺ μ⁻ ν

        ps = ThreeBodyPhaseSpace(m_parent, m1, m2, m3)
        daughters, _ = ps.sample(n_events=50, rng=np.random.default_rng(42))

        for i in range(len(daughters)):
            mass1 = four_vector_mass(daughters[i, 0])
            mass2 = four_vector_mass(daughters[i, 1])
            mass3 = four_vector_mass(daughters[i, 2])

            assert np.isclose(mass1, m1, rtol=1e-4)
            assert np.isclose(mass2, m2, rtol=1e-4)
            assert np.isclose(mass3, m3, atol=1e-4)

    def test_dalitz_plot_coverage(self):
        """Events should cover allowed Dalitz plot region."""
        ps = ThreeBodyPhaseSpace(2.0, 0.1, 0.1, 0.0)
        daughters, _ = ps.sample(n_events=1000, rng=np.random.default_rng(42))

        # Calculate Dalitz variables
        s12_values = []
        s13_values = []

        for i in range(len(daughters)):
            p1, p2, p3 = daughters[i]

            # s12 = (p1 + p2)²
            p12 = p1 + p2
            s12 = four_vector_mass(p12)**2

            # s13 = (p1 + p3)²
            p13 = p1 + p3
            s13 = four_vector_mass(p13)**2

            s12_values.append(s12)
            s13_values.append(s13)

        s12_values = np.array(s12_values)
        s13_values = np.array(s13_values)

        # Check that values are within bounds
        assert np.all(s12_values >= ps.s_min['12'])
        assert np.all(s12_values <= ps.s_max['12'])
        assert np.all(s13_values >= ps.s_min['13'])
        assert np.all(s13_values <= ps.s_max['13'])

    def test_with_matrix_element(self):
        """Matrix element weighting should affect distribution."""

        # Simple ME: prefer high s12
        def me_high_s12(s12, s13, s23, M, m1, m2, m3):
            return s12  # Linear in s12

        ps_uniform = ThreeBodyPhaseSpace(2.0, 0.1, 0.1, 0.0, matrix_element=None)
        ps_weighted = ThreeBodyPhaseSpace(2.0, 0.1, 0.1, 0.0, matrix_element=me_high_s12)

        daughters_uniform, _ = ps_uniform.sample(n_events=500, rng=np.random.default_rng(42))
        daughters_weighted, _ = ps_weighted.sample(n_events=500, rng=np.random.default_rng(43))

        # Calculate mean s12 for each
        def get_s12(daughters):
            s12_vals = []
            for i in range(len(daughters)):
                p12 = daughters[i, 0] + daughters[i, 1]
                s12_vals.append(four_vector_mass(p12)**2)
            return np.mean(s12_vals)

        s12_uniform = get_s12(daughters_uniform)
        s12_weighted = get_s12(daughters_weighted)

        # Weighted should have higher mean s12
        assert s12_weighted > s12_uniform


class TestConvenienceFunction:
    """Test sample_three_body_decay wrapper."""

    def test_basic_usage(self):
        """Basic sampling."""
        daughters = sample_three_body_decay(2.0, 0.1, 0.1, 0.0, n_events=10)

        assert daughters.shape == (10, 3, 4)

    def test_with_matrix_element(self):
        """Sampling with custom ME."""

        def me(s12, s13, s23, M, m1, m2, m3):
            return s12 * s13

        daughters = sample_three_body_decay(
            2.0, 0.1, 0.1, 0.0,
            n_events=10,
            matrix_element=me,
            rng=np.random.default_rng(42)
        )

        assert daughters.shape == (10, 3, 4)


class TestHNLMatrixElement:
    """Test HNL 3-body leptonic matrix element."""

    def test_me_positive(self):
        """Matrix element should be positive."""
        s12 = 1.0
        s13 = 1.0
        s23 = 2.0
        m_N = 2.0
        m_nu = 0.0
        m_l1 = m_l2 = 0.105

        me = hnl_three_body_leptonic_me(s12, s13, s23, m_N, m_nu, m_l1, m_l2)

        assert me >= 0

    def test_me_threshold_suppression(self):
        """ME should be suppressed near lepton pair threshold."""
        m_N = 1.0
        m_nu = 0.0
        m_l = 0.105

        # Near threshold
        s23_threshold = (2 * m_l)**2 + 0.01
        s12 = 0.5
        s13 = m_N**2 - s12 - s23_threshold + m_nu**2 + m_l**2 + m_l**2

        me_threshold = hnl_three_body_leptonic_me(s12, s13, s23_threshold,
                                                   m_N, m_nu, m_l, m_l)

        # Well above threshold
        s23_high = 0.8
        s13_high = m_N**2 - s12 - s23_high + m_nu**2 + m_l**2 + m_l**2

        me_high = hnl_three_body_leptonic_me(s12, s13_high, s23_high,
                                              m_N, m_nu, m_l, m_l)

        # Should be larger away from threshold
        assert me_high > me_threshold


class TestIntegration:
    """Integration tests with full decay chain."""

    def test_hnl_three_body_decay(self):
        """Complete 3-body decay sampling for HNL."""
        from llpdecay import HNL

        hnl = HNL(mass=2.0, Umu=1e-6, seed=42)

        parent_4vec = np.array([10.0, 3.0, 0.0, 9.5])

        # Sample 3-body channel if available
        channels = hnl.available_channels()
        three_body_channels = [ch for ch in channels
                               if 'nu_' in ch and ch not in ['nu_e_pi0', 'nu_mu_pi0', 'nu_tau_pi0']]

        if three_body_channels:
            ch = three_body_channels[0]
            daughters, sampled_ch = hnl.sample_decay(
                parent_4vec,
                channel=ch,
                n_events=10,
                return_channel=True
            )

            assert daughters.shape[0] == 10
            assert daughters.shape[1] == 3  # 3 daughters
            assert all(c == ch for c in sampled_ch)

            # Check conservation
            for i in range(10):
                total = np.sum(daughters[i], axis=0)
                assert np.allclose(total, parent_4vec, rtol=1e-5)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
