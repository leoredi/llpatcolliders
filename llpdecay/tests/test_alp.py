"""
Tests for ALP (Axion-Like Particle) model.

Validates decay channels, branching ratios, and kinematics.
"""

import numpy as np
import pytest
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from llpdecay import ALP
from llpdecay.core import four_vector_mass, invariant_mass_from_daughters


class TestALPInitialization:
    """Test ALP model initialization."""

    def test_valid_initialization(self):
        """Valid ALP creation."""
        alp = ALP(mass=0.5, g_agg=1e-5)
        assert alp.mass == 0.5
        assert alp.g_agg == 1e-5

    def test_photophilic_alp(self):
        """Photon-coupling dominated ALP."""
        alp = ALP(mass=0.5, g_agg=1e-4, f_a=1e10)
        assert alp.g_agg > 0
        assert all(c == 0 for c in alp.c_leptons.values())

    def test_leptophilic_alp(self):
        """Lepton-coupling dominated ALP."""
        alp = ALP(mass=1.0, g_agg=0, f_a=1e8, c_mu=1.0)
        assert alp.g_agg == 0
        assert alp.c_leptons['mu'] == 1.0

    def test_zero_couplings_warning(self):
        """Zero couplings should warn."""
        with pytest.warns(UserWarning, match="stable"):
            ALP(mass=1.0, g_agg=0, f_a=1e9)


class TestAvailableChannels:
    """Test channel accessibility."""

    def test_photophilic_channels(self):
        """Photophilic ALP has γγ channel."""
        alp = ALP(mass=0.1, g_agg=1e-5)
        channels = alp.available_channels()

        assert 'gamma_gamma' in channels

    def test_leptonic_channels(self):
        """Leptophilic ALP has lepton channels."""
        alp = ALP(mass=1.0, f_a=1e8, c_e=1.0, c_mu=1.0)
        channels = alp.available_channels()

        assert 'e_e' in channels  # a → e⁺ e⁻
        assert 'mu_mu' in channels  # a → μ⁺ μ⁻

    def test_hadronic_channels(self):
        """Hadronic ALP has meson channels."""
        alp = ALP(mass=1.0, f_a=1e8, c_u=1.0, c_d=1.0)
        channels = alp.available_channels()

        # Should have pion channel if above threshold
        if alp.mass > 2 * 0.140:
            assert 'pi_pi' in channels

    def test_threshold_behavior(self):
        """Channels below threshold shouldn't appear."""
        alp = ALP(mass=0.1, f_a=1e8, c_mu=1.0)  # Below μ⁺μ⁻ threshold
        channels = alp.available_channels()

        assert 'mu_mu' not in channels  # Should be absent


class TestBranchingRatios:
    """Test BR calculations."""

    def test_brs_sum_to_one(self):
        """BRs should be normalized."""
        alp = ALP(mass=0.5, g_agg=1e-5, f_a=1e8, c_e=1.0)
        brs = alp.branching_ratios()

        total = sum(brs.values())
        assert np.isclose(total, 1.0, rtol=1e-6)

    def test_brs_positive(self):
        """All BRs should be non-negative."""
        alp = ALP(mass=1.0, g_agg=1e-5, f_a=1e8, c_mu=1.0)
        brs = alp.branching_ratios()

        assert all(br >= 0 for br in brs.values())

    def test_dominant_channel(self):
        """Strongest coupling should dominate BRs."""
        # Very strong photon coupling
        alp = ALP(mass=0.5, g_agg=1e-3, f_a=1e10, c_e=0.01)
        brs = alp.branching_ratios()

        # γγ should dominate
        if 'gamma_gamma' in brs:
            assert brs['gamma_gamma'] > 0.9


class TestDecaySampling:
    """Test decay sampling."""

    def test_sample_single_decay(self):
        """Sample single decay."""
        alp = ALP(mass=0.5, g_agg=1e-5, seed=42)
        parent_4vec = np.array([10.0, 3.0, 0.0, 9.5])

        daughters = alp.sample_decay(parent_4vec)

        assert daughters.shape == (1, 2, 4)  # 1 event, 2 daughters, 4-vector

    def test_sample_multiple_decays(self):
        """Sample multiple decays."""
        alp = ALP(mass=0.5, g_agg=1e-5, seed=42)
        parent_4vec = np.array([10.0, 3.0, 0.0, 9.5])

        n_events = 100
        daughters = alp.sample_decay(parent_4vec, n_events=n_events)

        assert daughters.shape[0] == n_events

    def test_energy_momentum_conservation(self):
        """Decays should conserve energy-momentum."""
        alp = ALP(mass=0.5, g_agg=1e-5, seed=42)
        parent_4vec = np.array([10.0, 3.0, 0.0, 9.5])

        daughters = alp.sample_decay(parent_4vec, n_events=50)

        for i in range(len(daughters)):
            total = np.sum(daughters[i], axis=0)
            assert np.allclose(total, parent_4vec, rtol=1e-5)

    def test_mass_reconstruction(self):
        """Invariant mass should match ALP mass."""
        alp = ALP(mass=0.5, g_agg=1e-5, seed=42)
        parent_4vec = np.array([10.0, 3.0, 0.0, 9.5])

        daughters = alp.sample_decay(parent_4vec, n_events=50)
        m_recon = invariant_mass_from_daughters(daughters)

        assert np.allclose(m_recon, alp.mass, rtol=1e-5)

    def test_specific_channel(self):
        """Sample specific channel."""
        alp = ALP(mass=1.0, f_a=1e8, c_mu=1.0, seed=42)
        parent_4vec = np.array([10.0, 3.0, 0.0, 9.5])

        daughters, channels = alp.sample_decay(
            parent_4vec,
            channel='mu_mu',
            n_events=10,
            return_channel=True
        )

        assert all(ch == 'mu_mu' for ch in channels)

    def test_isotropic_decay(self):
        """Scalar ALP should decay isotropically."""
        alp = ALP(mass=0.5, g_agg=1e-5, seed=42)
        parent_4vec = np.array([2.0, 0.0, 0.0, 0.0])  # At rest

        daughters = alp.sample_decay(parent_4vec, n_events=1000)

        # Extract angles of daughter 1
        p1 = daughters[:, 0, 1:4]
        p_mag = np.linalg.norm(p1, axis=1)
        cos_theta = p1[:, 2] / p_mag

        # Should be uniform in cos θ
        mean = np.mean(cos_theta)
        std = np.std(cos_theta)

        assert np.abs(mean) < 0.1  # Mean ≈ 0
        assert np.abs(std - 1.0/np.sqrt(3)) < 0.1  # Std ≈ 0.577


class TestLifetimeCalculations:
    """Test lifetime and width."""

    def test_total_width_positive(self):
        """Width should be positive."""
        alp = ALP(mass=0.5, g_agg=1e-5)
        width = alp.total_width()

        assert width > 0

    def test_width_scales_with_coupling(self):
        """Width should scale with coupling squared."""
        alp1 = ALP(mass=0.5, g_agg=1e-5)
        alp2 = ALP(mass=0.5, g_agg=2e-5)

        width1 = alp1.total_width()
        width2 = alp2.total_width()

        # Γ ∝ g²
        ratio = width2 / width1
        assert np.isclose(ratio, 4.0, rtol=0.1)

    def test_ctau_positive(self):
        """Decay length should be positive."""
        alp = ALP(mass=0.5, g_agg=1e-5)
        ctau = alp.ctau()

        assert ctau > 0

    def test_ctau_decreases_with_coupling(self):
        """Larger coupling → shorter lifetime."""
        alp_weak = ALP(mass=0.5, g_agg=1e-6)
        alp_strong = ALP(mass=0.5, g_agg=1e-4)

        ctau_weak = alp_weak.ctau()
        ctau_strong = alp_strong.ctau()

        assert ctau_weak > ctau_strong


class TestDaughterPDGs:
    """Test PDG ID retrieval."""

    def test_photon_pdgs(self):
        """γγ channel PDG IDs."""
        alp = ALP(mass=0.5, g_agg=1e-5)
        pdgs = alp.get_daughter_pdgs('gamma_gamma')

        assert pdgs == [22, 22]

    def test_lepton_pdgs(self):
        """Lepton pair PDG IDs."""
        alp = ALP(mass=1.0, f_a=1e8, c_mu=1.0)
        pdgs = alp.get_daughter_pdgs('mu_mu')

        assert abs(pdgs[0]) == 13
        assert abs(pdgs[1]) == 13
        assert pdgs[0] * pdgs[1] < 0  # Opposite charges

    def test_charged_count(self):
        """Test charged particle counting."""
        alp = ALP(mass=1.0, f_a=1e8, c_mu=1.0)

        # γγ has 0 charged
        assert alp.get_charged_count('gamma_gamma') == 0

        # μ⁺μ⁻ has 2 charged
        assert alp.get_charged_count('mu_mu') == 2


class TestComparison:
    """Compare ALP with HNL."""

    def test_different_angular_distributions(self):
        """ALP (scalar) vs HNL (fermion) angular distributions."""
        # ALP is scalar → isotropic
        # HNL can be polarized → anisotropic

        from llpdecay import HNL

        alp = ALP(mass=0.5, g_agg=1e-5, seed=42)
        hnl = HNL(mass=0.5, Umu=1e-6, seed=43)

        # Both at rest
        parent_rest = np.array([0.5, 0.0, 0.0, 0.0])

        # This is just a structural test - both should work
        daughters_alp = alp.sample_decay(parent_rest, n_events=10)
        daughters_hnl = hnl.sample_decay(parent_rest, n_events=10)

        assert daughters_alp.shape[2] == 4
        assert daughters_hnl.shape[2] == 4


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
