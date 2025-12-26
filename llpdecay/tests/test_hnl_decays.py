"""
Tests for HNL model implementation.

Validates decay channel selection, branching ratios, and decay sampling.
"""

import numpy as np
import pytest
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from llpdecay import HNL
from llpdecay.core import four_vector_mass, invariant_mass_from_daughters


class TestHNLInitialization:
    """Test HNL model initialization."""

    def test_valid_initialization(self):
        """Valid HNL should initialize correctly."""
        hnl = HNL(mass=2.0, Umu=1e-6)
        assert hnl.mass == 2.0
        assert hnl.Umu == 1e-6
        assert hnl.Ue == 0.0
        assert hnl.Utau == 0.0

    def test_invalid_mass(self):
        """Negative mass should raise error."""
        with pytest.raises(ValueError, match="positive"):
            HNL(mass=-1.0, Umu=1e-6)

    def test_zero_mixing(self):
        """Zero mixing should issue warning."""
        with pytest.warns(UserWarning, match="stable"):
            HNL(mass=1.0, Ue=0, Umu=0, Utau=0)

    def test_flavor_fractions(self):
        """Flavor fractions should sum to 1."""
        hnl = HNL(mass=2.0, Ue=1e-6, Umu=2e-6, Utau=3e-6)
        fractions = hnl.flavor_fractions

        assert np.isclose(sum(fractions.values()), 1.0)
        assert np.isclose(fractions['e'], 1.0/6.0)
        assert np.isclose(fractions['mu'], 2.0/6.0)
        assert np.isclose(fractions['tau'], 3.0/6.0)


class TestAvailableChannels:
    """Test kinematic channel accessibility."""

    def test_low_mass_channels(self):
        """Low mass HNL has limited channels."""
        hnl = HNL(mass=0.2, Umu=1e-6)  # Below pion threshold
        channels = hnl.available_channels()

        # Should have only leptonic channels (3-body)
        assert len(channels) > 0
        # Should NOT have pion channels
        assert 'mu_pi' not in channels
        assert 'e_pi' not in channels

    def test_intermediate_mass_channels(self):
        """Intermediate mass HNL has pion channels."""
        hnl = HNL(mass=1.0, Umu=1e-6)
        channels = hnl.available_channels()

        # Should have charged pion channels
        assert 'mu_pi' in channels
        assert 'e_pi' in channels

        # Should NOT have kaon channels yet (threshold ~0.6 GeV)
        # Actually K+ mass is 0.494 GeV, so at 1 GeV it should be available
        assert 'mu_K' in channels

    def test_high_mass_channels(self):
        """High mass HNL has many channels."""
        hnl = HNL(mass=3.0, Umu=1e-6)
        channels = hnl.available_channels()

        # Should have most channels
        assert 'mu_pi' in channels
        assert 'mu_K' in channels
        assert 'mu_rho' in channels


class TestBranchingRatios:
    """Test branching ratio calculations."""

    def test_brs_sum_to_one(self):
        """Branching ratios should be normalized."""
        hnl = HNL(mass=2.0, Umu=1e-6)
        brs = hnl.branching_ratios()

        total = sum(brs.values())
        assert np.isclose(total, 1.0, rtol=1e-6)

    def test_brs_positive(self):
        """All branching ratios should be positive."""
        hnl = HNL(mass=2.0, Umu=1e-6)
        brs = hnl.branching_ratios()

        assert all(br >= 0 for br in brs.values())

    def test_flavor_dependence(self):
        """BRs should depend on mixing flavor."""
        # Electron-only mixing
        hnl_e = HNL(mass=2.0, Ue=1e-6, Umu=0, Utau=0)
        brs_e = hnl_e.branching_ratios()

        # Muon-only mixing
        hnl_mu = HNL(mass=2.0, Ue=0, Umu=1e-6, Utau=0)
        brs_mu = hnl_mu.branching_ratios()

        # Electron channels should dominate for Ue mixing
        if 'e_pi' in brs_e and 'mu_pi' in brs_e:
            # e_pi should be larger for electron mixing
            assert brs_e.get('e_pi', 0) > brs_e.get('mu_pi', 0)

        # Muon channels should dominate for Umu mixing
        if 'mu_pi' in brs_mu and 'e_pi' in brs_mu:
            assert brs_mu.get('mu_pi', 0) > brs_mu.get('e_pi', 0)

    def test_majorana_vs_dirac(self):
        """Majorana should have larger total width (factor ~2)."""
        mass = 2.0
        Umu = 1e-6

        hnl_majorana = HNL(mass=mass, Umu=Umu, is_majorana=True)
        hnl_dirac = HNL(mass=mass, Umu=Umu, is_majorana=False)

        width_majorana = hnl_majorana.total_width()
        width_dirac = hnl_dirac.total_width()

        # Majorana should decay faster (larger width)
        assert width_majorana > width_dirac
        # Approximately factor of 2
        assert np.isclose(width_majorana / width_dirac, 2.0, rtol=0.5)


class TestDecaySampling:
    """Test decay sampling functionality."""

    def test_sample_single_decay(self):
        """Sample single decay."""
        hnl = HNL(mass=2.0, Umu=1e-6, seed=42)
        parent_4vec = np.array([10.0, 3.0, 0.0, 9.5])

        daughters = hnl.sample_decay(parent_4vec)

        # Should return 2 daughters (for 2-body decay)
        assert daughters.shape[0] == 1  # 1 event
        assert daughters.shape[1] == 2  # 2 daughters
        assert daughters.shape[2] == 4  # 4-vectors

    def test_sample_multiple_decays(self):
        """Sample multiple decays."""
        hnl = HNL(mass=2.0, Umu=1e-6, seed=42)
        parent_4vec = np.array([10.0, 3.0, 0.0, 9.5])

        n_events = 100
        daughters = hnl.sample_decay(parent_4vec, n_events=n_events)

        assert daughters.shape[0] == n_events

    def test_energy_momentum_conservation(self):
        """Decay should conserve energy-momentum."""
        hnl = HNL(mass=2.0, Umu=1e-6, seed=42)
        parent_4vec = np.array([10.0, 3.0, 0.0, 9.5])

        n_events = 50
        daughters = hnl.sample_decay(parent_4vec, n_events=n_events)

        for i in range(n_events):
            # Sum daughter 4-vectors
            total_4vec = np.sum(daughters[i], axis=0)

            # Should match parent
            assert np.allclose(total_4vec, parent_4vec, rtol=1e-5)

    def test_mass_reconstruction(self):
        """Invariant mass should match HNL mass."""
        hnl = HNL(mass=2.0, Umu=1e-6, seed=42)
        parent_4vec = np.array([10.0, 3.0, 0.0, 9.5])

        n_events = 50
        daughters = hnl.sample_decay(parent_4vec, n_events=n_events)

        m_recon = invariant_mass_from_daughters(daughters)

        assert np.allclose(m_recon, hnl.mass, rtol=1e-5)

    def test_specific_channel(self):
        """Sample specific decay channel."""
        hnl = HNL(mass=2.0, Umu=1e-6, seed=42)
        parent_4vec = np.array([10.0, 3.0, 0.0, 9.5])

        daughters, channels = hnl.sample_decay(
            parent_4vec,
            channel='mu_pi',
            n_events=10,
            return_channel=True
        )

        # All channels should be mu_pi
        assert all(ch == 'mu_pi' for ch in channels)

    def test_channel_distribution(self):
        """Sampled channels should follow BRs."""
        hnl = HNL(mass=2.0, Umu=1e-6, seed=42)
        parent_4vec = np.array([10.0, 3.0, 0.0, 9.5])

        n_events = 10000
        _, channels = hnl.sample_decay(
            parent_4vec,
            n_events=n_events,
            return_channel=True
        )

        # Count channels
        channel_counts = {}
        for ch in channels:
            channel_counts[ch] = channel_counts.get(ch, 0) + 1

        # Get expected BRs
        brs = hnl.branching_ratios()

        # Check that sampled frequencies match BRs
        for ch, count in channel_counts.items():
            freq = count / n_events
            expected = brs[ch]
            # Should be within ~3 sigma (for large n_events)
            sigma = np.sqrt(expected * (1 - expected) / n_events)
            assert np.abs(freq - expected) < 5 * sigma, \
                f"Channel {ch}: freq={freq:.3f}, expected={expected:.3f}"

    def test_batch_parents(self):
        """Sample decays for multiple parent 4-vectors."""
        hnl = HNL(mass=2.0, Umu=1e-6, seed=42)

        # Different parent momenta
        parents = np.array([
            [10.0, 3.0, 0.0, 9.5],
            [15.0, 5.0, 2.0, 13.0],
            [20.0, 10.0, 5.0, 15.0],
        ])

        daughters = hnl.sample_decay(parents)

        assert daughters.shape[0] == 3  # 3 events

        # Check conservation for each
        for i in range(3):
            total = np.sum(daughters[i], axis=0)
            assert np.allclose(total, parents[i], rtol=1e-5)


class TestLifetimeCalculations:
    """Test lifetime and decay length calculations."""

    def test_total_width_positive(self):
        """Total width should be positive."""
        hnl = HNL(mass=2.0, Umu=1e-6)
        width = hnl.total_width()

        assert width > 0

    def test_width_scales_with_mixing(self):
        """Width should scale with |U|²."""
        mass = 2.0

        hnl1 = HNL(mass=mass, Umu=1e-6)
        hnl2 = HNL(mass=mass, Umu=4e-6)

        width1 = hnl1.total_width()
        width2 = hnl2.total_width()

        # Width ∝ |U|², so width2 / width1 ≈ 4
        ratio = width2 / width1
        assert np.isclose(ratio, 4.0, rtol=0.1)

    def test_ctau_positive(self):
        """Decay length should be positive."""
        hnl = HNL(mass=2.0, Umu=1e-6)
        ctau = hnl.ctau()

        assert ctau > 0

    def test_ctau_decreases_with_mixing(self):
        """Larger mixing → shorter lifetime."""
        mass = 2.0

        hnl_weak = HNL(mass=mass, Umu=1e-7)
        hnl_strong = HNL(mass=mass, Umu=1e-5)

        ctau_weak = hnl_weak.ctau()
        ctau_strong = hnl_strong.ctau()

        # Stronger mixing → shorter decay length
        assert ctau_weak > ctau_strong

    def test_decay_probability(self):
        """Test decay probability calculation."""
        hnl = HNL(mass=2.0, Umu=1e-6)
        parent_4vec = np.array([10.0, 3.0, 0.0, 9.5])

        # Probability to decay between 10m and 50m
        P = hnl.decay_probability(parent_4vec, x_min=10.0, x_max=50.0)

        # Should be between 0 and 1
        assert 0 <= P <= 1


class TestDaughterPDGs:
    """Test PDG ID retrieval."""

    def test_mu_pi_pdgs(self):
        """mu_pi channel should have correct PDG IDs."""
        hnl = HNL(mass=2.0, Umu=1e-6)
        pdgs = hnl.get_daughter_pdgs('mu_pi')

        assert len(pdgs) == 2
        # Should be muon (±13) and pion (±211)
        assert abs(pdgs[0]) == 13
        assert abs(pdgs[1]) == 211

    def test_charged_count(self):
        """Test charged particle counting."""
        hnl = HNL(mass=2.0, Umu=1e-6)

        # mu_pi has 2 charged particles
        assert hnl.get_charged_count('mu_pi') == 2

        # nu_pi0 has 0 charged (invisible)
        assert hnl.get_charged_count('nu_mu_pi0') == 0


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
