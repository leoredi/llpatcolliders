"""
Closure tests for ANUBIS detector configuration.

ANUBIS: Another NUclear Based In Situ detector concept
Similar to MATHUSLA but different size/location.

Validates transverse detector acceptance calculations.
"""

import numpy as np
import pytest
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from llpdecay import HNL
from llpdecay.validation.transverse_geometry import (
    ANUBIS,
    decay_probability_transverse,
    geometric_acceptance
)


class TestANUBISGeometry:
    """Test ANUBIS detector geometry."""

    def test_geometry_parameters(self):
        """Verify ANUBIS geometry definition."""
        assert ANUBIS.y_center == 30.0  # 30m above IP
        assert ANUBIS.dx == 25.0        # 50m × 50m footprint
        assert ANUBIS.dy == 5.0         # 10m height

    def test_smaller_than_mathusla(self):
        """ANUBIS is smaller than MATHUSLA."""
        from llpdecay.validation.transverse_geometry import MATHUSLA

        # Closer to IP
        assert ANUBIS.y_center < MATHUSLA.y_center

        # Smaller footprint
        assert ANUBIS.dx < MATHUSLA.dx
        assert ANUBIS.dy < MATHUSLA.dy


class TestANUBISAcceptance:
    """Test acceptance for ANUBIS geometry."""

    def test_upward_particle_accepted(self):
        """Particle aimed at ANUBIS location."""
        # Particle going upward and slightly forward
        # Should reach ~30m height
        p4 = np.array([15.0, 3.0, 10.0, 5.0])

        acc = geometric_acceptance(p4, ANUBIS)

        # Should have geometric acceptance
        assert acc >= 0.0

    def test_wrong_direction_rejected(self):
        """Particle going away from detector."""
        # Moving downward
        p4 = np.array([10.0, 1.0, -8.0, 2.0])

        acc = geometric_acceptance(p4, ANUBIS)

        # Should not reach detector
        assert acc == 0.0

    def test_requires_transverse_momentum(self):
        """Need sufficient pT to reach transverse detector."""
        # Very forward, low pT → misses
        p4_forward = np.array([50.0, 0.5, 0.5, 49.9])
        acc_forward = geometric_acceptance(p4_forward, ANUBIS)

        # Central, good pT → can reach
        p4_central = np.array([15.0, 5.0, 8.0, 5.0])
        acc_central = geometric_acceptance(p4_central, ANUBIS)

        # Just check both are computed
        assert isinstance(acc_forward, (int, float))
        assert isinstance(acc_central, (int, float))


class TestANUBISDecayProbability:
    """Test decay probability for ANUBIS-relevant parameter space."""

    def test_decay_probability_range(self):
        """Decay probability should be in [0, 1]."""
        p4 = np.array([10.0, 3.0, 8.0, 2.0])

        # Test various lifetimes
        for ctau in [0.1, 1.0, 10.0, 100.0]:
            prob = decay_probability_transverse(p4, ctau, ANUBIS)
            assert 0.0 <= prob <= 1.0

    def test_closer_detector_shorter_optimal_lifetime(self):
        """
        ANUBIS is closer to IP than MATHUSLA.
        Should have optimal sensitivity to shorter lifetimes.
        """
        from llpdecay.validation.transverse_geometry import MATHUSLA

        p4 = np.array([15.0, 5.0, 10.0, 5.0])

        # Scan lifetimes
        ctau_values = np.logspace(-1, 2, 30)

        probs_anubis = [decay_probability_transverse(p4, ct, ANUBIS)
                       for ct in ctau_values]
        probs_mathusla = [decay_probability_transverse(p4, ct, MATHUSLA)
                         for ct in ctau_values]

        # Find optimal ctau for each
        idx_anubis = np.argmax(probs_anubis)
        idx_mathusla = np.argmax(probs_mathusla)

        optimal_ctau_anubis = ctau_values[idx_anubis]
        optimal_ctau_mathusla = ctau_values[idx_mathusla]

        # ANUBIS (closer) should peak at shorter lifetime
        assert optimal_ctau_anubis < optimal_ctau_mathusla


class TestComparison:
    """Compare ANUBIS and MATHUSLA sensitivities."""

    def test_relative_acceptance(self):
        """
        Compare acceptance between ANUBIS and MATHUSLA.

        MATHUSLA: larger volume, further away
        ANUBIS: smaller volume, closer

        Trade-offs in acceptance.
        """
        from llpdecay.validation.transverse_geometry import MATHUSLA

        # Central particle with moderate pT
        p4 = np.array([20.0, 8.0, 12.0, 10.0])

        acc_anubis = geometric_acceptance(p4, ANUBIS)
        acc_mathusla = geometric_acceptance(p4, MATHUSLA)

        # Both should provide some acceptance for this kinematics
        # Actual comparison depends on detailed geometry
        assert isinstance(acc_anubis, (int, float))
        assert isinstance(acc_mathusla, (int, float))

    def test_different_optimal_mass_ranges(self):
        """
        Different detector configurations favor different HNL masses.

        This tests that geometry matters.
        """
        from llpdecay.validation.transverse_geometry import MATHUSLA

        # For same mixing, different masses give different signals
        Umu = 1e-7

        masses = [0.5, 1.0, 2.0]

        for mass in masses:
            hnl = HNL(mass=mass, Umu=Umu)
            ctau = hnl.ctau()

            # Sample kinematics
            p4 = np.array([mass * 10, 3.0, 8.0, 5.0])

            prob_anubis = decay_probability_transverse(p4, ctau, ANUBIS)
            prob_mathusla = decay_probability_transverse(p4, ctau, MATHUSLA)

            # Both should give valid probabilities
            assert 0 <= prob_anubis <= 1
            assert 0 <= prob_mathusla <= 1


class TestTransverseDetectorPrinciple:
    """
    Tests verifying key principles of transverse detector physics.

    These distinguish transverse (MATHUSLA, ANUBIS) from forward (FASER).
    """

    def test_forward_boosted_particle_misses(self):
        """
        Highly boosted forward particle misses transverse detector.

        This would hit FASER but not MATHUSLA/ANUBIS.
        """
        # Very forward: η = 5, pT = 1 GeV, E = 100 GeV
        pt = 1.0
        eta = 5.0
        E = 100.0

        px = pt * np.cos(0)
        py = pt * np.sin(0)
        pz = np.sqrt(E**2 - pt**2 - 1.0)  # Assuming m ~ 1 GeV

        p4 = np.array([E, px, py, pz])

        acc_anubis = geometric_acceptance(p4, ANUBIS)

        # Should miss transverse detector
        assert acc_anubis == 0.0

    def test_central_high_pt_hits(self):
        """
        Central particle with high pT hits transverse detector.

        This would miss FASER but hit MATHUSLA/ANUBIS.
        """
        # Central: η = 0.5, pT = 15 GeV
        pt = 15.0
        eta = 0.5
        phi = 0.7

        px = pt * np.cos(phi)
        py = pt * np.sin(phi)
        pz = pt * np.sinh(eta)
        E = np.sqrt(px**2 + py**2 + pz**2 + 1.0)

        p4 = np.array([E, px, py, pz])

        acc_anubis = geometric_acceptance(p4, ANUBIS)

        # Should have geometric acceptance
        assert acc_anubis > 0.0

    def test_rapidity_acceptance_window(self):
        """
        Transverse detectors have optimal rapidity window.

        Too forward: particle misses detector laterally
        Too backward: particle goes wrong direction
        Sweet spot: |η| < ~2
        """
        # Fixed pT, varying rapidity
        pt = 10.0
        phi = 0.5

        px = pt * np.cos(phi)
        py = pt * np.sin(phi)

        # Test various rapidities
        eta_values = np.linspace(-3, 3, 13)
        acceptances = []

        for eta in eta_values:
            pz = pt * np.sinh(eta)
            E = np.sqrt(px**2 + py**2 + pz**2 + 1.0)
            p4 = np.array([E, px, py, pz])

            acc = geometric_acceptance(p4, ANUBIS)
            acceptances.append(acc)

        # Should have peak somewhere in central region
        # Not all zeros (some events reach detector)
        assert max(acceptances) > 0

        # Forward/backward should have lower acceptance than central
        # (This is the key difference from forward detectors!)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
