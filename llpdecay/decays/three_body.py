"""
Three-body decay phase space sampling.

Implements proper 3-body kinematics using Dalitz plot variables and
importance sampling with matrix element weights.

References:
- PDG Review of Kinematics, Section 49.2
- James, F. (1968) Monte Carlo Phase Space, CERN Yellow Report 68-15
- Byckling & Kajantie (1973) Particle Kinematics
"""

import numpy as np
from typing import Tuple, Optional, Callable
import warnings

from ..core import kallen, four_vector_mass, create_rng


class ThreeBodyPhaseSpace:
    """
    Three-body decay phase space generator.

    For decay M → 1 + 2 + 3, generates events distributed according to
    phase space and optional matrix element weight.

    Parameters
    ----------
    m_parent : float
        Parent mass (GeV)
    m1, m2, m3 : float
        Daughter masses (GeV)
    matrix_element : callable, optional
        Function ME(s12, s13, s23, m_parent, m1, m2, m3) -> weight
        If None, use uniform phase space

    Attributes
    ----------
    m_parent, m1, m2, m3 : float
        Particle masses
    s_min, s_max : dict
        Kinematic limits for Dalitz variables
    """

    def __init__(
        self,
        m_parent: float,
        m1: float,
        m2: float,
        m3: float,
        matrix_element: Optional[Callable] = None
    ):
        self.m_parent = m_parent
        self.m1 = m1
        self.m2 = m2
        self.m3 = m3
        self.matrix_element = matrix_element

        # Check kinematic accessibility
        if m_parent < m1 + m2 + m3:
            raise ValueError(
                f"Kinematically forbidden: m_parent={m_parent:.4f} < "
                f"m1+m2+m3={m1+m2+m3:.4f}"
            )

        # Calculate Dalitz plot boundaries
        self._compute_dalitz_limits()

    def _compute_dalitz_limits(self):
        """
        Calculate kinematic boundaries in Dalitz plot.

        For invariant masses s_ij = (p_i + p_j)², the allowed region is
        determined by energy-momentum conservation.
        """
        M = self.m_parent
        m1, m2, m3 = self.m1, self.m2, self.m3

        # s_12 ranges
        s12_min = (m1 + m2)**2
        s12_max = (M - m3)**2

        # s_13 ranges
        s13_min = (m1 + m3)**2
        s13_max = (M - m2)**2

        # s_23 ranges
        s23_min = (m2 + m3)**2
        s23_max = (M - m1)**2

        self.s_min = {
            '12': s12_min,
            '13': s13_min,
            '23': s23_min,
        }
        self.s_max = {
            '12': s12_max,
            '13': s13_max,
            '23': s23_max,
        }

    def s13_limits(self, s12: float) -> Tuple[float, float]:
        """
        Calculate s_13 limits for given s_12.

        The Dalitz plot boundary is determined by requiring all energies
        in the 12 rest frame to be positive.

        Parameters
        ----------
        s12 : float
            Invariant mass squared of particles 1 and 2

        Returns
        -------
        s13_min, s13_max : float
            Allowed range of s_13 for this s_12
        """
        M = self.m_parent
        m1, m2, m3 = self.m1, self.m2, self.m3

        # Energy of particle 3 in 12 rest frame
        E3_12 = (M**2 - s12 - m3**2) / (2 * np.sqrt(s12))

        # Momentum of particle 3 in 12 rest frame
        p3_12_sq = E3_12**2 - m3**2
        if p3_12_sq < 0:
            return None, None
        p3_12 = np.sqrt(p3_12_sq)

        # Energy of particle 1 in 12 rest frame
        E1_12 = (s12 + m1**2 - m2**2) / (2 * np.sqrt(s12))

        # Momentum of particle 1 in 12 rest frame
        p1_12_sq = E1_12**2 - m1**2
        if p1_12_sq < 0:
            return None, None
        p1_12 = np.sqrt(p1_12_sq)

        # s_13 = m1² + m3² + 2(E1*E3 - p1*p3*cosθ) in 12 frame
        # Extrema occur at cosθ = ±1
        s13_min = m1**2 + m3**2 + 2 * (E1_12 * E3_12 - p1_12 * p3_12)
        s13_max = m1**2 + m3**2 + 2 * (E1_12 * E3_12 + p1_12 * p3_12)

        return s13_min, s13_max

    def sample(
        self,
        n_events: int = 1,
        rng: Optional[np.random.Generator] = None
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Sample three-body decays in parent rest frame.

        Uses importance sampling with matrix element weights if provided,
        otherwise samples uniform phase space.

        Parameters
        ----------
        n_events : int
            Number of events to generate
        rng : np.random.Generator, optional
            Random number generator

        Returns
        -------
        four_vectors : np.ndarray, shape (n_events, 3, 4)
            4-vectors of daughters [E, px, py, pz] in parent rest frame
        weights : np.ndarray, shape (n_events,)
            Event weights (1.0 for uniform sampling, ME value otherwise)

        Notes
        -----
        The algorithm:
        1. Sample s_12 uniformly in [s_12_min, s_12_max]
        2. For each s_12, sample s_13 in allowed range
        3. Sample azimuthal angle φ uniformly
        4. Construct 4-vectors from Dalitz variables
        5. Apply matrix element weight if provided
        """
        if rng is None:
            rng = create_rng()

        four_vectors = []
        weights = []

        # Maximum weight for rejection sampling (if using ME)
        max_weight = 1.0
        if self.matrix_element is not None:
            # Estimate max by sampling
            max_weight = self._estimate_max_weight(rng)

        n_generated = 0
        attempts = 0
        max_attempts = n_events * 1000  # Safety limit

        while n_generated < n_events and attempts < max_attempts:
            attempts += 1

            # Sample Dalitz variables
            s12 = rng.uniform(self.s_min['12'], self.s_max['12'])

            # Get s13 limits for this s12
            s13_min, s13_max = self.s13_limits(s12)
            if s13_min is None:
                continue

            s13 = rng.uniform(s13_min, s13_max)

            # Calculate s23 from constraint
            s23 = self.m_parent**2 + self.m1**2 + self.m2**2 + self.m3**2 - s12 - s13

            # Sample azimuthal angle
            phi = rng.uniform(0, 2 * np.pi)

            # Construct 4-vectors
            try:
                p1, p2, p3 = self._construct_four_vectors(s12, s13, phi)
            except ValueError:
                continue

            # Calculate weight
            if self.matrix_element is not None:
                weight = self.matrix_element(s12, s13, s23, self.m_parent,
                                            self.m1, self.m2, self.m3)
                # Rejection sampling
                if rng.uniform(0, max_weight) > weight:
                    continue
            else:
                weight = 1.0

            four_vectors.append(np.array([p1, p2, p3]))
            weights.append(weight)
            n_generated += 1

        if n_generated < n_events:
            warnings.warn(
                f"Only generated {n_generated}/{n_events} events after "
                f"{attempts} attempts. Adjust ME or check kinematics."
            )

        return np.array(four_vectors), np.array(weights)

    def _construct_four_vectors(
        self,
        s12: float,
        s13: float,
        phi: float
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Construct daughter 4-vectors from Dalitz variables.

        Strategy:
        1. Build particles 1,2 in their rest frame (12 frame)
        2. Build particle 3 in 12 frame
        3. Boost everything to parent rest frame

        Parameters
        ----------
        s12, s13 : float
            Dalitz plot invariant masses
        phi : float
            Azimuthal angle

        Returns
        -------
        p1, p2, p3 : np.ndarray, shape (4,)
            Daughter 4-vectors in parent rest frame
        """
        M = self.m_parent
        m1, m2, m3 = self.m1, self.m2, self.m3

        # --- Build in 12 rest frame ---
        m12 = np.sqrt(s12)

        # Energy and momentum of 1,2 in their rest frame
        E1_12 = (s12 + m1**2 - m2**2) / (2 * m12)
        E2_12 = (s12 + m2**2 - m1**2) / (2 * m12)

        p_12_sq = kallen(s12, m1**2, m2**2) / (4 * s12)
        if p_12_sq < 0:
            raise ValueError("Negative momentum squared in 12 frame")
        p_12 = np.sqrt(p_12_sq)

        # Place 1 along +z in 12 frame
        p1_12 = np.array([E1_12, 0, 0, p_12])
        p2_12 = np.array([E2_12, 0, 0, -p_12])

        # Particle 3 in 12 frame
        E3_12 = (M**2 - s12 - m3**2) / (2 * m12)
        p3_12_mag_sq = E3_12**2 - m3**2
        if p3_12_mag_sq < 0:
            raise ValueError("Negative momentum squared for particle 3")
        p3_12_mag = np.sqrt(p3_12_mag_sq)

        # Direction of particle 3 from s13
        # s_13 = m1² + m3² + 2(E1*E3 - p1*p3*cosθ)
        cos_theta_13 = (s13 - m1**2 - m3**2 - 2*E1_12*E3_12) / (-2*p_12*p3_12_mag)

        # Clamp to [-1, 1] for numerical stability
        cos_theta_13 = np.clip(cos_theta_13, -1.0, 1.0)
        sin_theta_13 = np.sqrt(1 - cos_theta_13**2)

        # 3-momentum of particle 3 in 12 frame
        p3_12_x = p3_12_mag * sin_theta_13 * np.cos(phi)
        p3_12_y = p3_12_mag * sin_theta_13 * np.sin(phi)
        p3_12_z = p3_12_mag * cos_theta_13

        p3_12 = np.array([E3_12, p3_12_x, p3_12_y, p3_12_z])

        # --- Boost to parent rest frame ---
        # The 12 system has 4-momentum (E_12, p_12_vec) in parent frame
        E_12 = (M**2 + s12 - m3**2) / (2 * M)
        p_12_vec_mag_sq = E_12**2 - s12
        if p_12_vec_mag_sq < 0:
            raise ValueError("Negative 12 momentum in parent frame")
        p_12_vec_mag = np.sqrt(p_12_vec_mag_sq)

        # 12 system moving along -z in parent frame (convention)
        p_12_parent = np.array([E_12, 0, 0, -p_12_vec_mag])

        # Boost
        from ..core import boost_to_lab
        p1 = boost_to_lab(p1_12, p_12_parent)
        p2 = boost_to_lab(p2_12, p_12_parent)
        p3 = boost_to_lab(p3_12, p_12_parent)

        return p1, p2, p3

    def _estimate_max_weight(
        self,
        rng: np.random.Generator,
        n_sample: int = 10000
    ) -> float:
        """
        Estimate maximum matrix element for rejection sampling.

        Samples random points in Dalitz plot and finds maximum weight.

        Parameters
        ----------
        rng : np.random.Generator
            Random generator
        n_sample : int
            Number of points to sample

        Returns
        -------
        float
            Estimated maximum weight (with 10% safety margin)
        """
        max_weight = 0.0

        for _ in range(n_sample):
            s12 = rng.uniform(self.s_min['12'], self.s_max['12'])
            s13_min, s13_max = self.s13_limits(s12)

            if s13_min is None:
                continue

            s13 = rng.uniform(s13_min, s13_max)
            s23 = self.m_parent**2 + self.m1**2 + self.m2**2 + self.m3**2 - s12 - s13

            weight = self.matrix_element(s12, s13, s23, self.m_parent,
                                        self.m1, self.m2, self.m3)
            max_weight = max(max_weight, weight)

        return max_weight * 1.1  # 10% safety margin


def sample_three_body_decay(
    m_parent: float,
    m1: float,
    m2: float,
    m3: float,
    n_events: int = 1,
    matrix_element: Optional[Callable] = None,
    rng: Optional[np.random.Generator] = None
) -> np.ndarray:
    """
    Sample three-body decay M → 1 + 2 + 3 in parent rest frame.

    Convenience function wrapping ThreeBodyPhaseSpace.

    Parameters
    ----------
    m_parent : float
        Parent mass (GeV)
    m1, m2, m3 : float
        Daughter masses (GeV)
    n_events : int
        Number of events to generate
    matrix_element : callable, optional
        Matrix element function(s12, s13, s23, M, m1, m2, m3) -> weight
    rng : np.random.Generator, optional
        Random generator

    Returns
    -------
    np.ndarray, shape (n_events, 3, 4)
        Daughter 4-vectors in parent rest frame

    Examples
    --------
    >>> # Uniform phase space
    >>> daughters = sample_three_body_decay(2.0, 0.1, 0.1, 0.0, n_events=100)

    >>> # With matrix element
    >>> def ME(s12, s13, s23, M, m1, m2, m3):
    ...     return (s12 - m1**2 - m2**2) * (s13 - m1**2 - m3**2)
    >>> daughters = sample_three_body_decay(2.0, 0.1, 0.1, 0.0,
    ...                                      n_events=100, matrix_element=ME)
    """
    ps = ThreeBodyPhaseSpace(m_parent, m1, m2, m3, matrix_element=matrix_element)
    four_vectors, weights = ps.sample(n_events, rng)
    return four_vectors


# =============================================================================
# Matrix elements for specific 3-body decays
# =============================================================================

def hnl_three_body_leptonic_me(
    s12: float,
    s13: float,
    s23: float,
    m_N: float,
    m_nu: float,
    m_l1: float,
    m_l2: float
) -> float:
    """
    Matrix element squared for HNL → ν ℓ⁺ ℓ⁻.

    This decay proceeds via virtual W* or Z* exchange. The full calculation
    is complex, so we use a simplified approximation valid in the massless
    neutrino limit.

    From Gorbunov & Shaposhnikov, arXiv:hep-ph/9911364, the differential
    decay rate for N → ν e⁺ e⁻ is:

    dΓ ∝ G_F² |U|² × f(s_ee, s_νe⁺, s_νe⁻)

    where f includes propagator and spinor structure.

    Parameters
    ----------
    s12 : float
        Invariant mass² of ν and ℓ⁺ (m_nu=m1, m_l1=m2)
    s13 : float
        Invariant mass² of ν and ℓ⁻ (m_nu=m1, m_l2=m3)
    s23 : float
        Invariant mass² of ℓ⁺ and ℓ⁻
    m_N : float
        HNL mass
    m_nu : float
        Neutrino mass (≈0)
    m_l1, m_l2 : float
        Charged lepton masses

    Returns
    -------
    float
        Matrix element squared (proportional to decay rate)

    References
    ----------
    - Gorbunov & Shaposhnikov, Phys.Rev.D75:083010, arXiv:hep-ph/9911364
    - Atre et al., JHEP 0905:030, arXiv:0901.3589, Eq. (43)
    """
    from ..core import G_FERMI

    # Simplified matrix element (leading order approximation)
    # Full ME includes W/Z propagators and helicity structure

    # V-A structure gives preference to certain regions of phase space
    # Approximate as product of propagators

    # W propagator for charged current contribution
    # Effective squared propagator: 1/(s_ij - m_W²)²
    # In low-mass limit s_ij << m_W², this becomes constant

    # For simplicity, use a phenomenological form that captures the
    # essential features:
    # - Suppression near lepton pair threshold (helicity)
    # - Enhancement from W/Z resonances (if m_N large enough)

    # Phase space factors
    x_ee = s23 / m_N**2  # ℓ⁺ℓ⁻ invariant mass fraction

    # Helicity suppression near threshold
    threshold_factor = x_ee * (1 - x_ee)

    # For massless leptons, the ME is symmetric in the neutrino direction
    # Full ME: |M|² ∝ s_νe⁺ * s_νe⁻ / (product of propagators)

    # Simplified: assume constant matrix element modulated by threshold
    me_squared = threshold_factor

    # Include mass corrections for charged leptons
    if m_l1 > 0 or m_l2 > 0:
        # Helicity suppression for massive leptons
        # (1 - m_l²/s_parent²)² factors
        beta_l1_sq = 1 - m_l1**2 / m_N**2
        beta_l2_sq = 1 - m_l2**2 / m_N**2
        me_squared *= beta_l1_sq * beta_l2_sq

    # Normalize to reasonable scale
    # This is a relative weight; absolute normalization doesn't matter
    # for phase space sampling
    return me_squared
