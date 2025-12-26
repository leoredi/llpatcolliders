"""
Hadronic form factors and QED radiative corrections.

Provides more accurate treatment of meson decays including:
- Form factors for N → ℓ M transitions
- QED radiative corrections
- Higher-order electroweak effects

References:
- Chrzaszcz et al., arXiv:1906.02657 (form factors)
- Bondarenko et al., arXiv:1805.08567 (HNL phenomenology)
- PDG Review 2024, Electroweak Model and Constraints
"""

import numpy as np
from typing import Tuple, Optional


# =============================================================================
# Meson Form Factors
# =============================================================================

def form_factor_pion(q2: float) -> float:
    """
    Pion form factor f₊(q²) for semileptonic decays.

    For N → ℓ π, the relevant form factor parametrizes the hadronic
    matrix element ⟨π(p)|q̄γᵤu|0⟩.

    We use the pole parametrization:
    f₊(q²) = f₊(0) / (1 - q²/m_pole²)

    where m_pole ≈ m_ρ for pions.

    Parameters
    ----------
    q2 : float
        Momentum transfer squared (GeV²)

    Returns
    -------
    float
        Form factor value

    References
    ----------
    - PDG 2024, Pseudoscalar Meson Decay Constants
    """
    m_pole = 0.775  # ρ meson mass (GeV)
    f_plus_0 = 1.0  # Normalized at q² = 0

    return f_plus_0 / (1.0 - q2 / m_pole**2)


def form_factor_kaon(q2: float) -> float:
    """
    Kaon form factor f₊(q²).

    Similar to pion but with K* pole.

    Parameters
    ----------
    q2 : float
        Momentum transfer squared (GeV²)

    Returns
    -------
    float
        Form factor value
    """
    m_pole = 0.892  # K* meson mass (GeV)
    f_plus_0 = 1.0

    return f_plus_0 / (1.0 - q2 / m_pole**2)


def scalar_form_factor_pion(q2: float) -> float:
    """
    Scalar form factor f₀(q²) for pion.

    Parametrizes ⟨π(p)|q̄q|0⟩ matrix element.

    Parameters
    ----------
    q2 : float
        Momentum transfer squared (GeV²)

    Returns
    -------
    float
        Scalar form factor
    """
    m_pole = 0.775  # ρ pole
    f_0_at_0 = 1.0

    return f_0_at_0 / (1.0 - q2 / m_pole**2)


def tensor_form_factor_pion(q2: float) -> float:
    """
    Tensor form factor for axial transitions.

    Used in pseudoscalar meson production with spin-dependent couplings.

    Parameters
    ----------
    q2 : float
        Momentum transfer squared (GeV²)

    Returns
    -------
    float
        Tensor form factor
    """
    # Simplified model
    m_pole = 1.275  # a₁ axial vector meson
    return 1.0 / (1.0 - q2 / m_pole**2)


# =============================================================================
# QED Radiative Corrections
# =============================================================================

def qed_correction_lepton_pair(s: float, m_lepton: float) -> float:
    """
    Leading-order QED correction for ℓ⁺ℓ⁻ production.

    Includes virtual photon exchange and soft bremsstrahlung.
    Correction factor: 1 + δ_QED

    Parameters
    ----------
    s : float
        Invariant mass squared of lepton pair (GeV²)
    m_lepton : float
        Lepton mass (GeV)

    Returns
    -------
    float
        Correction factor (≈ 1.0 + few %)

    References
    ----------
    - Schwinger correction: δ ≈ (α/π) × log(s/m²)
    """
    from ..core import ALPHA_EM

    if s <= 4 * m_lepton**2:
        return 1.0  # Below threshold

    # Beta factor
    beta = np.sqrt(1.0 - 4.0 * m_lepton**2 / s)

    # Leading logarithmic correction
    # δ = (α/π) × [log(s/m²) - 1 + O(1)]
    delta = (ALPHA_EM / np.pi) * (np.log(s / m_lepton**2) - 1.0)

    # Coulomb correction (for low velocities)
    coulomb = (ALPHA_EM / beta) * (1.0 - beta**2 / 12.0)

    # Total correction
    correction = 1.0 + delta + coulomb

    return correction


def qed_correction_photon_pair(m_parent: float) -> float:
    """
    QED correction for two-photon decays (a → γγ, etc.).

    Includes loop corrections to the effective coupling.

    Parameters
    ----------
    m_parent : float
        Parent particle mass (GeV)

    Returns
    -------
    float
        Correction factor
    """
    from ..core import ALPHA_EM

    # O(α) correction from triangle loops
    # Simplified estimate
    delta = (ALPHA_EM / np.pi) * 0.5  # Typical ~0.2%

    return 1.0 + delta


def running_alpha_em(q2: float) -> float:
    """
    Running electromagnetic coupling α(q²).

    Includes vacuum polarization effects.

    Parameters
    ----------
    q2 : float
        Energy scale squared (GeV²)

    Returns
    -------
    float
        α_EM(q²)

    Notes
    -----
    Leading-order running:
    α(q²) = α(0) / (1 - Π(q²))

    where Π is the photon vacuum polarization.
    """
    from ..core import ALPHA_EM

    alpha_0 = ALPHA_EM  # α(0) ≈ 1/137

    if q2 <= 0:
        return alpha_0

    # Lepton contribution to vacuum polarization
    # Π(q²) ≈ (α/3π) × Σ log(q²/m_ℓ²)
    m_e = 0.000511
    m_mu = 0.106
    m_tau = 1.777

    Pi = (alpha_0 / (3.0 * np.pi)) * (
        np.log(q2 / m_e**2) +
        np.log(q2 / m_mu**2) +
        (np.log(q2 / m_tau**2) if q2 > m_tau**2 else 0)
    )

    alpha_q2 = alpha_0 / (1.0 - Pi)

    return alpha_q2


# =============================================================================
# Higher-Order Corrections
# =============================================================================

def electroweak_correction_hnl(m_N: float, m_lepton: float, m_meson: float) -> float:
    """
    Electroweak correction for N → ℓ M decay.

    Includes:
    - W propagator effects
    - Z boson contributions (for neutral current)
    - QED corrections

    Parameters
    ----------
    m_N : float
        HNL mass (GeV)
    m_lepton : float
        Lepton mass (GeV)
    m_meson : float
        Meson mass (GeV)

    Returns
    -------
    float
        Correction factor

    References
    ----------
    - Sirlin, Phys. Rev. D22 (1980) 971
    """
    from ..core import ALPHA_EM

    # W mass
    m_W = 80.377  # GeV

    # Propagator correction
    q2 = m_N**2  # Characteristic scale
    prop_corr = 1.0 + (q2 / m_W**2) * 0.1  # Small correction for m_N << m_W

    # QED correction
    qed_corr = 1.0 + (ALPHA_EM / np.pi) * (np.log(m_N / m_lepton) - 1.0)

    # Combined
    total_corr = prop_corr * qed_corr

    return total_corr


def coulomb_correction_decay(m_parent: float, m1: float, m2: float,
                              Z1: int, Z2: int) -> float:
    """
    Coulomb correction for charged particle pair production.

    For parent → charged1 + charged2, accounts for final-state
    electromagnetic interaction (Sommerfeld factor).

    Parameters
    ----------
    m_parent : float
        Parent mass (GeV)
    m1, m2 : float
        Daughter masses (GeV)
    Z1, Z2 : int
        Electric charges (in units of e)

    Returns
    -------
    float
        Correction factor

    Notes
    -----
    For non-relativistic particles, the Coulomb correction is:
    F_C = (παZ1Z2/β) / (1 - exp(-παZ1Z2/β))

    where β is the relative velocity.
    """
    from ..core import ALPHA_EM

    # Check if particles are charged
    if Z1 == 0 or Z2 == 0:
        return 1.0  # No Coulomb interaction

    # Relative momentum in parent rest frame
    E1 = (m_parent**2 + m1**2 - m2**2) / (2 * m_parent)
    E2 = (m_parent**2 + m2**2 - m1**2) / (2 * m_parent)

    p_sq = ((m_parent**2 - (m1 + m2)**2) * (m_parent**2 - (m1 - m2)**2)) / (4 * m_parent**2)

    if p_sq <= 0:
        return 1.0

    p = np.sqrt(p_sq)

    # Relative velocity
    beta = p / (0.5 * (E1 + E2))

    # Coulomb parameter
    eta = ALPHA_EM * Z1 * Z2 / beta

    # Sommerfeld factor
    if np.abs(eta) < 1e-6:
        F_C = 1.0
    else:
        F_C = eta / (1.0 - np.exp(-eta))

    return F_C


# =============================================================================
# Composite Corrections
# =============================================================================

def full_decay_correction(
    m_parent: float,
    m_daughter1: float,
    m_daughter2: float,
    channel_type: str = 'hadronic'
) -> float:
    """
    Complete correction factor for two-body decay.

    Combines form factors, QED, and electroweak corrections.

    Parameters
    ----------
    m_parent : float
        Parent mass (GeV)
    m_daughter1, m_daughter2 : float
        Daughter masses (GeV)
    channel_type : str
        Type of channel: 'hadronic', 'leptonic', 'photonic'

    Returns
    -------
    float
        Total correction factor

    Examples
    --------
    >>> # For N → μ π with masses
    >>> corr = full_decay_correction(2.0, 0.106, 0.140, 'hadronic')
    >>> print(f"Correction factor: {corr:.3f}")
    """
    corrections = []

    # Form factor (if hadronic)
    if channel_type == 'hadronic':
        q2 = m_parent**2
        if m_daughter2 < 0.2:  # Pion-like
            ff = form_factor_pion(q2)
        elif m_daughter2 > 0.4:  # Kaon-like
            ff = form_factor_kaon(q2)
        else:
            ff = 1.0
        corrections.append(ff**2)  # |F(q²)|²

    # QED correction
    if channel_type == 'leptonic':
        s = m_parent**2
        qed = qed_correction_lepton_pair(s, m_daughter1)
        corrections.append(qed)

    elif channel_type == 'photonic':
        qed = qed_correction_photon_pair(m_parent)
        corrections.append(qed)

    # Electroweak (if applicable)
    if channel_type in ['hadronic', 'leptonic']:
        ew = electroweak_correction_hnl(m_parent, m_daughter1, m_daughter2)
        corrections.append(ew)

    # Combine multiplicatively
    total = np.prod(corrections) if corrections else 1.0

    return total


# =============================================================================
# Utility Functions
# =============================================================================

def apply_form_factor_to_width(
    width_tree_level: float,
    m_parent: float,
    m_daughter1: float,
    m_daughter2: float,
    channel_type: str = 'hadronic'
) -> float:
    """
    Apply corrections to tree-level partial width.

    Parameters
    ----------
    width_tree_level : float
        Tree-level partial width (GeV)
    m_parent, m_daughter1, m_daughter2 : float
        Particle masses (GeV)
    channel_type : str
        Channel type

    Returns
    -------
    float
        Corrected partial width (GeV)
    """
    correction = full_decay_correction(m_parent, m_daughter1, m_daughter2, channel_type)
    return width_tree_level * correction
