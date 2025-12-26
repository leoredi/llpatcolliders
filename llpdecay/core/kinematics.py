"""
Relativistic kinematics for LLP decays.

Handles Lorentz transformations, phase space calculations, and momentum sampling.
"""

import numpy as np
from typing import Tuple, Optional


def kallen(a: float, b: float, c: float) -> float:
    """
    Källén triangle function λ(a,b,c).

    λ(a,b,c) = a² + b² + c² - 2ab - 2bc - 2ca

    Used in two-body decay kinematics and other phase space calculations.

    Parameters
    ----------
    a, b, c : float
        Input values (typically squared masses)

    Returns
    -------
    float
        Value of Källén function
    """
    return a**2 + b**2 + c**2 - 2*a*b - 2*b*c - 2*c*a


def four_vector_mass(p: np.ndarray) -> np.ndarray:
    """
    Calculate invariant mass from 4-vector(s).

    m² = E² - px² - py² - pz²

    Parameters
    ----------
    p : np.ndarray, shape (..., 4)
        4-vector(s) [E, px, py, pz]

    Returns
    -------
    np.ndarray
        Invariant mass(es)
    """
    p = np.atleast_2d(p)
    m2 = p[..., 0]**2 - p[..., 1]**2 - p[..., 2]**2 - p[..., 3]**2
    return np.sqrt(np.maximum(m2, 0))  # Protect against numerical errors


def boost_to_lab(p_rest: np.ndarray, parent_4vec: np.ndarray) -> np.ndarray:
    """
    Boost 4-vector(s) from parent rest frame to lab frame.

    Uses the standard Lorentz boost transformation. Handles both single
    events and arrays of events efficiently.

    Parameters
    ----------
    p_rest : np.ndarray, shape (..., 4)
        4-vectors in parent rest frame [E, px, py, pz]
    parent_4vec : np.ndarray, shape (4,)
        Parent 4-vector in lab frame [E, px, py, pz]

    Returns
    -------
    np.ndarray, shape (..., 4)
        4-vectors in lab frame

    Examples
    --------
    >>> parent = np.array([10.0, 3.0, 0.0, 9.5])
    >>> daughters_rest = np.array([[1.0, 0.5, 0.0, 0.866],
    ...                            [1.0, -0.5, 0.0, -0.866]])
    >>> daughters_lab = boost_to_lab(daughters_rest, parent)
    """
    p_rest = np.atleast_2d(p_rest)
    original_shape = p_rest.shape

    E_parent = parent_4vec[0]
    p_parent = parent_4vec[1:4]
    m_parent = four_vector_mass(parent_4vec)

    # Boost parameters
    gamma = E_parent / m_parent
    beta_vec = p_parent / E_parent
    beta = np.linalg.norm(beta_vec)

    # No boost needed if parent at rest
    if beta < 1e-10:
        return p_rest

    n_hat = beta_vec / beta  # Boost direction

    # Extract energy and 3-momentum
    E_rest = p_rest[..., 0]
    p_rest_3 = p_rest[..., 1:4]

    # Decompose momentum into parallel and perpendicular components
    # p_parallel = (p · n̂) n̂
    p_dot_n = np.sum(p_rest_3 * n_hat, axis=-1, keepdims=True)
    p_parallel = p_dot_n * n_hat
    p_perp = p_rest_3 - p_parallel

    # Apply Lorentz boost
    E_lab = gamma * (E_rest + beta * p_dot_n.squeeze(-1))
    p_parallel_lab = gamma * (p_parallel + beta * E_rest[..., np.newaxis] * n_hat)
    p_lab_3 = p_parallel_lab + p_perp

    # Reconstruct 4-vector
    result = np.concatenate([E_lab[..., np.newaxis], p_lab_3], axis=-1)
    return result.reshape(original_shape)


def two_body_decay_momenta(m_parent: float, m1: float, m2: float) -> Tuple[float, float, float]:
    """
    Calculate daughter momenta magnitude in parent rest frame for two-body decay.

    In the rest frame of the parent, the two daughters are emitted back-to-back
    with equal momentum magnitude p. This function calculates p and the energies.

    Parameters
    ----------
    m_parent : float
        Parent particle mass
    m1, m2 : float
        Daughter particle masses

    Returns
    -------
    p_mag : float
        Momentum magnitude of both daughters (back-to-back)
    E1, E2 : float
        Energies of daughters

    Raises
    ------
    ValueError
        If decay is kinematically forbidden (m_parent < m1 + m2)

    Notes
    -----
    Energy-momentum conservation gives:
        E1 + E2 = m_parent
        |p1| = |p2| = p

    Solving these with E² = p² + m²:
        p = √λ(m_parent², m1², m2²) / (2 m_parent)
        E1 = (m_parent² + m1² - m2²) / (2 m_parent)
        E2 = (m_parent² + m2² - m1²) / (2 m_parent)
    """
    lam = kallen(m_parent**2, m1**2, m2**2)

    if lam < 0:
        raise ValueError(
            f"Kinematically forbidden decay: m_parent={m_parent:.4f} GeV, "
            f"m1={m1:.4f} GeV, m2={m2:.4f} GeV (sum={m1+m2:.4f} GeV)"
        )

    p_mag = np.sqrt(lam) / (2 * m_parent)
    E1 = (m_parent**2 + m1**2 - m2**2) / (2 * m_parent)
    E2 = (m_parent**2 + m2**2 - m1**2) / (2 * m_parent)

    return p_mag, E1, E2


def sample_two_body_decay(
    m_parent: float,
    m1: float,
    m2: float,
    n_events: int = 1,
    polarization: float = 0.0,
    alpha: float = 1.0,
    rng: Optional[np.random.Generator] = None
) -> np.ndarray:
    """
    Sample two-body decay in parent rest frame.

    Generates 4-vectors for decay products with proper kinematics and optional
    angular asymmetry from parent polarization.

    Parameters
    ----------
    m_parent, m1, m2 : float
        Masses of parent and daughters (GeV)
    n_events : int
        Number of decays to sample
    polarization : float, default=0.0
        Parent polarization along z-axis, P ∈ [-1, 1]
        - P = 0: isotropic decay
        - P > 0: prefer forward emission
        - P < 0: prefer backward emission
    alpha : float, default=1.0
        Asymmetry parameter for the decay channel
        Angular distribution: dΓ/dcosθ ∝ (1 + P·α·cosθ)
        - α = 1: maximal correlation (e.g., π → μ ν)
        - α = 0: no correlation
        - α = -1: anti-correlation
    rng : np.random.Generator, optional
        Random number generator. If None, creates new default generator.

    Returns
    -------
    np.ndarray, shape (n_events, 2, 4)
        Daughter 4-vectors [E, px, py, pz] in parent rest frame
        - result[:, 0, :] = first daughter
        - result[:, 1, :] = second daughter (back-to-back)

    Examples
    --------
    >>> # Isotropic pion decay
    >>> daughters = sample_two_body_decay(0.13957, 0.105658, 0.0, n_events=100)
    >>> assert daughters.shape == (100, 2, 4)

    >>> # Polarized decay
    >>> daughters = sample_two_body_decay(
    ...     0.13957, 0.105658, 0.0,
    ...     n_events=1000,
    ...     polarization=1.0,  # Fully polarized
    ...     alpha=1.0  # Maximal correlation
    ... )
    """
    if rng is None:
        rng = np.random.default_rng()

    p_mag, E1, E2 = two_body_decay_momenta(m_parent, m1, m2)

    # Sample azimuthal angle (always uniform)
    phi = rng.uniform(0, 2*np.pi, n_events)

    # Sample polar angle with optional polarization effect
    if abs(polarization * alpha) < 1e-10:
        # Isotropic: uniform in cos(θ)
        cos_theta = rng.uniform(-1, 1, n_events)
    else:
        # Polarized: dΓ/dcosθ ∝ (1 + P·α·cosθ)
        # Use rejection sampling for simplicity
        P_alpha = polarization * alpha
        max_weight = 1 + abs(P_alpha)

        cos_theta = np.zeros(n_events)
        for i in range(n_events):
            while True:
                x = rng.uniform(-1, 1)
                weight = 1 + P_alpha * x
                if rng.uniform(0, max_weight) < weight:
                    cos_theta[i] = x
                    break

    sin_theta = np.sqrt(1 - cos_theta**2)

    # Construct 3-momenta
    px = p_mag * sin_theta * np.cos(phi)
    py = p_mag * sin_theta * np.sin(phi)
    pz = p_mag * cos_theta

    # Build 4-vectors for both daughters (back-to-back)
    p1 = np.stack([np.full(n_events, E1), px, py, pz], axis=-1)
    p2 = np.stack([np.full(n_events, E2), -px, -py, -pz], axis=-1)

    return np.stack([p1, p2], axis=1)


def invariant_mass_from_daughters(daughters: np.ndarray) -> np.ndarray:
    """
    Calculate invariant mass of multi-particle system.

    M² = (ΣE)² - (Σp)²

    Parameters
    ----------
    daughters : np.ndarray, shape (..., n_daughters, 4)
        4-vectors of daughter particles

    Returns
    -------
    np.ndarray
        Invariant mass of the system
    """
    # Sum 4-vectors
    total_4vec = np.sum(daughters, axis=-2)
    return four_vector_mass(total_4vec)


def transverse_momentum(p: np.ndarray) -> np.ndarray:
    """
    Calculate transverse momentum pT = √(px² + py²).

    Parameters
    ----------
    p : np.ndarray, shape (..., 4)
        4-vector(s) [E, px, py, pz]

    Returns
    -------
    np.ndarray
        Transverse momentum
    """
    p = np.atleast_2d(p)
    return np.sqrt(p[..., 1]**2 + p[..., 2]**2)


def rapidity(p: np.ndarray) -> np.ndarray:
    """
    Calculate rapidity y = 0.5 * ln((E + pz) / (E - pz)).

    Parameters
    ----------
    p : np.ndarray, shape (..., 4)
        4-vector(s) [E, px, py, pz]

    Returns
    -------
    np.ndarray
        Rapidity
    """
    p = np.atleast_2d(p)
    E = p[..., 0]
    pz = p[..., 3]
    return 0.5 * np.log((E + pz) / (E - pz))


def pseudorapidity(p: np.ndarray) -> np.ndarray:
    """
    Calculate pseudorapidity η = -ln(tan(θ/2)).

    Parameters
    ----------
    p : np.ndarray, shape (..., 4)
        4-vector(s) [E, px, py, pz]

    Returns
    -------
    np.ndarray
        Pseudorapidity
    """
    p = np.atleast_2d(p)
    pt = transverse_momentum(p)
    pz = p[..., 3]
    theta = np.arctan2(pt, pz)
    return -np.log(np.tan(theta / 2))
