"""
Transverse detector geometry and acceptance.

For detectors placed perpendicular to the beam axis (MATHUSLA, ANUBIS).
Different from forward detectors (FASER, SHiP) - acceptance depends on
transverse momentum and rapidity.
"""

import numpy as np
from typing import Tuple, Optional
from dataclasses import dataclass


@dataclass
class TransverseDetectorGeometry:
    """
    Geometry for transverse LLP detector.

    Detector placed perpendicular to beam axis (e.g., above interaction point).

    Parameters
    ----------
    x_center, y_center, z_center : float
        Center position (m) in detector coordinates
        - x: horizontal transverse
        - y: vertical (typically upward from IP)
        - z: beam direction
    dx, dy, dz : float
        Half-widths of fiducial volume (m)
    decay_volume_min : float
        Minimum distance from IP to start of decay volume (m)
    decay_volume_max : float
        Maximum distance from IP for decay volume (m)
    """
    x_center: float
    y_center: float
    z_center: float
    dx: float
    dy: float
    dz: float
    decay_volume_min: float = 5.0
    decay_volume_max: float = 200.0

    def fiducial_bounds(self) -> dict:
        """Get fiducial volume boundaries."""
        return {
            'x_min': self.x_center - self.dx,
            'x_max': self.x_center + self.dx,
            'y_min': self.y_center - self.dy,
            'y_max': self.y_center + self.dy,
            'z_min': self.z_center - self.dz,
            'z_max': self.z_center + self.dz,
        }

    def is_in_fiducial(self, x: float, y: float, z: float) -> bool:
        """Check if point is in fiducial volume."""
        bounds = self.fiducial_bounds()
        return (bounds['x_min'] <= x <= bounds['x_max'] and
                bounds['y_min'] <= y <= bounds['y_max'] and
                bounds['z_min'] <= z <= bounds['z_max'])


# MATHUSLA geometry (arXiv:1901.04346)
MATHUSLA = TransverseDetectorGeometry(
    x_center=0.0,
    y_center=100.0,  # 100m above IP
    z_center=0.0,
    dx=100.0,  # 200m × 200m footprint
    dy=10.0,   # 20m height
    dz=100.0,
    decay_volume_min=5.0,
    decay_volume_max=100.0  # Up to detector entrance
)

# ANUBIS geometry (similar to MATHUSLA)
ANUBIS = TransverseDetectorGeometry(
    x_center=0.0,
    y_center=30.0,  # ~30m above IP
    z_center=0.0,
    dx=25.0,   # 50m × 50m
    dy=5.0,    # 10m height
    dz=25.0,
    decay_volume_min=5.0,
    decay_volume_max=50.0
)


def decay_position(production_vertex: np.ndarray,
                   particle_4vec: np.ndarray,
                   decay_length: float) -> np.ndarray:
    """
    Calculate decay position given production vertex and decay length.

    Parameters
    ----------
    production_vertex : np.ndarray, shape (3,)
        Production vertex (x, y, z) in meters
    particle_4vec : np.ndarray, shape (4,)
        Particle 4-vector [E, px, py, pz] in GeV
    decay_length : float
        Proper decay length c·τ in meters

    Returns
    -------
    np.ndarray, shape (3,)
        Decay position (x, y, z) in meters
    """
    E = particle_4vec[0]
    p_vec = particle_4vec[1:4]
    p_mag = np.linalg.norm(p_vec)

    # Boost factor γβ
    gamma = E / np.sqrt(E**2 - p_mag**2)
    beta = p_mag / E

    # Decay length in lab frame
    L_lab = gamma * beta * decay_length

    # Direction unit vector
    direction = p_vec / p_mag

    # Decay position
    decay_pos = production_vertex + L_lab * direction

    return decay_pos


def geometric_acceptance(particle_4vec: np.ndarray,
                        geometry: TransverseDetectorGeometry,
                        production_vertex: np.ndarray = None) -> float:
    """
    Calculate geometric acceptance for LLP to reach detector and decay inside.

    For transverse detectors, particle must:
    1. Have trajectory intersecting detector volume
    2. Decay probability in decay volume leading to detector

    Parameters
    ----------
    particle_4vec : np.ndarray, shape (4,) or (n, 4)
        LLP 4-vector(s) [E, px, py, pz] in GeV
    geometry : TransverseDetectorGeometry
        Detector geometry
    production_vertex : np.ndarray, shape (3,), optional
        Production point (x, y, z) in meters. Default is (0, 0, 0).

    Returns
    -------
    float or np.ndarray
        Acceptance fraction [0, 1]
    """
    if production_vertex is None:
        production_vertex = np.array([0.0, 0.0, 0.0])

    particle_4vec = np.atleast_2d(particle_4vec)
    n_events = len(particle_4vec)

    acceptance = np.zeros(n_events)

    for i in range(n_events):
        p4 = particle_4vec[i]
        E = p4[0]
        p_vec = p4[1:4]
        p_mag = np.linalg.norm(p_vec)

        # Direction unit vector
        direction = p_vec / p_mag

        # Find intersection with detector volume
        # Parametric line: r(t) = r0 + t * direction
        # Find t where line enters/exits detector

        bounds = geometry.fiducial_bounds()

        # Find t for each boundary plane
        t_x_min = (bounds['x_min'] - production_vertex[0]) / direction[0] if direction[0] != 0 else np.inf
        t_x_max = (bounds['x_max'] - production_vertex[0]) / direction[0] if direction[0] != 0 else np.inf
        t_y_min = (bounds['y_min'] - production_vertex[1]) / direction[1] if direction[1] != 0 else np.inf
        t_y_max = (bounds['y_max'] - production_vertex[1]) / direction[1] if direction[1] != 0 else np.inf
        t_z_min = (bounds['z_min'] - production_vertex[2]) / direction[2] if direction[2] != 0 else np.inf
        t_z_max = (bounds['z_max'] - production_vertex[2]) / direction[2] if direction[2] != 0 else np.inf

        # Entry and exit t values
        t_entry = max(min(t_x_min, t_x_max), min(t_y_min, t_y_max), min(t_z_min, t_z_max), 0)
        t_exit = min(max(t_x_min, t_x_max), max(t_y_min, t_y_max), max(t_z_min, t_z_max))

        if t_entry >= t_exit or t_exit < 0:
            # No intersection
            acceptance[i] = 0.0
            continue

        # Particle intersects detector
        # For now, return 1.0 if trajectory goes through
        # In reality, need to check decay probability in volume
        acceptance[i] = 1.0

    return acceptance.squeeze() if n_events == 1 else acceptance


def decay_probability_transverse(particle_4vec: np.ndarray,
                                 ctau: float,
                                 geometry: TransverseDetectorGeometry,
                                 production_vertex: np.ndarray = None) -> float:
    """
    Probability that LLP decays in the decay volume leading to detector.

    P = P(survive to d_min) × P(decay before d_max)
      = exp(-d_min/λ) - exp(-d_max/λ)

    where λ = γβ·cτ is the boosted decay length.

    Parameters
    ----------
    particle_4vec : np.ndarray, shape (4,) or (n, 4)
        LLP 4-vector(s) [E, px, py, pz]
    ctau : float
        Proper decay length c·τ in meters
    geometry : TransverseDetectorGeometry
        Detector geometry
    production_vertex : np.ndarray, shape (3,), optional
        Production vertex in meters

    Returns
    -------
    float or np.ndarray
        Decay probability in relevant volume
    """
    if production_vertex is None:
        production_vertex = np.array([0.0, 0.0, 0.0])

    particle_4vec = np.atleast_2d(particle_4vec)
    n_events = len(particle_4vec)

    prob = np.zeros(n_events)

    for i in range(n_events):
        p4 = particle_4vec[i]
        E = p4[0]
        p_vec = p4[1:4]
        p_mag = np.linalg.norm(p_vec)

        # Boosted decay length
        m = np.sqrt(E**2 - p_mag**2)
        gamma = E / m
        beta = p_mag / E
        lambda_lab = gamma * beta * ctau

        # Distance to detector (minimum decay distance)
        d_min = geometry.decay_volume_min
        d_max = geometry.decay_volume_max

        # Decay probability in [d_min, d_max]
        prob[i] = np.exp(-d_min / lambda_lab) - np.exp(-d_max / lambda_lab)

    return prob.squeeze() if n_events == 1 else prob


def signal_yield_transverse(particle_4vecs: np.ndarray,
                           ctau: float,
                           geometry: TransverseDetectorGeometry,
                           weights: np.ndarray = None,
                           min_charged: int = 2) -> float:
    """
    Calculate signal yield for transverse detector.

    N_signal = Σ w_i × P_decay × P_acceptance × P_reconstruction

    Parameters
    ----------
    particle_4vecs : np.ndarray, shape (n, 4)
        LLP 4-vectors from production
    ctau : float
        Proper decay length (m)
    geometry : TransverseDetectorGeometry
        Detector geometry
    weights : np.ndarray, shape (n,), optional
        Event weights (e.g., cross-section × BR)
    min_charged : int
        Minimum charged tracks for reconstruction

    Returns
    -------
    float
        Expected signal yield
    """
    n_events = len(particle_4vecs)

    if weights is None:
        weights = np.ones(n_events)

    # Geometric acceptance
    acc = geometric_acceptance(particle_4vecs, geometry)

    # Decay probability
    p_decay = decay_probability_transverse(particle_4vecs, ctau, geometry)

    # Reconstruction efficiency (simplified: assume 1.0 if visible)
    # In reality, depends on decay channel and kinematics
    eff_reco = np.ones(n_events)  # Placeholder

    # Total signal
    signal = np.sum(weights * acc * p_decay * eff_reco)

    return signal
