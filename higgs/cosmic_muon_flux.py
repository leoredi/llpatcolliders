"""
Cosmic-muon depth (post-overburden) momentum/energy spectrum for the PX56 tunnel.

Surface flux: modified Gaisser parameterisation with the Chirkin cos(theta*)
Earth-curvature correction. Energy loss through the rock overburden is the
continuous-slowing-down approximation dE/dx = -(a + b E), so a muon arriving at
depth with E_depth started at the surface with
    E_surface = (E_depth + a/b) exp(b X) - a/b,   X = slant depth [g/cm^2].

This module is the single source for the depth spectrum, used by both
`cosmic_decay_check.py` (`--spectrum depth`) and `background_trident_update.py`.
The default cosmic-decay spectrum is the measured one read from a ROOT file
(`--spectrum empirical`); this CSDA model is the analytic fallback.
"""
import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d

from grendel_geometry import Y_POSITION

# Rock overburden properties
ROCK_RHO      = 2.65     # g/cm^3 (standard rock)
ROCK_A_LOSS   = 2.0e-3   # GeV / (g/cm^2), ionization
ROCK_B_LOSS   = 4.4e-6   # 1 / (g/cm^2), radiative
SURFACE_TO_IP = 100.0    # m, surface elevation above the CMS IP
OVERBURDEN    = SURFACE_TO_IP - Y_POSITION   # ~78 m vertical


class CosmicMuonFlux:
    """
    Cosmic muon flux at the tunnel depth.

    Surface flux: modified Gaisser parameterisation with Chirkin cos(theta*)
    correction for Earth curvature. Energy loss in rock: dE/dx = -(a + bE),
    E_surface = (E_depth + a/b) exp(b X) - a/b with X the slant depth [g/cm^2].
    """

    def __init__(self, overburden_m=OVERBURDEN, rock_rho=ROCK_RHO):
        self.overburden_m = overburden_m
        self.rock_rho = rock_rho
        self.a = ROCK_A_LOSS   # GeV / (g/cm^2)
        self.b = ROCK_B_LOSS   # 1 / (g/cm^2)

    def slant_depth(self, cos_theta):
        """Slant depth in g/cm^2 for a muon at zenith angle theta."""
        return self.rock_rho * self.overburden_m / max(cos_theta, 0.05) * 100  # m->cm

    def min_surface_energy(self, cos_theta):
        """Minimum surface energy (GeV) to reach the tunnel at angle theta."""
        X = self.slant_depth(cos_theta)
        return (self.a / self.b) * (np.exp(self.b * X) - 1)

    def energy_at_depth(self, E_surface, cos_theta):
        """Muon energy at tunnel depth given surface energy and angle."""
        X = self.slant_depth(cos_theta)
        E_depth = (E_surface + self.a / self.b) * np.exp(-self.b * X) - self.a / self.b
        return max(E_depth, 0.0)

    def surface_energy(self, E_depth, cos_theta):
        """Surface energy required to arrive at the tunnel with E_depth."""
        X = self.slant_depth(cos_theta)
        return (E_depth + self.a / self.b) * np.exp(self.b * X) - self.a / self.b

    @staticmethod
    def cos_theta_star(cos_theta):
        """Chirkin cos(theta*) correction for Earth curvature."""
        p = [0.102573, -0.068287, 0.958633, 0.0407253, 0.817285]
        cs2 = (cos_theta**2 + p[0]**2 + p[1] * cos_theta**p[2]
               + p[3] * cos_theta**p[4])
        return np.sqrt(cs2 / (1 + p[0]**2 + p[1] + p[3]))

    def surface_flux(self, E, cos_theta):
        """Gaisser surface flux dI/dE (cm^-2 s^-1 sr^-1 GeV^-1), modified
        Gaisser formula with the cos(theta*) correction."""
        if E < 1.0:
            return 0.0
        cs = self.cos_theta_star(cos_theta)
        return 0.14 * E**(-2.7) * (
            1.0 / (1.0 + 1.1 * E * cs / 115.0) +
            0.054 / (1.0 + 1.1 * E * cs / 850.0))

    def flux_at_depth(self, E_depth, cos_theta):
        """Differential flux at tunnel depth (cm^-2 s^-1 sr^-1 GeV^-1):
        dI/dE_depth = dI/dE_surface * dE_surface/dE_depth."""
        E_surf = self.surface_energy(E_depth, cos_theta)
        if E_surf <= 0 or E_depth <= 0:
            return 0.0
        X = self.slant_depth(cos_theta)
        dEdE = np.exp(self.b * X)  # dE_surface / dE_depth
        return self.surface_flux(E_surf, cos_theta) * dEdE

    def integrated_flux_at_depth(self, cos_theta, E_min=1.0, E_max=2000.0):
        """Total muon flux at depth for a given zenith angle (cm^-2 s^-1 sr^-1)."""
        E_min_eff = max(E_min, 0.5)  # muons below ~0.5 GeV don't reach
        result, _ = quad(lambda E: self.flux_at_depth(E, cos_theta),
                         E_min_eff, E_max, limit=200)
        return result

    def sample_energy_at_depth(self, cos_theta, n, rng, E_max=2000.0):
        """Inverse-CDF sampling of muon energies at tunnel depth."""
        E_min = max(self.min_surface_energy(cos_theta) * 0.01, 0.5)
        E_grid = np.geomspace(E_min, E_max, 500)
        pdf = np.array([self.flux_at_depth(E, cos_theta) for E in E_grid])
        pdf = np.maximum(pdf, 0)
        if pdf.sum() == 0:
            return np.full(n, E_min)
        cdf = np.cumsum(pdf)
        cdf /= cdf[-1]
        inv = interp1d(cdf, E_grid, bounds_error=False,
                       fill_value=(E_min, E_max))
        return inv(rng.uniform(0, 1, n))
