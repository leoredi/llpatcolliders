"""
Decay channel definitions for LLPs.

Defines the final state particles, quantum numbers, and selection rules
for various decay modes.
"""

from dataclasses import dataclass
from typing import List, Tuple
import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from llpdecay.core import MASS, PDG, get_mass, get_pdg, is_charged


@dataclass
class DecayChannel:
    """
    Definition of a decay channel.

    Attributes
    ----------
    name : str
        Channel identifier (e.g., 'mu_pi')
    daughters : List[Tuple[str, int]]
        List of (particle_name, pdg_id) for each daughter
        Negative PDG IDs indicate antiparticles
    n_body : int
        Number of particles in final state
    charged_count : int
        Number of electrically charged particles
    asymmetry_param : float
        Angular asymmetry parameter α for polarized decays
        dΓ/dcosθ ∝ (1 + P·α·cosθ) where P is polarization
    """
    name: str
    daughters: List[Tuple[str, int]]
    n_body: int
    charged_count: int
    asymmetry_param: float = 0.0

    def daughter_masses(self) -> List[float]:
        """Get masses of daughter particles."""
        return [get_mass(d[0]) for d in self.daughters]

    def daughter_pdgs(self) -> List[int]:
        """Get PDG IDs of daughter particles."""
        return [d[1] for d in self.daughters]

    def daughter_names(self) -> List[str]:
        """Get names of daughter particles."""
        return [d[0] for d in self.daughters]

    def threshold_mass(self) -> float:
        """Minimum parent mass for kinematic accessibility."""
        return sum(self.daughter_masses())

    def is_visible(self, pt_threshold: float = 0.5) -> bool:
        """
        Check if channel produces detectable signature.

        Parameters
        ----------
        pt_threshold : float
            Minimum pT requirement (GeV)

        Returns
        -------
        bool
            True if channel has charged tracks that can be detected
        """
        return self.charged_count >= 2  # Need at least 2 charged tracks

    def __repr__(self) -> str:
        particles = ' '.join([d[0] for d in self.daughters])
        return f"DecayChannel({self.name}: {particles}, {self.n_body}-body, {self.charged_count} charged)"


# ============================================================================
# HNL Decay Channels
# ============================================================================

HNL_CHANNELS = {
    # ------------------------------------------------------------------------
    # Charged current: N → ℓ⁻ π⁺ (and charge conjugate)
    # Dominant modes for GeV-scale HNLs
    # ------------------------------------------------------------------------
    'e_pi': DecayChannel(
        name='e_pi',
        daughters=[('electron', -PDG['electron']), ('pi_charged', PDG['pi_charged'])],
        n_body=2,
        charged_count=2,
        asymmetry_param=1.0,  # V-A coupling
    ),

    'mu_pi': DecayChannel(
        name='mu_pi',
        daughters=[('muon', -PDG['muon']), ('pi_charged', PDG['pi_charged'])],
        n_body=2,
        charged_count=2,
        asymmetry_param=1.0,
    ),

    'tau_pi': DecayChannel(
        name='tau_pi',
        daughters=[('tau', -PDG['tau']), ('pi_charged', PDG['pi_charged'])],
        n_body=2,
        charged_count=2,
        asymmetry_param=1.0,
    ),

    # ------------------------------------------------------------------------
    # Charged current with kaon: N → ℓ⁻ K⁺
    # Cabibbo-suppressed relative to pion modes
    # ------------------------------------------------------------------------
    'e_K': DecayChannel(
        name='e_K',
        daughters=[('electron', -PDG['electron']), ('K_charged', PDG['K_charged'])],
        n_body=2,
        charged_count=2,
        asymmetry_param=1.0,
    ),

    'mu_K': DecayChannel(
        name='mu_K',
        daughters=[('muon', -PDG['muon']), ('K_charged', PDG['K_charged'])],
        n_body=2,
        charged_count=2,
        asymmetry_param=1.0,
    ),

    'tau_K': DecayChannel(
        name='tau_K',
        daughters=[('tau', -PDG['tau']), ('K_charged', PDG['K_charged'])],
        n_body=2,
        charged_count=2,
        asymmetry_param=1.0,
    ),

    # ------------------------------------------------------------------------
    # Charged current with rho: N → ℓ⁻ ρ⁺
    # Important at higher masses (> 1 GeV)
    # ------------------------------------------------------------------------
    'e_rho': DecayChannel(
        name='e_rho',
        daughters=[('electron', -PDG['electron']), ('rho_charged', PDG['rho_charged'])],
        n_body=2,
        charged_count=2,
        asymmetry_param=1.0,
    ),

    'mu_rho': DecayChannel(
        name='mu_rho',
        daughters=[('muon', -PDG['muon']), ('rho_charged', PDG['rho_charged'])],
        n_body=2,
        charged_count=2,
        asymmetry_param=1.0,
    ),

    'tau_rho': DecayChannel(
        name='tau_rho',
        daughters=[('tau', -PDG['tau']), ('rho_charged', PDG['rho_charged'])],
        n_body=2,
        charged_count=2,
        asymmetry_param=1.0,
    ),

    # ------------------------------------------------------------------------
    # Neutral current: N → ν π⁰
    # INVISIBLE - neutrino escapes, π⁰ decays to photons
    # ------------------------------------------------------------------------
    'nu_e_pi0': DecayChannel(
        name='nu_e_pi0',
        daughters=[('nu_e', PDG['nu_e']), ('pi_neutral', PDG['pi_neutral'])],
        n_body=2,
        charged_count=0,  # Invisible!
        asymmetry_param=0.0,  # Neutrino undetected
    ),

    'nu_mu_pi0': DecayChannel(
        name='nu_mu_pi0',
        daughters=[('nu_mu', PDG['nu_mu']), ('pi_neutral', PDG['pi_neutral'])],
        n_body=2,
        charged_count=0,
        asymmetry_param=0.0,
    ),

    'nu_tau_pi0': DecayChannel(
        name='nu_tau_pi0',
        daughters=[('nu_tau', PDG['nu_tau']), ('pi_neutral', PDG['pi_neutral'])],
        n_body=2,
        charged_count=0,
        asymmetry_param=0.0,
    ),

    # ------------------------------------------------------------------------
    # 3-body leptonic: N → ν ℓ⁺ ℓ⁻
    # Proceeds via virtual W* or Z*
    # Phase 2 implementation (requires 3-body phase space)
    # ------------------------------------------------------------------------
    'nu_e_e': DecayChannel(
        name='nu_e_e',
        daughters=[
            ('nu_e', PDG['nu_e']),
            ('electron', PDG['electron']),
            ('electron', -PDG['electron'])
        ],
        n_body=3,
        charged_count=2,
        asymmetry_param=0.0,
    ),

    'nu_mu_mu': DecayChannel(
        name='nu_mu_mu',
        daughters=[
            ('nu_mu', PDG['nu_mu']),
            ('muon', PDG['muon']),
            ('muon', -PDG['muon'])
        ],
        n_body=3,
        charged_count=2,
        asymmetry_param=0.0,
    ),

    'nu_tau_tau': DecayChannel(
        name='nu_tau_tau',
        daughters=[
            ('nu_tau', PDG['nu_tau']),
            ('tau', PDG['tau']),
            ('tau', -PDG['tau'])
        ],
        n_body=3,
        charged_count=2,
        asymmetry_param=0.0,
    ),

    # Mixed flavor 3-body
    'nu_e_mu': DecayChannel(
        name='nu_e_mu',
        daughters=[
            ('nu_e', PDG['nu_e']),
            ('electron', PDG['electron']),
            ('muon', -PDG['muon'])
        ],
        n_body=3,
        charged_count=2,
        asymmetry_param=0.0,
    ),

    'nu_mu_e': DecayChannel(
        name='nu_mu_e',
        daughters=[
            ('nu_mu', PDG['nu_mu']),
            ('muon', PDG['muon']),
            ('electron', -PDG['electron'])
        ],
        n_body=3,
        charged_count=2,
        asymmetry_param=0.0,
    ),
}


# ============================================================================
# Helper functions
# ============================================================================

def get_channel(channel_name: str) -> DecayChannel:
    """
    Get decay channel by name.

    Parameters
    ----------
    channel_name : str
        Name of the channel

    Returns
    -------
    DecayChannel
        Channel object

    Raises
    ------
    KeyError
        If channel not found
    """
    if channel_name not in HNL_CHANNELS:
        raise KeyError(f"Unknown channel: {channel_name}. Available: {list(HNL_CHANNELS.keys())}")
    return HNL_CHANNELS[channel_name]


def list_channels(min_mass: float = 0.0, max_mass: float = float('inf'),
                  visible_only: bool = False) -> List[str]:
    """
    List available decay channels.

    Parameters
    ----------
    min_mass : float
        Minimum parent mass (GeV)
    max_mass : float
        Maximum parent mass (GeV)
    visible_only : bool
        If True, only return channels with charged tracks

    Returns
    -------
    list of str
        Channel names
    """
    channels = []
    for name, ch in HNL_CHANNELS.items():
        threshold = ch.threshold_mass()
        if min_mass <= threshold <= max_mass:
            if not visible_only or ch.is_visible():
                channels.append(name)
    return channels


def group_by_flavor(channels: List[str]) -> dict:
    """
    Group channels by lepton flavor.

    Parameters
    ----------
    channels : list of str
        Channel names

    Returns
    -------
    dict
        {flavor: [channel_names]} where flavor is 'e', 'mu', or 'tau'
    """
    groups = {'e': [], 'mu': [], 'tau': [], 'mixed': []}

    for ch_name in channels:
        if ch_name.startswith('e_'):
            groups['e'].append(ch_name)
        elif ch_name.startswith('mu_'):
            groups['mu'].append(ch_name)
        elif ch_name.startswith('tau_'):
            groups['tau'].append(ch_name)
        elif ch_name.startswith('nu_'):
            # Neutral current - assign to flavor based on neutrino
            if 'nu_e' in ch_name:
                groups['e'].append(ch_name)
            elif 'nu_mu' in ch_name:
                groups['mu'].append(ch_name)
            elif 'nu_tau' in ch_name:
                groups['tau'].append(ch_name)
        else:
            groups['mixed'].append(ch_name)

    return groups
