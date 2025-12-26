"""
Axion-Like Particle (ALP) model implementation.

Implements pseudoscalar ALP with configurable couplings to photons,
leptons, and quarks.

References:
- Bauer et al., arXiv:1708.00443 (JHEP 2017)
- Calibbi et al., arXiv:2006.04795 (JHEP 2020)
- Marciano et al., Phys.Rev.D94:115033, arXiv:1607.01022
"""

import numpy as np
from typing import Dict, Optional, Tuple, List, Union
import warnings

from .base import LLPModel
from ..core import (
    ALPHA_EM, HBAR_C_GEV_M,
    sample_two_body_decay,
    boost_to_lab,
    MASS, PDG
)


class ALP(LLPModel):
    """
    Axion-Like Particle model.

    A pseudoscalar particle coupling to SM via effective operators.
    Main couplings:
    - Photons: a γγ → g_aγγ a F F̃
    - Leptons: a ℓ ℓ → (c_ℓ/f_a) a ℓ̄ γ⁵ ℓ
    - Quarks: a qq → (c_q/f_a) a q̄ γ⁵ q

    Parameters
    ----------
    mass : float
        ALP mass in GeV
    g_agg : float, default=0
        Coupling to photons (GeV⁻¹)
    f_a : float, default=1e9
        ALP decay constant (GeV)
    c_e, c_mu, c_tau : float, default=0
        Dimensionless lepton Yukawa-like couplings
    c_u, c_d, c_s : float, default=0
        Dimensionless quark couplings
    seed : int, optional
        Random seed

    Attributes
    ----------
    mass : float
        ALP mass in GeV
    g_agg : float
        Photon coupling
    f_a : float
        Decay constant
    c_leptons : dict
        Lepton couplings {flavor: coupling}
    c_quarks : dict
        Quark couplings {flavor: coupling}

    Examples
    --------
    >>> # Photophilic ALP
    >>> alp = ALP(mass=0.5, g_agg=1e-5)  # Couples mainly to photons
    >>> brs = alp.branching_ratios()
    >>> print(f"BR(a → γγ) = {brs.get('gamma_gamma', 0):.3f}")

    >>> # Leptophilic ALP
    >>> alp = ALP(mass=1.0, f_a=1e8, c_e=1.0, c_mu=1.0)
    >>> daughters = alp.sample_decay(parent_4vec)
    """

    def __init__(
        self,
        mass: float,
        g_agg: float = 0.0,
        f_a: float = 1e9,
        c_e: float = 0.0,
        c_mu: float = 0.0,
        c_tau: float = 0.0,
        c_u: float = 0.0,
        c_d: float = 0.0,
        c_s: float = 0.0,
        seed: Optional[int] = None
    ):
        super().__init__(mass)

        self.g_agg = g_agg
        self.f_a = f_a

        self.c_leptons = {
            'e': c_e,
            'mu': c_mu,
            'tau': c_tau,
        }

        self.c_quarks = {
            'u': c_u,
            'd': c_d,
            's': c_s,
        }

        self._br_cache: Optional[Dict[str, float]] = None
        self._total_width_cache: Optional[float] = None
        self._rng = np.random.default_rng(seed)

        # Check if any couplings are non-zero
        if self.g_agg == 0 and all(c == 0 for c in self.c_leptons.values()) \
           and all(c == 0 for c in self.c_quarks.values()):
            warnings.warn("All couplings are zero - ALP is stable!")

    def available_channels(self) -> List[str]:
        """
        List decay channels kinematically accessible.

        Returns
        -------
        list of str
            Channel names
        """
        channels = []

        # Photon decays
        if self.g_agg != 0:
            # a → γγ (always kinematically allowed)
            channels.append('gamma_gamma')

        # Leptonic decays
        for flavor, coupling in self.c_leptons.items():
            if coupling == 0:
                continue

            m_lepton = MASS[flavor]
            if self.mass > 2 * m_lepton:
                channels.append(f'{flavor}_{flavor}')  # e.g., 'e_e'

        # Hadronic decays (mesons)
        # a → π⁺ π⁻
        if self.mass > 2 * MASS['pi_charged']:
            if self.c_quarks['u'] != 0 or self.c_quarks['d'] != 0:
                channels.append('pi_pi')

        # a → K⁺ K⁻
        if self.mass > 2 * MASS['K_charged']:
            if self.c_quarks['u'] != 0 or self.c_quarks['s'] != 0:
                channels.append('K_K')

        return channels

    def branching_ratios(self) -> Dict[str, float]:
        """
        Calculate branching ratios for all channels.

        Returns
        -------
        dict
            {channel_name: BR}

        Notes
        -----
        Uses leading-order partial widths from effective field theory.
        """
        if self._br_cache is not None:
            return self._br_cache

        partial_widths = {}
        m_a = self.mass

        channels = self.available_channels()

        for ch in channels:
            if ch == 'gamma_gamma':
                # Γ(a → γγ) = (g_aγγ² m_a³)/(64π)
                width = (self.g_agg**2 * m_a**3) / (64.0 * np.pi)
                partial_widths[ch] = width

            elif ch.endswith('_' + ch.split('_')[0]):
                # Leptonic: a → ℓ⁺ ℓ⁻
                flavor = ch.split('_')[0]
                m_lepton = MASS[flavor]

                if m_a <= 2 * m_lepton:
                    continue

                # Γ(a → ℓ ℓ̄) = (c_ℓ² m_a)/(8π f_a²) × (1 - 4m_ℓ²/m_a²)^(1/2)
                c_l = self.c_leptons[flavor]
                beta = np.sqrt(1 - 4*m_lepton**2 / m_a**2)

                width = (c_l**2 * m_a) / (8.0 * np.pi * self.f_a**2) * beta
                partial_widths[ch] = width

            elif ch == 'pi_pi':
                # a → π⁺ π⁻ (simplified)
                m_pi = MASS['pi_charged']
                if m_a <= 2 * m_pi:
                    continue

                # Effective coupling (mix of u and d couplings)
                c_eff = np.sqrt(self.c_quarks['u']**2 + self.c_quarks['d']**2)

                beta = np.sqrt(1 - 4*m_pi**2 / m_a**2)
                width = (c_eff**2 * m_a) / (16.0 * np.pi * self.f_a**2) * beta**3

                partial_widths[ch] = width

            elif ch == 'K_K':
                # a → K⁺ K⁻
                m_K = MASS['K_charged']
                if m_a <= 2 * m_K:
                    continue

                c_eff = np.sqrt(self.c_quarks['u']**2 + self.c_quarks['s']**2)

                beta = np.sqrt(1 - 4*m_K**2 / m_a**2)
                width = (c_eff**2 * m_a) / (16.0 * np.pi * self.f_a**2) * beta**3

                partial_widths[ch] = width

        # Normalize to BRs
        total_width = sum(partial_widths.values())
        if total_width <= 0:
            n_ch = len(channels)
            return {ch: 1.0/n_ch for ch in channels} if n_ch > 0 else {}

        self._total_width_cache = total_width
        brs = {k: v / total_width for k, v in partial_widths.items()}
        self._br_cache = brs

        return brs

    def total_width(self) -> float:
        """
        Get total decay width.

        Returns
        -------
        float
            Total width in GeV
        """
        if self._total_width_cache is not None:
            return self._total_width_cache

        self.branching_ratios()  # Computes and caches width
        return self._total_width_cache

    def sample_decay(
        self,
        parent_4vec: np.ndarray,
        channel: Optional[str] = None,
        n_events: int = 1,
        return_channel: bool = False
    ) -> Union[np.ndarray, Tuple[np.ndarray, List[str]]]:
        """
        Sample ALP decay(s).

        Parameters
        ----------
        parent_4vec : np.ndarray, shape (4,) or (n_events, 4)
            ALP 4-vector(s) in lab frame [E, px, py, pz]
        channel : str, optional
            Specific decay channel
        n_events : int
            Number of decays
        return_channel : bool
            If True, also return channel names

        Returns
        -------
        daughters : np.ndarray, shape (n_events, n_daughters, 4)
            Daughter 4-vectors
        channels : list of str (if return_channel=True)
            Channel names
        """
        parent_4vec = np.atleast_2d(parent_4vec)
        n_input = len(parent_4vec)

        if n_input == 1 and n_events > 1:
            parent_4vec = np.tile(parent_4vec, (n_events, 1))
        else:
            n_events = n_input

        # Select channels
        if channel is not None:
            channels = [channel] * n_events
        else:
            brs = self.branching_ratios()
            if not brs:
                raise ValueError("No decay channels available")

            channel_names = list(brs.keys())
            probs = [brs[c] for c in channel_names]
            indices = self._rng.choice(len(channel_names), size=n_events, p=probs)
            channels = [channel_names[i] for i in indices]

        # Sample decays
        all_daughters = []

        for i in range(n_events):
            ch = channels[i]

            # Get daughter masses
            if ch == 'gamma_gamma':
                m1, m2 = 0.0, 0.0
            elif '_' in ch:
                parts = ch.split('_')
                if parts[0] == parts[1]:  # Lepton pair
                    m1 = m2 = MASS[parts[0]]
                elif ch == 'pi_pi':
                    m1 = m2 = MASS['pi_charged']
                elif ch == 'K_K':
                    m1 = m2 = MASS['K_charged']
                else:
                    raise ValueError(f"Unknown channel: {ch}")
            else:
                raise ValueError(f"Unknown channel: {ch}")

            # Two-body decay (all ALP decays are 2-body)
            daughters_rest = sample_two_body_decay(
                self.mass, m1, m2,
                n_events=1,
                polarization=0.0,  # Scalar particle → isotropic
                rng=self._rng
            )[0]

            # Boost to lab
            daughters_lab = boost_to_lab(daughters_rest, parent_4vec[i])
            all_daughters.append(daughters_lab)

        result = np.array(all_daughters)

        if return_channel:
            return result, channels
        return result

    def get_daughter_pdgs(self, channel: str) -> List[int]:
        """Get PDG IDs for decay channel."""
        if channel == 'gamma_gamma':
            return [PDG['photon'], PDG['photon']]
        elif channel in ['e_e', 'mu_mu', 'tau_tau']:
            flavor = channel.split('_')[0]
            return [PDG[flavor], -PDG[flavor]]  # ℓ⁺ ℓ⁻
        elif channel == 'pi_pi':
            return [PDG['pi_charged'], -PDG['pi_charged']]
        elif channel == 'K_K':
            return [PDG['K_charged'], -PDG['K_charged']]
        else:
            raise ValueError(f"Unknown channel: {channel}")

    def get_charged_count(self, channel: str) -> int:
        """Get number of charged particles."""
        if channel == 'gamma_gamma':
            return 0  # Photons are neutral
        else:
            return 2  # All other channels have 2 charged particles

    def __repr__(self) -> str:
        return (f"ALP(m={self.mass:.3f} GeV, g_aγγ={self.g_agg:.2e}, "
                f"f_a={self.f_a:.2e} GeV)")
