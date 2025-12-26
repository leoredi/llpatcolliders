"""
Heavy Neutral Lepton (HNL) model implementation.

Implements Dirac or Majorana HNL with configurable mixing angles.
Supports production-independent decay sampling with proper matrix elements.

References:
- Bondarenko et al., arXiv:1805.08567 (JHEP 2019)
- Gorbunov & Shaposhnikov, arXiv:hep-ph/9911364
- Atre et al., arXiv:0901.3589
"""

import numpy as np
from typing import Dict, Optional, Tuple, List, Union
import warnings

from .base import LLPModel
from ..core import (
    G_FERMI, F_PI, F_K,
    V_UD, V_US,
    sample_two_body_decay,
    boost_to_lab,
    four_vector_mass,
    kallen
)
from ..decays import HNL_CHANNELS, DecayChannel


class HNL(LLPModel):
    """
    Heavy Neutral Lepton model.

    Implements production-independent decay sampling with proper kinematics
    and matrix elements. Can compute branching ratios internally or import
    from HNLCalc.

    Parameters
    ----------
    mass : float
        HNL mass in GeV
    Ue : float, default=0
        Mixing angle squared |U_e|²
    Umu : float, default=0
        Mixing angle squared |U_μ|²
    Utau : float, default=0
        Mixing angle squared |U_τ|²
    is_majorana : bool, default=True
        If True, include Majorana modes (LNV decays)
    use_simple_brs : bool, default=True
        If True, use simplified BR formulas. If False, try to use HNLCalc.
    seed : int, optional
        Random seed for reproducibility

    Attributes
    ----------
    mass : float
        HNL mass in GeV
    Ue, Umu, Utau : float
        Mixing angles squared
    is_majorana : bool
        Majorana vs Dirac nature

    Examples
    --------
    >>> # Muon-mixing dominated HNL
    >>> hnl = HNL(mass=2.0, Ue=0, Umu=1e-6, Utau=0)
    >>> print(hnl.available_channels())
    ['e_pi', 'mu_pi', 'e_K', 'mu_K', ...]

    >>> # Sample a decay
    >>> parent_4vec = np.array([10.0, 3.0, 0.5, 9.3])
    >>> daughters, channel = hnl.sample_decay(parent_4vec, return_channel=True)
    >>> print(f"Decayed via {channel}: {daughters.shape}")

    >>> # Get branching ratios
    >>> brs = hnl.branching_ratios()
    >>> for ch, br in sorted(brs.items(), key=lambda x: -x[1])[:5]:
    ...     print(f"{ch}: {br:.3f}")
    """

    def __init__(
        self,
        mass: float,
        Ue: float = 0.0,
        Umu: float = 0.0,
        Utau: float = 0.0,
        is_majorana: bool = True,
        use_simple_brs: bool = True,
        seed: Optional[int] = None
    ):
        super().__init__(mass)

        self.Ue = Ue
        self.Umu = Umu
        self.Utau = Utau
        self.is_majorana = is_majorana
        self.use_simple_brs = use_simple_brs

        if self.U2_total <= 0:
            warnings.warn("Total mixing is zero - HNL is stable!")

        self._br_cache: Optional[Dict[str, float]] = None
        self._total_width_cache: Optional[float] = None
        self._rng = np.random.default_rng(seed)

    @property
    def U2_total(self) -> float:
        """Total mixing |U|² = |U_e|² + |U_μ|² + |U_τ|²"""
        return self.Ue + self.Umu + self.Utau

    @property
    def flavor_fractions(self) -> Dict[str, float]:
        """Fraction of mixing per flavor."""
        total = self.U2_total
        if total == 0:
            return {'e': 0.0, 'mu': 0.0, 'tau': 0.0}
        return {
            'e': self.Ue / total,
            'mu': self.Umu / total,
            'tau': self.Utau / total,
        }

    def available_channels(self) -> List[str]:
        """
        List decay channels kinematically accessible at this mass.

        Returns
        -------
        list of str
            Channel names where m_N > sum(m_daughters)
        """
        channels = []
        for name, ch in HNL_CHANNELS.items():
            # Check kinematic threshold
            if self.mass > ch.threshold_mass():
                # For 3-body, only include if not too heavy
                if ch.n_body == 3 and self.mass > 5.0:
                    # 3-body becomes subleading at high mass
                    continue
                channels.append(name)
        return channels

    def branching_ratios(self, use_hnlcalc: bool = False) -> Dict[str, float]:
        """
        Get branching ratios for all accessible channels.

        Parameters
        ----------
        use_hnlcalc : bool
            If True, import from HNLCalc. If False, use internal calculation.

        Returns
        -------
        dict
            {channel_name: BR} for accessible channels, normalized to sum to 1
        """
        if self._br_cache is not None and not use_hnlcalc:
            return self._br_cache

        if use_hnlcalc and not self.use_simple_brs:
            self._br_cache = self._get_brs_from_hnlcalc()
        else:
            self._br_cache = self._calculate_brs_internal()

        return self._br_cache

    def _calculate_brs_internal(self) -> Dict[str, float]:
        """
        Calculate branching ratios using simplified formulas.

        Uses leading-order matrix elements from literature.
        Good approximation for 0.1 < m_N < 5 GeV.

        References
        ----------
        - Bondarenko et al., arXiv:1805.08567, Eqs. (3.1)-(3.8)
        """
        partial_widths = {}
        m_N = self.mass

        available = self.available_channels()

        for ch_name in available:
            ch = HNL_CHANNELS[ch_name]

            if ch.n_body == 2:
                # Two-body decay
                width = self._partial_width_two_body(ch)
                if width > 0:
                    partial_widths[ch_name] = width

            elif ch.n_body == 3:
                # Three-body decay
                width = self._partial_width_three_body(ch)
                if width > 0:
                    partial_widths[ch_name] = width

        # Normalize to branching ratios
        total_width = sum(partial_widths.values())
        if total_width <= 0:
            # No decays available - return uniform (shouldn't happen)
            n_channels = len(available)
            return {ch: 1.0/n_channels for ch in available} if n_channels > 0 else {}

        # Cache total width
        self._total_width_cache = total_width

        brs = {k: v / total_width for k, v in partial_widths.items()}
        return brs

    def _partial_width_two_body(self, channel: DecayChannel) -> float:
        """
        Calculate partial width for two-body decay.

        Γ(N → ℓ M) where M is a meson (π, K, ρ, ...)

        Formula from Bondarenko et al., arXiv:1805.08567, Eq. (3.1):
        Γ = (G_F² f_M² |V_qq'|² |U_α|² m_N³) / (16π) × λ^(1/2) × (1 - m_ℓ²/m_N²)²

        where λ = λ(m_N², m_ℓ², m_M²) is the Källén function.
        """
        m_N = self.mass
        m1, m2 = channel.daughter_masses()

        # Determine which mixing applies
        U2 = self._get_mixing_for_channel(channel)
        if U2 <= 0:
            return 0.0

        # Determine decay constant and CKM element
        f_M, V_ckm = self._get_meson_params(channel)

        # Phase space factor
        lam = kallen(m_N**2, m1**2, m2**2)
        if lam <= 0:
            return 0.0

        # For charged current: N → ℓ⁻ M⁺
        # Determine lepton and meson
        if channel.name.startswith('nu_'):
            # Neutral current N → ν M⁰
            # Different formula - suppressed by factor
            # Approximate as ~0.5 × charged current
            is_nc = True
        else:
            is_nc = False

        # Width formula
        if is_nc:
            # Neutral current (simplified)
            prefactor = (G_FERMI**2 * f_M**2 * U2 * m_N**3) / (32.0 * np.pi)
        else:
            # Charged current
            prefactor = (G_FERMI**2 * f_M**2 * V_ckm**2 * U2 * m_N**3) / (16.0 * np.pi)

        # Phase space and helicity suppression
        ps_factor = np.sqrt(lam) / m_N**2

        # Helicity suppression for lepton (when applicable)
        if not is_nc:
            # m1 is the lepton mass
            helicity_factor = (1.0 - m1**2 / m_N**2)**2
        else:
            helicity_factor = 1.0

        width = prefactor * ps_factor * helicity_factor

        # Majorana factor: can decay to ℓ⁺ and ℓ⁻ (doubles rate)
        if self.is_majorana and not is_nc:
            width *= 2.0

        return width

    def _partial_width_three_body(self, channel: DecayChannel) -> float:
        """
        Calculate partial width for three-body decay.

        Γ(N → ν ℓ ℓ') proceeds via virtual W* or Z*.

        Simplified formula (from Gorbunov & Shaposhnikov):
        Γ ≈ (G_F² m_N⁵ |U|²) / (192 π³) × C

        where C is a channel-dependent coefficient (~1 for most modes).
        """
        m_N = self.mass

        # Get mixing
        U2 = self._get_mixing_for_channel(channel)
        if U2 <= 0:
            return 0.0

        # Three-body phase space (simplified)
        # From Gorbunov & Shaposhnikov, arXiv:hep-ph/9911364
        C_factor = 1.0  # Order unity for most channels

        width = (G_FERMI**2 * m_N**5 * U2) / (192.0 * np.pi**3) * C_factor

        # Majorana doubles the rate
        if self.is_majorana:
            width *= 2.0

        return width

    def _get_mixing_for_channel(self, channel: DecayChannel) -> float:
        """Determine which mixing angle applies to this channel."""
        name = channel.name

        # Charged current: flavor determined by lepton
        if name.startswith('e_'):
            return self.Ue
        elif name.startswith('mu_'):
            return self.Umu
        elif name.startswith('tau_'):
            return self.Utau

        # Neutral current: sum over flavors
        elif name.startswith('nu_'):
            # For NC, mixing depends on which neutrino
            if 'nu_e' in name:
                return self.Ue
            elif 'nu_mu' in name:
                return self.Umu
            elif 'nu_tau' in name:
                return self.Utau
            else:
                # Mixed or summed
                return self.U2_total

        else:
            # Default: use total mixing
            return self.U2_total

    def _get_meson_params(self, channel: DecayChannel) -> Tuple[float, float]:
        """
        Get meson decay constant and CKM element for channel.

        Returns
        -------
        f_M : float
            Meson decay constant in GeV
        V_ckm : float
            Relevant CKM matrix element magnitude
        """
        name = channel.name

        if 'pi' in name:
            return F_PI, V_UD
        elif 'K' in name:
            return F_K, V_US
        elif 'rho' in name:
            # ρ meson decay constant (approximate)
            return 0.220, V_UD  # f_ρ ≈ 220 MeV
        else:
            # Default to pion
            return F_PI, V_UD

    def _get_brs_from_hnlcalc(self) -> Dict[str, float]:
        """
        Import branching ratios from HNLCalc.

        TODO: Interface with analysis_pbc/HNLCalc
        """
        raise NotImplementedError(
            "HNLCalc integration not yet implemented. "
            "Set use_simple_brs=True to use internal calculation."
        )

    def total_width(self) -> float:
        """
        Get total decay width in GeV.

        Returns
        -------
        float
            Total width Γ_total in GeV
        """
        if self._total_width_cache is not None:
            return self._total_width_cache

        # Compute BRs (which caches total width)
        self.branching_ratios()

        return self._total_width_cache

    def sample_decay(
        self,
        parent_4vec: np.ndarray,
        channel: Optional[str] = None,
        polarization: float = 0.0,
        n_events: int = 1,
        return_channel: bool = False
    ) -> Union[np.ndarray, Tuple[np.ndarray, List[str]]]:
        """
        Sample HNL decay(s) and boost to lab frame.

        Parameters
        ----------
        parent_4vec : np.ndarray, shape (4,) or (n_events, 4)
            HNL 4-vector(s) in lab frame [E, px, py, pz]
        channel : str, optional
            Specific decay channel. If None, sample according to BRs.
        polarization : float, default=0.0
            HNL polarization along momentum direction, P ∈ [-1, 1]
        n_events : int, default=1
            Number of decays to sample (if parent_4vec is single event)
        return_channel : bool, default=False
            If True, also return the sampled channel name(s)

        Returns
        -------
        daughters : np.ndarray, shape (n_events, n_daughters, 4)
            Daughter 4-vectors in lab frame
        channels : list of str (only if return_channel=True)
            Channel name for each event

        Examples
        --------
        >>> hnl = HNL(mass=2.0, Umu=1e-6)
        >>> parent = np.array([10.0, 3.0, 0.0, 9.5])
        >>> daughters, ch = hnl.sample_decay(parent, return_channel=True)
        >>> print(f"Channel: {ch}, daughters shape: {daughters.shape}")
        """
        parent_4vec = np.atleast_2d(parent_4vec)
        n_input = len(parent_4vec)

        if n_input == 1 and n_events > 1:
            parent_4vec = np.tile(parent_4vec, (n_events, 1))
        else:
            n_events = n_input

        # Select channel(s)
        if channel is not None:
            channels = [channel] * n_events
        else:
            brs = self.branching_ratios()
            channel_names = list(brs.keys())
            if len(channel_names) == 0:
                raise ValueError(f"No decay channels available for m_N = {self.mass} GeV")
            channel_probs = [brs[c] for c in channel_names]
            channel_indices = self._rng.choice(
                len(channel_names), size=n_events, p=channel_probs
            )
            channels = [channel_names[i] for i in channel_indices]

        # Sample decays
        all_daughters = []
        for i in range(n_events):
            ch = HNL_CHANNELS[channels[i]]

            if ch.n_body == 2:
                m1, m2 = ch.daughter_masses()
                daughters_rest = sample_two_body_decay(
                    self.mass, m1, m2,
                    n_events=1,
                    polarization=polarization,
                    alpha=ch.asymmetry_param,
                    rng=self._rng
                )[0]  # Shape: (2, 4)

                # Boost to lab frame
                daughters_lab = boost_to_lab(daughters_rest, parent_4vec[i])
                all_daughters.append(daughters_lab)

            elif ch.n_body == 3:
                # TODO: Implement 3-body phase space sampling (Phase 2)
                raise NotImplementedError(
                    f"3-body decay {ch.name} not yet implemented. "
                    "This is planned for Phase 2."
                )

        result = np.array(all_daughters)  # Shape: (n_events, n_daughters, 4)

        if return_channel:
            return result, channels
        return result

    def get_daughter_pdgs(self, channel: str) -> List[int]:
        """Get PDG IDs of daughters for a channel."""
        if channel not in HNL_CHANNELS:
            raise KeyError(f"Unknown channel: {channel}")
        return HNL_CHANNELS[channel].daughter_pdgs()

    def get_charged_count(self, channel: str) -> int:
        """Get number of charged particles in decay."""
        if channel not in HNL_CHANNELS:
            raise KeyError(f"Unknown channel: {channel}")
        return HNL_CHANNELS[channel].charged_count

    def __repr__(self) -> str:
        flavor_str = f"Ue²={self.Ue:.1e}, Uμ²={self.Umu:.1e}, Uτ²={self.Utau:.1e}"
        nature = "Majorana" if self.is_majorana else "Dirac"
        return f"HNL(m={self.mass:.3f} GeV, {flavor_str}, {nature})"
