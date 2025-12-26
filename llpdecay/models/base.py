"""
Abstract base class for LLP models.

Defines the interface that all LLP decay models must implement.
"""

from abc import ABC, abstractmethod
from typing import Dict, List, Optional, Tuple
import numpy as np


class LLPModel(ABC):
    """
    Abstract base class for Long-Lived Particle models.

    All concrete LLP implementations (HNL, ALP, etc.) should inherit from this
    class and implement the required methods.

    Attributes
    ----------
    mass : float
        LLP mass in GeV

    Methods to Implement
    --------------------
    available_channels() : List[str]
        Return list of kinematically accessible decay channels

    branching_ratios() : Dict[str, float]
        Return branching ratios for all channels

    sample_decay(parent_4vec, ...) : np.ndarray
        Sample decay and return daughter 4-vectors

    get_daughter_pdgs(channel) : List[int]
        Get PDG IDs of daughters for a channel

    get_charged_count(channel) : int
        Get number of charged particles in a channel
    """

    def __init__(self, mass: float):
        """
        Initialize LLP model.

        Parameters
        ----------
        mass : float
            LLP mass in GeV
        """
        if mass <= 0:
            raise ValueError(f"Mass must be positive, got {mass}")
        self.mass = mass

    @abstractmethod
    def available_channels(self) -> List[str]:
        """
        List decay channels kinematically accessible at this mass.

        Returns
        -------
        list of str
            Channel names that satisfy m_LLP > sum(m_daughters)
        """
        pass

    @abstractmethod
    def branching_ratios(self) -> Dict[str, float]:
        """
        Get branching ratios for all accessible channels.

        Returns
        -------
        dict
            {channel_name: BR} normalized to sum to 1.0
        """
        pass

    @abstractmethod
    def sample_decay(
        self,
        parent_4vec: np.ndarray,
        channel: Optional[str] = None,
        **kwargs
    ) -> np.ndarray:
        """
        Sample LLP decay(s) and boost to lab frame.

        Parameters
        ----------
        parent_4vec : np.ndarray, shape (4,) or (n_events, 4)
            LLP 4-vector(s) in lab frame [E, px, py, pz]
        channel : str, optional
            Specific decay channel. If None, sample according to BRs.
        **kwargs
            Additional model-specific parameters

        Returns
        -------
        np.ndarray, shape (n_events, n_daughters, 4)
            Daughter 4-vectors in lab frame
        """
        pass

    @abstractmethod
    def get_daughter_pdgs(self, channel: str) -> List[int]:
        """
        Get PDG IDs of daughter particles for a channel.

        Parameters
        ----------
        channel : str
            Decay channel name

        Returns
        -------
        list of int
            PDG IDs (negative for antiparticles)
        """
        pass

    @abstractmethod
    def get_charged_count(self, channel: str) -> int:
        """
        Get number of electrically charged particles in decay.

        Parameters
        ----------
        channel : str
            Decay channel name

        Returns
        -------
        int
            Number of charged particles (0 = invisible)
        """
        pass

    def total_width(self) -> float:
        """
        Get total decay width in GeV.

        Default implementation returns None. Override if model computes width.

        Returns
        -------
        float or None
            Total decay width in GeV
        """
        return None

    def lifetime(self) -> float:
        """
        Get mean lifetime in seconds.

        Requires total_width() to be implemented.

        Returns
        -------
        float or None
            Mean lifetime τ = ℏ/Γ in seconds

        Raises
        ------
        NotImplementedError
            If total_width() not implemented
        """
        from ..core import HBAR_GEV_S

        width = self.total_width()
        if width is None:
            raise NotImplementedError("total_width() must be implemented to compute lifetime")
        if width <= 0:
            return float('inf')  # Stable particle
        return HBAR_GEV_S / width

    def ctau(self) -> float:
        """
        Get mean decay length c·τ in meters.

        Requires total_width() to be implemented.

        Returns
        -------
        float or None
            Mean decay length in meters

        Raises
        ------
        NotImplementedError
            If total_width() not implemented
        """
        from ..core import HBAR_C_GEV_M

        width = self.total_width()
        if width is None:
            raise NotImplementedError("total_width() must be implemented to compute ctau")
        if width <= 0:
            return float('inf')
        return HBAR_C_GEV_M / width

    def decay_probability(
        self,
        parent_4vec: np.ndarray,
        x_min: float,
        x_max: float
    ) -> np.ndarray:
        """
        Calculate probability of decay in interval [x_min, x_max].

        P(x_min < x < x_max) = exp(-x_min/(γβcτ)) - exp(-x_max/(γβcτ))

        Parameters
        ----------
        parent_4vec : np.ndarray, shape (4,) or (n_events, 4)
            LLP 4-vector(s) [E, px, py, pz]
        x_min : float
            Start of decay region (meters)
        x_max : float
            End of decay region (meters)

        Returns
        -------
        np.ndarray
            Decay probability for each event

        Raises
        ------
        NotImplementedError
            If total_width() not implemented
        """
        parent_4vec = np.atleast_2d(parent_4vec)

        # Lorentz boost factor γβ = p/m
        E = parent_4vec[:, 0]
        p = np.sqrt(parent_4vec[:, 1]**2 + parent_4vec[:, 2]**2 + parent_4vec[:, 3]**2)
        gamma_beta = p / self.mass

        # Decay length in lab frame
        ctau_lab = self.ctau() * gamma_beta

        # Probability to decay in [x_min, x_max]
        P = np.exp(-x_min / ctau_lab) - np.exp(-x_max / ctau_lab)

        return P.squeeze() if len(parent_4vec) == 1 else P

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(mass={self.mass:.3f} GeV)"
