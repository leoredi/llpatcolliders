"""
LLP model implementations.

Provides concrete implementations of various LLP models:
- HNL: Heavy Neutral Lepton (Dirac or Majorana)
- ALP: Axion-Like Particle
"""

from .base import LLPModel
from .hnl import HNL
from .alp import ALP

__all__ = [
    'LLPModel',
    'HNL',
    'ALP',
]
