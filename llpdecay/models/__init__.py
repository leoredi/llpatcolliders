"""
LLP model implementations.

Provides concrete implementations of various LLP models:
- HNL: Heavy Neutral Lepton (Dirac or Majorana)
- ALP: Axion-Like Particle (Phase 2)
"""

from .base import LLPModel
from .hnl import HNL

__all__ = [
    'LLPModel',
    'HNL',
]
