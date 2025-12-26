"""
Decay channel definitions and sampling for LLPs.
"""

from .channels import (
    DecayChannel,
    HNL_CHANNELS,
    get_channel,
    list_channels,
    group_by_flavor
)

__all__ = [
    'DecayChannel',
    'HNL_CHANNELS',
    'get_channel',
    'list_channels',
    'group_by_flavor',
]
