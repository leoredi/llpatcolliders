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

from .three_body import (
    ThreeBodyPhaseSpace,
    sample_three_body_decay,
    hnl_three_body_leptonic_me
)

__all__ = [
    'DecayChannel',
    'HNL_CHANNELS',
    'get_channel',
    'list_channels',
    'group_by_flavor',
    'ThreeBodyPhaseSpace',
    'sample_three_body_decay',
    'hnl_three_body_leptonic_me',
]
