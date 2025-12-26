"""
Advanced features for llpdecay.

Includes:
- Hadronic form factors for meson transitions
- QED radiative corrections
- Electroweak higher-order effects
- Running couplings
"""

from .form_factors import (
    form_factor_pion,
    form_factor_kaon,
    scalar_form_factor_pion,
    tensor_form_factor_pion,
    qed_correction_lepton_pair,
    qed_correction_photon_pair,
    running_alpha_em,
    electroweak_correction_hnl,
    coulomb_correction_decay,
    full_decay_correction,
    apply_form_factor_to_width
)

__all__ = [
    'form_factor_pion',
    'form_factor_kaon',
    'scalar_form_factor_pion',
    'tensor_form_factor_pion',
    'qed_correction_lepton_pair',
    'qed_correction_photon_pair',
    'running_alpha_em',
    'electroweak_correction_hnl',
    'coulomb_correction_decay',
    'full_decay_correction',
    'apply_form_factor_to_width',
]
