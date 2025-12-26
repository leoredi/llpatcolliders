"""
Validation tools for llpdecay.

Provides interfaces to external codes and validation utilities:
- HNLCalc integration for BR cross-validation
- Transverse detector geometry (MATHUSLA, ANUBIS)
- Comparison tools and reporting
"""

from .hnlcalc_interface import (
    get_hnlcalc_branching_ratios,
    compare_branching_ratios,
    validate_total_width,
    print_comparison_table,
    HNLCalcValidator
)

from .transverse_geometry import (
    TransverseDetectorGeometry,
    MATHUSLA,
    ANUBIS,
    decay_position,
    geometric_acceptance,
    decay_probability_transverse,
    signal_yield_transverse
)

__all__ = [
    'get_hnlcalc_branching_ratios',
    'compare_branching_ratios',
    'validate_total_width',
    'print_comparison_table',
    'HNLCalcValidator',
    'TransverseDetectorGeometry',
    'MATHUSLA',
    'ANUBIS',
    'decay_position',
    'geometric_acceptance',
    'decay_probability_transverse',
    'signal_yield_transverse',
]
