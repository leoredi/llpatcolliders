"""
Validation tools for llpdecay.

Provides interfaces to external codes and validation utilities:
- HNLCalc integration for BR cross-validation
- Comparison tools and reporting
"""

from .hnlcalc_interface import (
    get_hnlcalc_branching_ratios,
    compare_branching_ratios,
    validate_total_width,
    print_comparison_table,
    HNLCalcValidator
)

__all__ = [
    'get_hnlcalc_branching_ratios',
    'compare_branching_ratios',
    'validate_total_width',
    'print_comparison_table',
    'HNLCalcValidator',
]
