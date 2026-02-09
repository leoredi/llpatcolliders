#!/usr/bin/env python3
"""
Cross-Section Validation for HNL EW Production

Validates MadGraph cross-sections against expected values from literature.

Expected cross-sections for pp → W/Z → ℓ N at √s = 14 TeV, |V_ℓN|² = 1:
- σ(W → ℓN) ≈ 10-20 nb (10,000-20,000 pb) for m_HNL = 15 GeV
- Decreases as m_HNL approaches m_W (80 GeV) due to phase space
- Reference: arXiv:1805.08567, Fig. 12 and Table 2

Usage:
    python tools/madgraph/validate_xsec.py production/madgraph_production/summary_HNL_EW_production.csv
    python tools/madgraph/validate_xsec.py production/madgraph_production/summary_HNL_EW_production.csv --warn-threshold 0.5
"""

import sys
import argparse
import pandas as pd
import numpy as np
from pathlib import Path


# Reference cross-sections from arXiv:1805.08567 (Bondarenko et al.)
# σ(pp → W → ℓN) at 13 TeV with |V|² = 1
# Scaled to 14 TeV by luminosity ratio
REFERENCE_XSEC = {
    # mass_GeV: (xsec_pb_min, xsec_pb_max)
    5:  (15000, 25000),   # Very low mass, large phase space
    10: (12000, 20000),   # Low mass regime
    15: (10000, 18000),   # Reference point from validation review
    20: (8000,  15000),   # Mid-mass
    30: (5000,  10000),   # Higher mass, phase space closing
    40: (3000,  7000),    # Approaching W mass
    50: (2000,  5000),    # Phase space suppression
    60: (1000,  3000),    # Near kinematic limit
    70: (500,   2000),    # Very suppressed
    80: (100,   1000),    # At W mass threshold
}


def interpolate_expected_range(mass_gev):
    """
    Interpolate expected cross-section range for given mass

    Args:
        mass_gev: HNL mass in GeV

    Returns:
        tuple: (xsec_min_pb, xsec_max_pb)
    """
    masses = sorted(REFERENCE_XSEC.keys())

    # Exact match
    if mass_gev in REFERENCE_XSEC:
        return REFERENCE_XSEC[mass_gev]

    # Below lowest mass
    if mass_gev < masses[0]:
        return REFERENCE_XSEC[masses[0]]

    # Above highest mass
    if mass_gev > masses[-1]:
        return REFERENCE_XSEC[masses[-1]]

    # Linear interpolation in log space
    for i in range(len(masses) - 1):
        m1, m2 = masses[i], masses[i + 1]
        if m1 <= mass_gev <= m2:
            xsec1_min, xsec1_max = REFERENCE_XSEC[m1]
            xsec2_min, xsec2_max = REFERENCE_XSEC[m2]

            # Linear interpolation
            frac = (mass_gev - m1) / (m2 - m1)
            xsec_min = xsec1_min + frac * (xsec2_min - xsec1_min)
            xsec_max = xsec1_max + frac * (xsec2_max - xsec1_max)

            return (xsec_min, xsec_max)

    # Fallback
    return (1000, 20000)


def validate_summary_csv(csv_path, warn_threshold=0.3):
    """
    Validate cross-sections in summary CSV

    Args:
        csv_path: Path to summary CSV file
        warn_threshold: Fraction outside expected range to trigger warning

    Returns:
        dict: Validation results
    """
    csv_path = Path(csv_path)
    if not csv_path.exists():
        print(f"ERROR: Summary CSV not found: {csv_path}")
        return None

    # Read summary CSV
    df = pd.read_csv(csv_path)

    print("="*80)
    print("HNL Cross-Section Validation")
    print("="*80)
    print(f"Summary CSV: {csv_path}")
    print(f"Total entries: {len(df)}")
    print()

    # Validate each row
    results = []
    warnings = []
    errors = []

    for idx, row in df.iterrows():
        mass = row['mass_hnl_GeV']
        flavour = row['flavour']
        xsec_pb = row['xsec_pb']
        k_factor = row.get('k_factor', 1.0)

        # Get expected range
        xsec_min, xsec_max = interpolate_expected_range(mass)

        # Compare
        in_range = xsec_min <= xsec_pb <= xsec_max
        ratio_to_min = xsec_pb / xsec_min if xsec_min > 0 else np.nan
        ratio_to_max = xsec_pb / xsec_max if xsec_max > 0 else np.nan

        result = {
            'mass': mass,
            'flavour': flavour,
            'xsec_pb': xsec_pb,
            'xsec_min': xsec_min,
            'xsec_max': xsec_max,
            'in_range': in_range,
            'ratio_to_min': ratio_to_min,
            'ratio_to_max': ratio_to_max,
        }
        results.append(result)

        # Check for issues
        if xsec_pb < xsec_min * warn_threshold:
            errors.append(f"m={mass:.1f} GeV ({flavour}): σ={xsec_pb:.1e} pb << expected [{xsec_min:.1e}, {xsec_max:.1e}] pb")
        elif not in_range:
            if xsec_pb < xsec_min:
                warnings.append(f"m={mass:.1f} GeV ({flavour}): σ={xsec_pb:.1e} pb < expected [{xsec_min:.1e}, {xsec_max:.1e}] pb")
            else:
                warnings.append(f"m={mass:.1f} GeV ({flavour}): σ={xsec_pb:.1e} pb > expected [{xsec_min:.1e}, {xsec_max:.1e}] pb")

    # Summary statistics
    results_df = pd.DataFrame(results)
    n_total = len(results_df)
    n_in_range = results_df['in_range'].sum()
    n_out_of_range = n_total - n_in_range

    print("VALIDATION SUMMARY")
    print("-"*80)
    print(f"In expected range:     {n_in_range}/{n_total} ({100*n_in_range/n_total:.1f}%)")
    print(f"Out of expected range: {n_out_of_range}/{n_total} ({100*n_out_of_range/n_total:.1f}%)")
    print()

    # Print warnings
    if warnings:
        print("WARNINGS:")
        print("-"*80)
        for w in warnings[:10]:  # Limit to 10
            print(f"  ⚠ {w}")
        if len(warnings) > 10:
            print(f"  ... and {len(warnings)-10} more")
        print()

    # Print errors
    if errors:
        print("ERRORS (cross-section too low):")
        print("-"*80)
        for e in errors:
            print(f"  ✗ {e}")
        print()

    # Detailed table for out-of-range points
    if n_out_of_range > 0 and n_out_of_range <= 20:
        print("OUT-OF-RANGE DETAILS:")
        print("-"*80)
        out_of_range = results_df[~results_df['in_range']]
        for _, row in out_of_range.iterrows():
            print(f"  m={row['mass']:.1f} GeV ({row['flavour']}): "
                  f"σ={row['xsec_pb']:.2e} pb, "
                  f"expected [{row['xsec_min']:.2e}, {row['xsec_max']:.2e}] pb, "
                  f"ratio={row['ratio_to_min']:.2f}×min, {row['ratio_to_max']:.2f}×max")
        print()

    # Check for common issues
    print("DIAGNOSTIC CHECKS:")
    print("-"*80)

    # Check 1: Are all cross-sections very low?
    median_xsec = results_df['xsec_pb'].median()
    if median_xsec < 100:
        print(f"  ✗ Median cross-section is very low ({median_xsec:.2e} pb)")
        print(f"    → Check param_card mixing parameters (should be V=1.0, not V²=1.0)")
    elif median_xsec > 100000:
        print(f"  ✗ Median cross-section is very high ({median_xsec:.2e} pb)")
        print(f"    → Check param_card mixing parameters (might be V²=1.0 instead of V=1.0)")
    else:
        print(f"  ✓ Median cross-section is reasonable ({median_xsec:.2e} pb)")

    # Check 2: Is K-factor applied?
    if 'k_factor' in df.columns:
        k_factors = df['k_factor'].unique()
        if len(k_factors) == 1 and k_factors[0] == 1.3:
            print(f"  ✓ K-factor applied (K={k_factors[0]})")
        else:
            print(f"  ⚠ K-factor not standard: {k_factors}")
    else:
        print(f"  ⚠ K-factor column missing")

    # Check 3: Mass dependence
    if len(results_df) > 5:
        mass_sorted = results_df.sort_values('mass')
        xsec_trend = mass_sorted['xsec_pb'].values
        masses_trend = mass_sorted['mass'].values

        # Cross-section should decrease with mass
        correlation = np.corrcoef(masses_trend, xsec_trend)[0, 1]
        if correlation < -0.5:
            print(f"  ✓ Cross-section decreases with mass (correlation={correlation:.2f})")
        else:
            print(f"  ⚠ Cross-section trend unusual (correlation={correlation:.2f})")

    print()
    print("="*80)

    # Final verdict
    if errors:
        print("STATUS: ✗ FAILED - Some cross-sections are too low")
        print("        Check param_card mixing parameters and process definition")
        return_code = 2
    elif n_out_of_range > n_total * 0.3:
        print(f"STATUS: ⚠ WARNING - {n_out_of_range}/{n_total} cross-sections out of expected range")
        print("        This may be acceptable depending on PDF/scale uncertainties")
        return_code = 1
    else:
        print("STATUS: ✓ PASSED - Cross-sections in expected range")
        return_code = 0

    print("="*80)

    return {
        'results': results_df,
        'n_in_range': n_in_range,
        'n_out_of_range': n_out_of_range,
        'warnings': warnings,
        'errors': errors,
        'return_code': return_code,
    }


def main():
    parser = argparse.ArgumentParser(
        description='Validate HNL cross-sections against expected values'
    )
    parser.add_argument(
        'csv_file',
        help='Path to summary CSV file'
    )
    parser.add_argument(
        '--warn-threshold',
        type=float,
        default=0.3,
        help='Fraction of expected minimum to trigger error (default: 0.3)'
    )

    args = parser.parse_args()

    result = validate_summary_csv(args.csv_file, args.warn_threshold)

    if result is None:
        return 1

    return result['return_code']


if __name__ == '__main__':
    sys.exit(main())
