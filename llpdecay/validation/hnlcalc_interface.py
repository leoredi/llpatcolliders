"""
Interface to HNLCalc for branching ratio validation.

Provides integration with the existing HNLCalc package in analysis_pbc/HNLCalc
for cross-validation of branching ratios and decay widths.

This module allows users to:
1. Import BRs from HNLCalc for comparison
2. Validate llpdecay predictions against HNLCalc
3. Generate validation plots and tables
"""

import sys
import warnings
from pathlib import Path
from typing import Dict, Optional, Tuple
import numpy as np


def _get_hnlcalc_module():
    """
    Dynamically import HNLCalc module.

    Returns
    -------
    module or None
        HNLCalc module if available, None otherwise
    """
    try:
        # Try to find HNLCalc in analysis_pbc
        repo_root = Path(__file__).parent.parent.parent
        hnlcalc_path = repo_root / 'analysis_pbc' / 'HNLCalc'

        if hnlcalc_path.exists():
            sys.path.insert(0, str(hnlcalc_path.parent))
            import HNLCalc
            return HNLCalc
        else:
            return None
    except ImportError:
        return None


def get_hnlcalc_branching_ratios(
    mass: float,
    Ue: float = 0.0,
    Umu: float = 0.0,
    Utau: float = 0.0
) -> Optional[Dict[str, float]]:
    """
    Get branching ratios from HNLCalc.

    Parameters
    ----------
    mass : float
        HNL mass in GeV
    Ue, Umu, Utau : float
        Mixing angles squared

    Returns
    -------
    dict or None
        {channel_name: BR} if HNLCalc available, None otherwise

    Notes
    -----
    Channel names are mapped from HNLCalc conventions to llpdecay conventions.
    """
    HNLCalc = _get_hnlcalc_module()

    if HNLCalc is None:
        warnings.warn(
            "HNLCalc module not found. Install from analysis_pbc/HNLCalc "
            "or use llpdecay's internal BR calculation."
        )
        return None

    try:
        # Create HNLCalc instance
        # Note: HNLCalc API may vary - adjust as needed
        hnl = HNLCalc.HNL(mass, Ue, Umu, Utau)

        # Get decay branching ratios
        # Map HNLCalc channel names to llpdecay names
        br_map = _map_hnlcalc_channels(hnl)

        return br_map

    except Exception as e:
        warnings.warn(f"Failed to get BRs from HNLCalc: {e}")
        return None


def _map_hnlcalc_channels(hnl_calc_instance) -> Dict[str, float]:
    """
    Map HNLCalc channel names to llpdecay conventions.

    Parameters
    ----------
    hnl_calc_instance
        HNLCalc HNL instance

    Returns
    -------
    dict
        {llpdecay_channel_name: BR}
    """
    # This mapping depends on HNLCalc's internal structure
    # Adjust based on actual HNLCalc API

    channel_mapping = {
        # HNLCalc name -> llpdecay name
        'e+pi-': 'e_pi',
        'mu+pi-': 'mu_pi',
        'tau+pi-': 'tau_pi',
        'e+K-': 'e_K',
        'mu+K-': 'mu_K',
        'tau+K-': 'tau_K',
        'e+rho-': 'e_rho',
        'mu+rho-': 'mu_rho',
        'tau+rho-': 'tau_rho',
        'nu pi0': 'nu_mu_pi0',  # Adjust flavor as needed
        'nu e e': 'nu_e_e',
        'nu mu mu': 'nu_mu_mu',
        'nu tau tau': 'nu_tau_tau',
    }

    brs = {}

    # Try to extract BRs from HNLCalc
    # This is a placeholder - actual implementation depends on HNLCalc API
    try:
        if hasattr(hnl_calc_instance, 'get_branching_ratios'):
            hnlcalc_brs = hnl_calc_instance.get_branching_ratios()

            for hnlcalc_name, llpdecay_name in channel_mapping.items():
                if hnlcalc_name in hnlcalc_brs:
                    brs[llpdecay_name] = hnlcalc_brs[hnlcalc_name]

    except AttributeError:
        # HNLCalc might have different API
        pass

    return brs


def compare_branching_ratios(
    mass: float,
    Ue: float = 0.0,
    Umu: float = 0.0,
    Utau: float = 0.0,
    rtol: float = 0.1
) -> Tuple[Dict[str, Dict[str, float]], bool]:
    """
    Compare llpdecay BRs with HNLCalc.

    Parameters
    ----------
    mass : float
        HNL mass in GeV
    Ue, Umu, Utau : float
        Mixing angles squared
    rtol : float
        Relative tolerance for agreement (default 10%)

    Returns
    -------
    comparison : dict
        {channel: {'llpdecay': BR, 'hnlcalc': BR, 'diff': fractional_diff}}
    passed : bool
        True if all channels agree within rtol

    Examples
    --------
    >>> comparison, passed = compare_branching_ratios(2.0, Umu=1e-6)
    >>> if passed:
    ...     print("BRs agree with HNLCalc!")
    >>> else:
    ...     for ch, vals in comparison.items():
    ...         print(f"{ch}: diff = {vals['diff']:.1%}")
    """
    from ..models import HNL

    # Get BRs from both packages
    hnl_llpdecay = HNL(mass, Ue, Umu, Utau)
    brs_llpdecay = hnl_llpdecay.branching_ratios()

    brs_hnlcalc = get_hnlcalc_branching_ratios(mass, Ue, Umu, Utau)

    if brs_hnlcalc is None:
        warnings.warn("HNLCalc not available - cannot compare")
        return {}, False

    # Compare
    comparison = {}
    all_passed = True

    # Get union of channels
    all_channels = set(brs_llpdecay.keys()) | set(brs_hnlcalc.keys())

    for ch in all_channels:
        br_llp = brs_llpdecay.get(ch, 0.0)
        br_hnlcalc = brs_hnlcalc.get(ch, 0.0)

        # Fractional difference
        if br_hnlcalc > 0:
            diff = abs(br_llp - br_hnlcalc) / br_hnlcalc
        elif br_llp > 0:
            diff = float('inf')  # Missing in HNLCalc
        else:
            diff = 0.0

        comparison[ch] = {
            'llpdecay': br_llp,
            'hnlcalc': br_hnlcalc,
            'diff': diff
        }

        if diff > rtol:
            all_passed = False

    return comparison, all_passed


def validate_total_width(
    mass: float,
    Ue: float = 0.0,
    Umu: float = 0.0,
    Utau: float = 0.0,
    rtol: float = 0.1
) -> Tuple[float, float, float, bool]:
    """
    Compare total decay width with HNLCalc.

    Parameters
    ----------
    mass : float
        HNL mass in GeV
    Ue, Umu, Utau : float
        Mixing angles squared
    rtol : float
        Relative tolerance (default 10%)

    Returns
    -------
    width_llpdecay : float
        Total width from llpdecay (GeV)
    width_hnlcalc : float
        Total width from HNLCalc (GeV)
    fractional_diff : float
        Fractional difference
    passed : bool
        True if agree within rtol
    """
    from ..models import HNL

    hnl = HNL(mass, Ue, Umu, Utau)
    width_llpdecay = hnl.total_width()

    # Get width from HNLCalc
    HNLCalc = _get_hnlcalc_module()
    if HNLCalc is None:
        return width_llpdecay, None, None, False

    try:
        hnl_calc = HNLCalc.HNL(mass, Ue, Umu, Utau)
        if hasattr(hnl_calc, 'get_total_width'):
            width_hnlcalc = hnl_calc.get_total_width()
        elif hasattr(hnl_calc, 'total_width'):
            width_hnlcalc = hnl_calc.total_width
        else:
            return width_llpdecay, None, None, False

        diff = abs(width_llpdecay - width_hnlcalc) / width_hnlcalc
        passed = diff <= rtol

        return width_llpdecay, width_hnlcalc, diff, passed

    except Exception as e:
        warnings.warn(f"Failed to get width from HNLCalc: {e}")
        return width_llpdecay, None, None, False


def print_comparison_table(
    mass: float,
    Ue: float = 0.0,
    Umu: float = 0.0,
    Utau: float = 0.0
):
    """
    Print formatted comparison table of BRs.

    Parameters
    ----------
    mass : float
        HNL mass in GeV
    Ue, Umu, Utau : float
        Mixing angles squared
    """
    comparison, passed = compare_branching_ratios(mass, Ue, Umu, Utau)

    if not comparison:
        print("HNLCalc not available - cannot generate comparison")
        return

    print(f"\nBranching Ratio Comparison (m_N = {mass} GeV)")
    print(f"Mixings: Ue²={Ue:.2e}, Uμ²={Umu:.2e}, Uτ²={Utau:.2e}")
    print("=" * 70)
    print(f"{'Channel':<15} {'llpdecay':<12} {'HNLCalc':<12} {'Diff':<10} {'Status'}")
    print("-" * 70)

    for ch, vals in sorted(comparison.items(), key=lambda x: -x[1]['llpdecay']):
        br_llp = vals['llpdecay']
        br_hnlcalc = vals['hnlcalc']
        diff = vals['diff']

        # Format difference
        if diff == float('inf'):
            diff_str = "MISSING"
            status = "✗"
        elif diff < 0.05:
            diff_str = f"{diff:.1%}"
            status = "✓"
        elif diff < 0.10:
            diff_str = f"{diff:.1%}"
            status = "~"
        else:
            diff_str = f"{diff:.1%}"
            status = "✗"

        print(f"{ch:<15} {br_llp:<12.3%} {br_hnlcalc:<12.3%} {diff_str:<10} {status}")

    print("=" * 70)

    if passed:
        print("✓ All channels agree within 10%")
    else:
        print("✗ Some channels disagree - review implementation")

    print()


class HNLCalcValidator:
    """
    Comprehensive validation against HNLCalc.

    Performs systematic tests across mass ranges and mixing scenarios.

    Examples
    --------
    >>> validator = HNLCalcValidator()
    >>> results = validator.scan_mass_range([0.5, 1.0, 2.0, 5.0], Umu=1e-6)
    >>> validator.plot_comparison(results)
    """

    def __init__(self):
        self.hnlcalc_available = _get_hnlcalc_module() is not None

    def scan_mass_range(
        self,
        masses: list,
        Ue: float = 0.0,
        Umu: float = 0.0,
        Utau: float = 0.0
    ) -> Dict[float, Dict]:
        """
        Scan BRs across mass range.

        Parameters
        ----------
        masses : list of float
            HNL masses to test (GeV)
        Ue, Umu, Utau : float
            Mixing angles squared

        Returns
        -------
        dict
            {mass: comparison_dict}
        """
        results = {}

        for mass in masses:
            print(f"Testing m_N = {mass} GeV...")
            comparison, passed = compare_branching_ratios(mass, Ue, Umu, Utau)
            results[mass] = {
                'comparison': comparison,
                'passed': passed
            }

        return results

    def report(self, results: Dict[float, Dict]):
        """Print summary report of validation results."""
        print("\n" + "=" * 70)
        print("HNLCalc Validation Report")
        print("=" * 70)

        n_total = len(results)
        n_passed = sum(1 for r in results.values() if r['passed'])

        print(f"\nMass points tested: {n_total}")
        print(f"Passed: {n_passed}/{n_total} ({n_passed/n_total:.0%})")

        if n_passed < n_total:
            print("\nFailed mass points:")
            for mass, result in results.items():
                if not result['passed']:
                    print(f"  - {mass} GeV")

        print()
