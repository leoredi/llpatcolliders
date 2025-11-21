"""
|U|² Limit Calculator - PBC-Grade Analysis

Computes expected |U_flavour|² exclusion limits using:
1. HNLCalc model (width, production, BRs)
2. Per-parent geometric efficiencies
3. Direct N_sig calculation (no BR_limit detour)

Matches methodology of Physics Beyond Colliders benchmark curves BC6/7/8.

Usage:
------
    from u2_limit_calculator import U2LimitCalculator

    calc = U2LimitCalculator(lumi_fb=3000.0)  # HL-LHC
    U2_limit = calc.find_U2_limit(mass_GeV=1.0, flavour='muon')
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, Optional, Tuple
from scipy.optimize import brentq
from scipy.interpolate import interp1d
import sys

# Add models to path
sys.path.append(str(Path(__file__).parent.parent / 'models'))
sys.path.append(str(Path(__file__).parent.parent / 'geometry'))

from hnl_model_hnlcalc import HNLModelHNLCalc
from per_parent_efficiency import PerParentEfficiency


class U2LimitCalculator:
    """
    Calculate |U|² exclusion limits from signal yield.

    For a given (mass, flavour), scans |U|² to find the value where
    the expected signal yield N_sig reaches the exclusion threshold.
    """

    def __init__(self,
                 lumi_fb: float = 3000.0,
                 n_limit: float = 3.0,
                 cl: float = 0.95):
        """
        Initialize limit calculator.

        Parameters:
        -----------
        lumi_fb : float
            Integrated luminosity in fb⁻¹ (default: 3 ab⁻¹ = 3000 fb⁻¹)
        n_limit : float
            Number of events for exclusion (default: 3 for background-free 95% CL)
        cl : float
            Confidence level (default: 0.95)
        """
        self.lumi_fb = lumi_fb
        self.n_limit = n_limit
        self.cl = cl

        # Initialize model
        self.model = HNLModelHNLCalc()

        print(f"U2LimitCalculator initialized:")
        print(f"  Luminosity: {lumi_fb:.1f} fb⁻¹")
        print(f"  Exclusion threshold: {n_limit:.1f} events @ {cl*100:.0f}% CL")

    def compute_Nsig(self,
                     mass_GeV: float,
                     flavour: str,
                     U2: float,
                     efficiency_map: Dict) -> float:
        """
        Compute expected signal yield N_sig(m, flavour, |U|²).

        N_sig = L × Σ_parents [σ_p × BR(p→ℓN, |U|²) × ε_geom(p, cτ(|U|²))]

        Parameters:
        -----------
        mass_GeV : float
            HNL mass in GeV
        flavour : str
            Lepton flavour (electron, muon, tau)
        U2 : float
            Mixing parameter |U_flavour|²
        efficiency_map : Dict
            Per-parent efficiency map: {parent_pdg: {ctau: eps_geom}}

        Returns:
        --------
        float
            Expected number of signal events
        """
        # 1. Compute cτ at this |U|²
        ctau0 = self.model.ctau0_m_for_U2_eq_1(mass_GeV, flavour)
        ctau = ctau0 / U2

        # 2. Sum contributions from all parents
        N_sig_total = 0.0

        for parent_pdg, eps_vs_ctau in efficiency_map.items():
            # Get σ and BR for this parent
            sigma_pb = self.model.sigma_parent(parent_pdg, mass_GeV, flavour)
            BR = self.model.BR_parent_to_HNL(parent_pdg, mass_GeV, flavour, U2)

            # Interpolate ε_geom at this cτ
            ctau_array = np.array(sorted(eps_vs_ctau.keys()))
            eps_array = np.array([eps_vs_ctau[c] for c in ctau_array])

            # Use log-log interpolation (efficiency often varies over many orders)
            log_interp = interp1d(np.log10(ctau_array), np.log10(eps_array + 1e-20),
                                 kind='linear', bounds_error=False, fill_value=-20)

            log_eps = log_interp(np.log10(ctau))
            eps_geom = 10**log_eps

            # Contribution from this parent
            N_parent = self.lumi_fb * sigma_pb * BR * eps_geom

            N_sig_total += N_parent

            # Debug print for first few parents
            if N_parent > 1e-10:
                parent_name = self.model.get_parent_name(parent_pdg)
                print(f"    {parent_name:10s} (PDG {parent_pdg:6d}): "
                      f"σ={sigma_pb:.2e} pb, BR={BR:.2e}, ε={eps_geom:.2e} "
                      f"→ N={N_parent:.2e}")

        return N_sig_total

    def find_U2_limit(self,
                     mass_GeV: float,
                     flavour: str,
                     efficiency_map: Dict,
                     U2_min: float = 1e-12,
                     U2_max: float = 1.0) -> Tuple[float, Dict]:
        """
        Find |U|² limit where N_sig = N_limit.

        Uses root-finding (Brent's method) to solve:
            N_sig(m, flavour, |U|²) = N_limit

        Parameters:
        -----------
        mass_GeV : float
            HNL mass in GeV
        flavour : str
            Lepton flavour
        efficiency_map : Dict
            Per-parent efficiency map
        U2_min : float
            Minimum |U|² to search (default: 1e-12)
        U2_max : float
            Maximum |U|² to search (default: 1.0)

        Returns:
        --------
        Tuple[float, Dict]
            (U2_limit, info_dict)
            info_dict contains:
                - N_sig_at_limit
                - ctau_at_limit
                - dominant_parents
        """
        print(f"\nFinding |U|² limit for m={mass_GeV} GeV, flavour={flavour}")

        # Define objective function: N_sig(U2) - N_limit
        def objective(U2):
            N = self.compute_Nsig(mass_GeV, flavour, U2, efficiency_map)
            return N - self.n_limit

        # Check bounds
        N_at_min = objective(U2_min) + self.n_limit
        N_at_max = objective(U2_max) + self.n_limit

        print(f"  N_sig at |U|²={U2_min:.2e}: {N_at_min:.2e}")
        print(f"  N_sig at |U|²={U2_max:.2e}: {N_at_max:.2e}")

        if N_at_min < self.n_limit and N_at_max < self.n_limit:
            # No exclusion possible in this range
            print(f"  WARNING: No exclusion - N_sig < {self.n_limit} everywhere")
            return U2_max, {'N_sig_at_limit': N_at_max, 'ctau_at_limit': np.nan}

        if N_at_min > self.n_limit:
            # Even at minimum |U|², we have too many events
            print(f"  WARNING: Excluded even at |U|²={U2_min}")
            return U2_min, {'N_sig_at_limit': N_at_min, 'ctau_at_limit': np.nan}

        # Find root
        try:
            U2_limit = brentq(objective, U2_min, U2_max, xtol=1e-14)

            # Compute info at limit
            ctau0 = self.model.ctau0_m_for_U2_eq_1(mass_GeV, flavour)
            ctau_at_limit = ctau0 / U2_limit

            info = {
                'N_sig_at_limit': self.n_limit,
                'ctau_at_limit': ctau_at_limit,
                'U2_limit': U2_limit
            }

            print(f"  ✓ Found limit: |U|² = {U2_limit:.3e}")
            print(f"    cτ = {ctau_at_limit:.3e} m")

            return U2_limit, info

        except ValueError as e:
            print(f"  ERROR: Root finding failed: {e}")
            return np.nan, {'error': str(e)}

    def scan_mass_grid(self,
                      mass_grid: np.ndarray,
                      flavour: str,
                      efficiency_map_dir: Path) -> pd.DataFrame:
        """
        Scan full mass grid to produce |U|² limit curve.

        Parameters:
        -----------
        mass_grid : np.ndarray
            Array of HNL masses in GeV
        flavour : str
            Lepton flavour
        efficiency_map_dir : Path
            Directory containing efficiency map pickles

        Returns:
        --------
        pd.DataFrame
            Columns: mass_GeV, flavour, U2_limit, ctau_at_limit, N_sig_at_limit
        """
        results = []

        for mass in mass_grid:
            print(f"\n{'='*70}")
            print(f"Mass point: {mass} GeV ({flavour})")
            print(f"{'='*70}")

            # Find efficiency map file
            # Expected: HNL_mass_1.0_muon_Meson_efficiency_map.pkl
            pkl_files = list(efficiency_map_dir.glob(f"HNL_mass_{mass}_"
                                                     f"{flavour}_*_efficiency_map.pkl"))

            if not pkl_files:
                print(f"WARNING: No efficiency map found for {mass} GeV {flavour}")
                results.append({
                    'mass_GeV': mass,
                    'flavour': flavour,
                    'U2_limit': np.nan,
                    'ctau_at_limit': np.nan,
                    'N_sig_at_limit': np.nan
                })
                continue

            # Load efficiency map
            efficiency_pkl = pkl_files[0]
            print(f"Loading efficiency map: {efficiency_pkl.name}")

            from per_parent_efficiency import PerParentEfficiency
            calc = PerParentEfficiency()
            efficiency_map = calc.load_efficiency_map(efficiency_pkl)

            # Find limit
            U2_limit, info = self.find_U2_limit(mass, flavour, efficiency_map)

            # Store result
            results.append({
                'mass_GeV': mass,
                'flavour': flavour,
                'U2_limit': U2_limit,
                'ctau_at_limit': info.get('ctau_at_limit', np.nan),
                'N_sig_at_limit': info.get('N_sig_at_limit', np.nan)
            })

        return pd.DataFrame(results)


def main():
    """Example usage: compute limits for a mass grid."""
    import argparse

    parser = argparse.ArgumentParser(description="Compute |U|² limits")
    parser.add_argument('--flavour', type=str, default='muon',
                       choices=['electron', 'muon', 'tau'],
                       help='Lepton flavour')
    parser.add_argument('--lumi', type=float, default=3000.0,
                       help='Luminosity in fb⁻¹')
    parser.add_argument('--efficiency-dir', type=Path, required=True,
                       help='Directory with efficiency maps')
    parser.add_argument('--output', type=Path, default='U2_limits.csv',
                       help='Output CSV file')

    args = parser.parse_args()

    # Mass grid (example - adjust based on available data)
    mass_grid = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0,
                         10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0])

    # Initialize calculator
    calc = U2LimitCalculator(lumi_fb=args.lumi)

    # Scan mass grid
    results = calc.scan_mass_grid(mass_grid, args.flavour, args.efficiency_dir)

    # Save results
    results.to_csv(args.output, index=False)
    print(f"\n{'='*70}")
    print(f"Results saved to {args.output}")
    print(f"{'='*70}")

    # Print summary
    print("\nSummary:")
    print(results.to_string(index=False))


if __name__ == "__main__":
    main()
