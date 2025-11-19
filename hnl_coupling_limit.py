"""
Convert BR vs lifetime exclusions to coupling^2 vs mass exclusion (the "money plot" for HNLs)

For Heavy Neutral Leptons (HNLs), the production BR and lifetime are coupled:
- BR(W → ℓ N) ∝ |U_ℓ|^2 × f(m_N/m_W)  [production]
- τ_N ∝ 1 / |U_ℓ|^2                     [decay]

This script:
1. Takes BR vs lifetime exclusion curves from decayProbPerEvent.py analysis
2. Uses HNL physics to convert to |U_ℓ|^2 vs mass exclusion
3. Creates the standard HNL coupling plot
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import argparse
import os

# Physical constants
HBAR_GEV_S = 6.582119569e-25  # ℏ in GeV·s
G_F = 1.1663787e-5  # Fermi constant in GeV^-2
M_W = 80.377  # W boson mass in GeV
M_MU = 0.10566  # Muon mass in GeV
M_TAU = 1.77686  # Tau mass in GeV
SPEED_OF_LIGHT = 299792458.0  # m/s

def phase_space_factor(m_N, m_W=M_W, m_lepton=M_MU):
    """
    Phase space factor for W → ℓ N decay

    For W → ℓ N, the phase space suppression is:
    f(m_N, m_W, m_ℓ) = [(1 - x_N^2)^2 - x_ℓ^2(1 + x_N^2)]
    where x_N = m_N/m_W and x_ℓ = m_ℓ/m_W

    This is proportional to the partial width ratio.
    """
    x_N = m_N / m_W
    x_l = m_lepton / m_W

    # Kinematic threshold
    if m_N + m_lepton >= m_W:
        return 0.0

    # Phase space factor (proportional to decay width)
    f = (1 - x_N**2)**2 - x_l**2 * (1 + x_N**2)

    return max(f, 0.0)

def br_w_to_hnl(coupling_sq, m_N, lepton='mu'):
    """
    Calculate BR(W → ℓ N) as a function of |U_ℓ|^2 and m_N

    The branching ratio is:
    BR(W → ℓ N) = |U_ℓ|^2 × f(m_N) / Γ_W_total

    For simplicity, we normalize to the W → ℓ ν branching ratio:
    BR(W → ℓ N) / BR(W → ℓ ν) ≈ |U_ℓ|^2 × f(m_N)

    Since BR(W → ℓ ν) ≈ 1/9 for each lepton family, we get:
    BR(W → ℓ N) ≈ (1/9) × |U_ℓ|^2 × f(m_N)

    Args:
        coupling_sq: |U_ℓ|^2
        m_N: HNL mass in GeV
        lepton: 'mu' or 'tau'

    Returns:
        Branching ratio BR(W → ℓ N)
    """
    m_lepton = M_MU if lepton == 'mu' else M_TAU

    # Phase space factor
    f = phase_space_factor(m_N, M_W, m_lepton)

    # Approximate normalization (relative to W → ℓ ν)
    # BR(W → ℓ ν) ≈ 1/9 for each generation
    br_w_lnu = 1.0 / 9.0

    # BR(W → ℓ N) = |U_ℓ|^2 × f × BR(W → ℓ ν)
    br = coupling_sq * f * br_w_lnu

    return br

def hnl_lifetime(coupling_sq, m_N):
    """
    Calculate HNL lifetime as a function of |U_ℓ|^2 and m_N

    The HNL decay width is:
    Γ_N ≈ C × |U_ℓ|^2 × G_F^2 × m_N^5

    Calibrated using CMS 2024 data: for m=10 GeV, |U|^2 = 5×10^-7 gives cτ = 17 mm

    This gives C ≈ 1.7×10^-3 (dimensionless constant accounting for phase space
    and decay channels)

    The lifetime is τ_N = ℏ / Γ_N

    Args:
        coupling_sq: |U_ℓ|^2
        m_N: HNL mass in GeV

    Returns:
        Lifetime in seconds
    """
    # Calibration: CMS 2024 result gives C ≈ 1.7×10^-3
    # τ = cτ/c = 0.017 m / 3×10^8 m/s = 5.67×10^-11 s
    # Γ = ℏ/τ = 6.582×10^-25 / 5.67×10^-11 = 1.16×10^-14 GeV
    # Γ = C × 5×10^-7 × (1.166×10^-5)^2 × 10^5 = C × 6.8×10^-12
    # C = 1.16×10^-14 / 6.8×10^-12 ≈ 1.7×10^-3

    C = 1.7e-3

    # Γ_N in GeV
    gamma_N = coupling_sq * C * G_F**2 * m_N**5

    # Lifetime τ = ℏ / Γ
    if gamma_N > 0:
        tau = HBAR_GEV_S / gamma_N  # in seconds
    else:
        tau = np.inf

    return tau

def load_br_vs_lifetime_data(mass_gev, scenario='mu'):
    """
    Load BR vs lifetime exclusion data for a given mass

    This function looks for existing analysis output from decayProbPerEvent.py

    Args:
        mass_gev: HNL mass in GeV
        scenario: 'mu' or 'tau'

    Returns:
        (lifetimes_array, br_limits_array) or None if data not available
    """
    # Check for the exclusion data CSV file
    if scenario == 'mu':
        exclusion_file = f"output/csv/hnlLL_m{mass_gev}GeVLLP_exclusion_data.csv"
    else:
        exclusion_file = f"output/csv/hnlTauLL_m{mass_gev}GeVLLP_exclusion_data.csv"

    if not os.path.exists(exclusion_file):
        print(f"Warning: Exclusion data not found: {exclusion_file}")
        print(f"Run: python decayProbPerEvent.py output/csv/{'hnlLL' if scenario == 'mu' else 'hnlTauLL'}_m{mass_gev}GeVLLP.csv")
        return None

    # Load the exclusion data
    df = pd.read_csv(exclusion_file)

    lifetimes = df['lifetime_s'].values
    br_limits = df['BR_limit'].values

    print(f"Loaded exclusion data for m={mass_gev} GeV ({scenario}): {len(lifetimes)} points")

    return lifetimes, br_limits

def compute_coupling_limit(mass_gev, br_limits, lifetimes, lepton='mu'):
    """
    Convert BR vs lifetime exclusion to coupling^2 limit

    For each coupling value, we calculate:
    - BR(|U|^2, m_N) - the expected BR for that coupling
    - τ(|U|^2, m_N) - the expected lifetime for that coupling

    The coupling is excluded if BR(|U|^2) > BR_limit(τ(|U|^2))

    Args:
        mass_gev: HNL mass in GeV
        br_limits: Array of BR limits
        lifetimes: Array of lifetimes (seconds)
        lepton: 'mu' or 'tau'

    Returns:
        |U_ℓ|^2 limit for this mass
    """
    # Create interpolation function for BR_limit(τ)
    # Use log-log interpolation since both axes are log-scale
    log_tau = np.log10(lifetimes)
    log_br = np.log10(br_limits)

    # Remove any infinite or nan values
    valid = np.isfinite(log_tau) & np.isfinite(log_br)
    if not np.any(valid):
        return np.nan

    br_limit_interp = interp1d(log_tau[valid], log_br[valid],
                               kind='linear', bounds_error=False,
                               fill_value=(log_br[valid][0], log_br[valid][-1]))

    # Scan over coupling values
    log_coupling_sq_range = np.linspace(-12, 0, 1000)  # |U|^2 from 10^-12 to 1
    coupling_sq_range = 10**log_coupling_sq_range

    excluded = []

    for coupling_sq in coupling_sq_range:
        # Calculate expected BR and lifetime for this coupling
        br_expected = br_w_to_hnl(coupling_sq, mass_gev, lepton)
        tau_expected = hnl_lifetime(coupling_sq, mass_gev)

        # Get BR limit at this lifetime
        log_tau_exp = np.log10(tau_expected)
        log_br_lim = br_limit_interp(log_tau_exp)
        br_lim = 10**log_br_lim

        # Check if excluded
        is_excluded = (br_expected > br_lim)
        excluded.append(is_excluded)

    # Find the boundary (smallest coupling that's excluded)
    excluded = np.array(excluded)
    if np.any(excluded):
        # Find the transition point
        excluded_indices = np.where(excluded)[0]
        if len(excluded_indices) > 0:
            boundary_idx = excluded_indices[0]
            coupling_sq_limit = coupling_sq_range[boundary_idx]
            return coupling_sq_limit

    # No exclusion found
    return np.nan

def create_coupling_mass_plot(masses, coupling_limits_mu=None, coupling_limits_tau=None,
                             save_path='output/images/hnl_coupling_vs_mass.png'):
    """
    Create the "money plot": |U_ℓ|^2 vs mass with experimental comparisons

    Args:
        masses: Array of HNL masses in GeV
        coupling_limits_mu: Array of |U_μ|^2 limits (or None)
        coupling_limits_tau: Array of |U_τ|^2 limits (or None)
        save_path: Output file path
    """
    fig, ax = plt.subplots(figsize=(10, 8))

    # Plot the milliQan sensitivity
    if coupling_limits_mu is not None:
        valid = np.isfinite(coupling_limits_mu) & (coupling_limits_mu > 0)
        if np.any(valid):
            ax.loglog(masses[valid], coupling_limits_mu[valid],
                     'b-', linewidth=3, label='milliQan (muon-coupled)', marker='o')

    if coupling_limits_tau is not None:
        valid = np.isfinite(coupling_limits_tau) & (coupling_limits_tau > 0)
        if np.any(valid):
            ax.loglog(masses[valid], coupling_limits_tau[valid],
                     'r-', linewidth=3, label='milliQan (tau-coupled)', marker='s')

    # Load and plot experimental limits (if available)
    # Note: The external files (MATHUSLA, CODEX, ANUBIS) are in cτ vs BR format
    # They would need to be converted to coupling vs mass format as well
    # For now, we'll just plot our sensitivity

    ax.set_xlabel('HNL Mass (GeV)', fontsize=14)
    ax.set_ylabel(r'$|U_\ell|^2$ (Mixing Parameter)', fontsize=14)
    ax.set_title('Heavy Neutral Lepton Sensitivity at HL-LHC', fontsize=16)
    ax.grid(True, which='both', alpha=0.3)
    ax.legend(fontsize=12)

    # Set reasonable axis limits
    ax.set_xlim(10, 100)
    ax.set_ylim(1e-10, 1e-2)

    plt.tight_layout()
    plt.savefig(save_path, dpi=150)
    print(f"Saved coupling vs mass plot: {save_path}")

    return fig, ax

def main():
    """
    Main function to create HNL coupling limit plot

    This requires that decayProbPerEvent.py has been run for all mass points first.
    """
    parser = argparse.ArgumentParser(
        description='Create HNL coupling vs mass exclusion plot from BR vs lifetime data'
    )
    parser.add_argument('--scenario', type=str, default='mu', choices=['mu', 'tau'],
                       help='HNL coupling scenario: mu or tau')
    parser.add_argument('--test', action='store_true',
                       help='Run in test mode with synthetic data')
    args = parser.parse_args()

    # Define mass points to analyze
    masses = np.array([15, 23, 31, 39, 47, 55, 63, 71])  # GeV

    print("="*60)
    print(f"HNL COUPLING LIMIT ANALYSIS ({args.scenario}-coupled)")
    print("="*60)

    if args.test:
        print("\nRunning in TEST MODE with synthetic data...")
        print("This demonstrates the conversion logic without requiring full analysis.")
        print("")

        # Create synthetic BR vs lifetime data for demonstration
        lifetimes_test = np.logspace(-9, -5, 50)  # seconds

        coupling_limits = []

        for mass in masses:
            # Synthetic BR limits (decreasing with lifetime, as expected)
            # This is just for demonstration
            br_limits_test = 1e-3 * (lifetimes_test / 1e-7)**(-0.5)

            # Compute coupling limit
            coupling_limit = compute_coupling_limit(
                mass, br_limits_test, lifetimes_test, lepton=args.scenario
            )

            coupling_limits.append(coupling_limit)
            print(f"Mass {mass} GeV: |U_{args.scenario}|^2 limit = {coupling_limit:.2e}")

        coupling_limits = np.array(coupling_limits)

        # Create the plot
        if args.scenario == 'mu':
            create_coupling_mass_plot(masses, coupling_limits_mu=coupling_limits,
                                    save_path=f'output/images/hnl_coupling_vs_mass_{args.scenario}_test.png')
        else:
            create_coupling_mass_plot(masses, coupling_limits_tau=coupling_limits,
                                    save_path=f'output/images/hnl_coupling_vs_mass_{args.scenario}_test.png')

    else:
        print("\nLoading real exclusion data and computing coupling limits...")
        print("")

        coupling_limits = []
        masses_with_data = []

        for mass in masses:
            data = load_br_vs_lifetime_data(mass, scenario=args.scenario)

            if data is None:
                print(f"  Skipping mass {mass} GeV (no data)")
                continue

            lifetimes, br_limits = data

            # Compute coupling limit
            coupling_limit = compute_coupling_limit(
                mass, br_limits, lifetimes, lepton=args.scenario
            )

            coupling_limits.append(coupling_limit)
            masses_with_data.append(mass)

            print(f"  Mass {mass} GeV: |U_{args.scenario}|^2 limit = {coupling_limit:.2e}")

        print("")

        if len(masses_with_data) == 0:
            print("No data available. Please run decayProbPerEvent.py for the mass points first.")
            print("Example:")
            print(f"  python decayProbPerEvent.py output/csv/{'hnlLL' if args.scenario == 'mu' else 'hnlTauLL'}_m31GeVLLP.csv")
            return

        coupling_limits = np.array(coupling_limits)
        masses_with_data = np.array(masses_with_data)

        # Create the plot
        if args.scenario == 'mu':
            create_coupling_mass_plot(masses_with_data, coupling_limits_mu=coupling_limits,
                                    save_path=f'output/images/hnl_coupling_vs_mass_{args.scenario}.png')
        else:
            create_coupling_mass_plot(masses_with_data, coupling_limits_tau=coupling_limits,
                                    save_path=f'output/images/hnl_coupling_vs_mass_{args.scenario}.png')

if __name__ == "__main__":
    main()
