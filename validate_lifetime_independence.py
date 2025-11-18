#!/usr/bin/env python3
"""
Lifetime Independence Validation Study

This script validates whether HNL kinematics and detector acceptance
depend on particle lifetime by comparing distributions from simulations
at three different tau values (1ns, 30ns, 1us) at m = 39 GeV.

Usage:
    python validate_lifetime_independence.py
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
from scipy import stats
import os

# Configuration
DATA_DIR = "output/validation"
OUTPUT_DIR = "output/validation"
MASS_POINT = "39GeV"

# File paths
FILES = {
    '1ns': f'{DATA_DIR}/hnl{MASS_POINT}_tau1nsLLP.csv',
    '30ns': f'{DATA_DIR}/hnl{MASS_POINT}_tau30nsLLP.csv',
    '1us': f'{DATA_DIR}/hnl{MASS_POINT}_tau1usLLP.csv'
}

# Detector geometry (simplified acceptance cuts)
# Based on detector at z=22m with radius ~1.54m
ETA_MIN = -5.0  # Very forward acceptance
ETA_MAX = 5.0
PT_MIN = 5.0    # Minimum pT in GeV

def load_data(filepath):
    """Load CSV data and clean column names"""
    df = pd.read_csv(filepath, sep=',\s*', engine='python')
    df.columns = df.columns.str.strip()
    return df

def apply_acceptance_cuts(df):
    """Apply detector acceptance cuts"""
    acceptance_mask = (
        (df['eta'] >= ETA_MIN) &
        (df['eta'] <= ETA_MAX) &
        (df['pt'] >= PT_MIN)
    )
    return acceptance_mask

def calculate_statistics(data_dict):
    """Calculate statistical comparisons between datasets"""

    lifetimes = list(data_dict.keys())
    variables = ['pt', 'eta', 'phi', 'momentum']

    results = {
        'ks_tests': {},
        'chi2_tests': {},
        'acceptance_rates': {},
        'event_counts': {},
        'particle_counts': {}
    }

    # Calculate acceptance rates
    for lifetime, df in data_dict.items():
        acceptance = apply_acceptance_cuts(df)
        results['acceptance_rates'][lifetime] = acceptance.sum() / len(df)
        results['event_counts'][lifetime] = df['event'].nunique()
        results['particle_counts'][lifetime] = len(df)

    # KS tests comparing all pairs
    for var in variables:
        results['ks_tests'][var] = {}
        for i in range(len(lifetimes)):
            for j in range(i+1, len(lifetimes)):
                tau1, tau2 = lifetimes[i], lifetimes[j]
                data1 = data_dict[tau1][var].values
                data2 = data_dict[tau2][var].values

                ks_stat, ks_pval = stats.ks_2samp(data1, data2)
                results['ks_tests'][var][f'{tau1}_vs_{tau2}'] = {
                    'statistic': ks_stat,
                    'pvalue': ks_pval
                }

    return results

def plot_kinematic_comparisons(data_dict, output_dir):
    """Generate comparison plots for kinematic variables"""

    variables = {
        'pt': {'label': r'$p_T$ [GeV]', 'bins': 50, 'range': (0, 150)},
        'eta': {'label': r'$\eta$', 'bins': 50, 'range': (-6, 6)},
        'phi': {'label': r'$\phi$ [rad]', 'bins': 50, 'range': (-np.pi, np.pi)},
        'momentum': {'label': r'$p$ [GeV]', 'bins': 50, 'range': (0, 2000)}
    }

    colors = {'1ns': '#1f77b4', '30ns': '#ff7f0e', '1us': '#2ca02c'}

    for var, config in variables.items():
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))

        # Plot 1: Overlaid histograms
        for lifetime, df in data_dict.items():
            ax1.hist(df[var], bins=config['bins'], range=config['range'],
                    alpha=0.5, label=f'τ = {lifetime} (N={len(df)})',
                    color=colors[lifetime], density=True)

        ax1.set_xlabel(config['label'], fontsize=12)
        ax1.set_ylabel('Normalized Counts', fontsize=12)
        ax1.set_title(f'{config["label"]} Distribution Comparison (m = {MASS_POINT})',
                     fontsize=14, fontweight='bold')
        ax1.legend(fontsize=10)
        ax1.grid(True, alpha=0.3)

        # Plot 2: Cumulative distributions for better visual comparison
        for lifetime, df in data_dict.items():
            sorted_data = np.sort(df[var])
            cumulative = np.arange(1, len(sorted_data) + 1) / len(sorted_data)
            ax2.plot(sorted_data, cumulative, label=f'τ = {lifetime}',
                    color=colors[lifetime], linewidth=2, alpha=0.8)

        ax2.set_xlabel(config['label'], fontsize=12)
        ax2.set_ylabel('Cumulative Probability', fontsize=12)
        ax2.set_title(f'{config["label"]} Cumulative Distribution',
                     fontsize=14, fontweight='bold')
        ax2.legend(fontsize=10)
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        output_file = os.path.join(output_dir, f'{var}_comparison.png')
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Saved: {output_file}")
        plt.close()

def plot_acceptance_summary(results, output_dir):
    """Generate bar chart of acceptance rates"""

    fig, ax = plt.subplots(figsize=(10, 6))

    lifetimes = list(results['acceptance_rates'].keys())
    rates = [results['acceptance_rates'][tau] * 100 for tau in lifetimes]

    colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
    bars = ax.bar(lifetimes, rates, color=colors, alpha=0.7, edgecolor='black')

    # Add value labels on bars
    for bar, rate in zip(bars, rates):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'{rate:.2f}%',
                ha='center', va='bottom', fontsize=12, fontweight='bold')

    ax.set_ylabel('Acceptance Rate (%)', fontsize=12)
    ax.set_xlabel('Particle Lifetime', fontsize=12)
    ax.set_title(f'Detector Acceptance vs Lifetime (m = {MASS_POINT})',
                fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_ylim(0, max(rates) * 1.2)

    plt.tight_layout()
    output_file = os.path.join(output_dir, 'acceptance_summary.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()

def plot_correlation_analysis(data_dict, output_dir):
    """Generate correlation plots between kinematic variables"""

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    lifetimes = list(data_dict.keys())
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c']

    for idx, (lifetime, color) in enumerate(zip(lifetimes, colors)):
        df = data_dict[lifetime]

        # Plot eta vs phi
        axes[idx].scatter(df['phi'], df['eta'], alpha=0.3, s=10, color=color)
        axes[idx].set_xlabel(r'$\phi$ [rad]', fontsize=12)
        axes[idx].set_ylabel(r'$\eta$', fontsize=12)
        axes[idx].set_title(f'τ = {lifetime}', fontsize=14, fontweight='bold')
        axes[idx].grid(True, alpha=0.3)
        axes[idx].set_xlim(-np.pi, np.pi)
        axes[idx].set_ylim(-6, 6)

    plt.suptitle(r'Angular Distribution ($\eta$ vs $\phi$)',
                fontsize=16, fontweight='bold', y=1.02)
    plt.tight_layout()
    output_file = os.path.join(output_dir, 'correlation_analysis.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()

def generate_report(results, output_dir):
    """Generate markdown summary report"""

    report_path = os.path.join(output_dir, 'LIFETIME_VALIDATION_RESULTS.md')

    with open(report_path, 'w') as f:
        f.write("# Lifetime Independence Validation Results\n\n")
        f.write(f"**Mass Point:** {MASS_POINT}\n\n")
        f.write(f"**Validation Date:** {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        f.write("## Executive Summary\n\n")

        # Analyze KS test results
        all_ks_pvals = []
        for var, tests in results['ks_tests'].items():
            for comparison, test_result in tests.items():
                all_ks_pvals.append(test_result['pvalue'])

        min_pval = min(all_ks_pvals)
        max_pval = max(all_ks_pvals)
        mean_pval = np.mean(all_ks_pvals)

        # Analyze acceptance rates
        accept_rates = list(results['acceptance_rates'].values())
        accept_variation = (max(accept_rates) - min(accept_rates)) / np.mean(accept_rates) * 100

        f.write(f"- **KS Test p-values:** min={min_pval:.4f}, max={max_pval:.4f}, mean={mean_pval:.4f}\n")
        f.write(f"- **Acceptance rate variation:** {accept_variation:.2f}%\n\n")

        # Determine if validation passes
        ks_pass = min_pval > 0.05
        accept_pass = accept_variation < 10.0

        if ks_pass and accept_pass:
            f.write("**VALIDATION STATUS: ✅ PASSED**\n\n")
            f.write("Kinematics are statistically independent of lifetime. "
                   "All KS test p-values > 0.05 and acceptance variation < 10%.\n\n")
        else:
            f.write("**VALIDATION STATUS: ⚠️ NEEDS REVIEW**\n\n")
            if not ks_pass:
                f.write(f"- KS test failed: minimum p-value = {min_pval:.4f} < 0.05\n")
            if not accept_pass:
                f.write(f"- Acceptance variation too large: {accept_variation:.2f}% >= 10%\n")
            f.write("\n")

        f.write("## Dataset Statistics\n\n")
        f.write("| Lifetime | Events | Particles | Acceptance Rate |\n")
        f.write("|----------|--------|-----------|----------------|\n")
        for lifetime in ['1ns', '30ns', '1us']:
            events = results['event_counts'][lifetime]
            particles = results['particle_counts'][lifetime]
            accept = results['acceptance_rates'][lifetime] * 100
            f.write(f"| {lifetime} | {events:,} | {particles:,} | {accept:.2f}% |\n")

        f.write("\n## Statistical Tests\n\n")
        f.write("### Kolmogorov-Smirnov Tests\n\n")
        f.write("Tests the null hypothesis that two samples come from the same distribution. "
               "p-value > 0.05 indicates distributions are statistically identical.\n\n")

        for var in ['pt', 'eta', 'phi', 'momentum']:
            f.write(f"#### {var.upper()}\n\n")
            f.write("| Comparison | KS Statistic | p-value | Result |\n")
            f.write("|------------|--------------|---------|--------|\n")

            for comparison, test_result in results['ks_tests'][var].items():
                stat = test_result['statistic']
                pval = test_result['pvalue']
                result = "✅ Same" if pval > 0.05 else "⚠️ Different"
                f.write(f"| {comparison} | {stat:.4f} | {pval:.4f} | {result} |\n")
            f.write("\n")

        f.write("## Interpretation\n\n")

        if ks_pass and accept_pass:
            f.write("### ✅ Kinematics are Lifetime-Independent\n\n")
            f.write("The validation study confirms that HNL kinematics at production "
                   "do not depend on the particle lifetime (tau0 parameter in PYTHIA). "
                   "This means:\n\n")
            f.write("1. **Option 1 is valid:** We can generate kinematics at one lifetime "
                   "and re-weight for different lifetimes in post-processing\n")
            f.write("2. **Option 2 remains valid:** Brute-force approach with separate "
                   "simulations per lifetime point still works\n")
            f.write("3. **Option 3 is unnecessary:** We don't need to modify PYTHIA code "
                   "for custom lifetime handling\n\n")
            f.write("### Recommendation\n\n")
            f.write("**Proceed with Option 1** (generate once, re-weight in analysis) as it provides:\n")
            f.write("- Maximum computational efficiency\n")
            f.write("- Statistical validity (confirmed by this study)\n")
            f.write("- Clean separation between event generation and physics analysis\n\n")
        else:
            f.write("### ⚠️ Lifetime Dependence Detected\n\n")
            f.write("The validation study suggests some statistical differences between "
                   "lifetime scenarios. This could indicate:\n\n")
            f.write("1. PYTHIA's decay-in-flight handling affects production kinematics\n")
            f.write("2. Statistical fluctuations (consider increasing event count)\n")
            f.write("3. Detector acceptance effects\n\n")
            f.write("### Recommendation\n\n")
            f.write("**Further investigation needed:**\n")
            f.write("- Increase statistics to 50k-100k events per lifetime\n")
            f.write("- Examine which kinematic variables show dependence\n")
            f.write("- Consider Option 2 (brute force) if dependence is real\n\n")

        f.write("## Plots Generated\n\n")
        f.write("1. `pt_comparison.png` - Transverse momentum distributions\n")
        f.write("2. `eta_comparison.png` - Pseudorapidity distributions\n")
        f.write("3. `phi_comparison.png` - Azimuthal angle distributions\n")
        f.write("4. `momentum_comparison.png` - Total momentum distributions\n")
        f.write("5. `acceptance_summary.png` - Detector acceptance rates\n")
        f.write("6. `correlation_analysis.png` - Angular correlations\n\n")

        f.write("## Next Steps\n\n")

        if ks_pass and accept_pass:
            f.write("1. Implement Option 1 workflow:\n")
            f.write("   - Generate large statistics (1M events) at single reference lifetime\n")
            f.write("   - Modify `decayProbPerEvent.py` to scan over lifetime range\n")
            f.write("   - Generate exclusion limits with lifetime re-weighting\n")
            f.write("2. Validate results against Option 2 at a few mass points\n")
            f.write("3. Document the workflow in repository README\n")
        else:
            f.write("1. Run extended validation with 50k-100k events per lifetime\n")
            f.write("2. Investigate source of kinematic dependence\n")
            f.write("3. If dependence persists, proceed with Option 2 (brute force)\n")

    print(f"\nSaved: {report_path}")

def main():
    """Main validation workflow"""

    print("="*70)
    print("  Lifetime Independence Validation Study")
    print("  Mass Point: m = " + MASS_POINT)
    print("="*70)
    print()

    # Load data
    print("Loading datasets...")
    data_dict = {}
    for lifetime, filepath in FILES.items():
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"Data file not found: {filepath}")
        data_dict[lifetime] = load_data(filepath)
        print(f"  {lifetime}: {len(data_dict[lifetime])} particles "
              f"from {data_dict[lifetime]['event'].nunique()} events")
    print()

    # Calculate statistics
    print("Calculating statistical tests...")
    results = calculate_statistics(data_dict)
    print("  Completed KS tests and acceptance calculations")
    print()

    # Generate plots
    print("Generating comparison plots...")
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    plot_kinematic_comparisons(data_dict, OUTPUT_DIR)
    plot_acceptance_summary(results, OUTPUT_DIR)
    plot_correlation_analysis(data_dict, OUTPUT_DIR)
    print()

    # Generate report
    print("Generating validation report...")
    generate_report(results, OUTPUT_DIR)
    print()

    # Print summary to console
    print("="*70)
    print("  VALIDATION SUMMARY")
    print("="*70)
    print()
    print("Acceptance Rates:")
    for lifetime in ['1ns', '30ns', '1us']:
        rate = results['acceptance_rates'][lifetime] * 100
        print(f"  τ = {lifetime:4s}: {rate:6.2f}%")
    print()

    print("KS Test Results (minimum p-values):")
    for var in ['pt', 'eta', 'phi', 'momentum']:
        min_pval = min([test['pvalue'] for test in results['ks_tests'][var].values()])
        status = "✅" if min_pval > 0.05 else "⚠️"
        print(f"  {var:10s}: p = {min_pval:.4f} {status}")
    print()

    print("All results saved to:", OUTPUT_DIR)
    print("="*70)

if __name__ == "__main__":
    main()
