#!/usr/bin/env python3
"""
analyze_hnl_output.py

Quick analysis of HNL production output for validation.
Produces basic kinematic distributions and sanity checks.

Usage:
    python analyze_hnl_output.py <input.csv> [--plots]
"""

import sys
import pandas as pd
import numpy as np

def analyze_csv(filename, make_plots=False):
    """Analyze HNL production CSV output."""
    
    print(f"\n{'='*60}")
    print(f"Analyzing: {filename}")
    print('='*60)
    
    # Load data
    try:
        df = pd.read_csv(filename)
    except Exception as e:
        print(f"Error reading file: {e}")
        return
    
    print(f"\nTotal HNLs: {len(df)}")
    print(f"Unique events: {df['event'].nunique()}")
    print(f"HNLs per event: {len(df) / df['event'].nunique():.2f}")
    
    # Basic kinematics
    print(f"\n--- Kinematics ---")
    print(f"pT:   min={df['pt'].min():.3f}, max={df['pt'].max():.3f}, mean={df['pt'].mean():.3f} GeV")
    print(f"eta:  min={df['eta'].min():.3f}, max={df['eta'].max():.3f}, mean={df['eta'].mean():.3f}")
    print(f"E:    min={df['E'].min():.3f}, max={df['E'].max():.3f}, mean={df['E'].mean():.3f} GeV")
    print(f"mass: min={df['mass'].min():.4f}, max={df['mass'].max():.4f}, mean={df['mass'].mean():.4f} GeV")
    
    # Boost factor
    print(f"\n--- Boost Factor (γ = E/m) ---")
    print(f"γ:    min={df['boost_gamma'].min():.1f}, max={df['boost_gamma'].max():.1f}, mean={df['boost_gamma'].mean():.1f}")
    print(f"βγ:   mean={np.sqrt(df['boost_gamma'].mean()**2 - 1):.1f}")
    
    # Production vertex
    print(f"\n--- Production Vertex ---")
    r_prod = np.sqrt(df['prod_x_mm']**2 + df['prod_y_mm']**2)
    print(f"r_T:  min={r_prod.min():.3f}, max={r_prod.max():.3f}, mean={r_prod.mean():.3f} mm")
    print(f"z:    min={df['prod_z_mm'].min():.3f}, max={df['prod_z_mm'].max():.3f}, mean={df['prod_z_mm'].mean():.3f} mm")
    
    # Parent particles
    print(f"\n--- Parent Particles ---")
    parent_counts = df['parent_pdg'].value_counts()
    
    # PDG name lookup (partial)
    pdg_names = {
        321: 'K+', -321: 'K-',
        411: 'D+', -411: 'D-',
        421: 'D0', -421: 'D0bar',
        431: 'Ds+', -431: 'Ds-',
        511: 'B0', -511: 'B0bar',
        521: 'B+', -521: 'B-',
        531: 'Bs', -531: 'Bsbar',
        541: 'Bc+', -541: 'Bc-',
        24: 'W+', -24: 'W-',
        23: 'Z',
        15: 'tau-', -15: 'tau+',
        0: 'Unknown'
    }
    
    for pdg, count in parent_counts.head(10).items():
        name = pdg_names.get(pdg, f'PDG:{pdg}')
        frac = 100 * count / len(df)
        print(f"  {name:10s} ({pdg:6d}): {count:6d} ({frac:5.1f}%)")
    
    # HNL ID distribution (should be symmetric)
    print(f"\n--- HNL Charge ---")
    hnl_counts = df['hnl_id'].value_counts()
    for hnl_id, count in hnl_counts.items():
        print(f"  ID {hnl_id}: {count} ({100*count/len(df):.1f}%)")
    
    # Angular coverage for detectors
    print(f"\n--- Angular Coverage ---")
    eta_central = ((df['eta'] > -2.5) & (df['eta'] < 2.5)).sum()
    eta_forward = (np.abs(df['eta']) > 2.5).sum()
    print(f"  |η| < 2.5 (central): {eta_central} ({100*eta_central/len(df):.1f}%)")
    print(f"  |η| > 2.5 (forward): {eta_forward} ({100*eta_forward/len(df):.1f}%)")
    
    # Transverse momentum for MATHUSLA-like acceptance
    pt_high = (df['pt'] > 10).sum()
    print(f"  pT > 10 GeV:         {pt_high} ({100*pt_high/len(df):.1f}%)")
    
    if make_plots:
        try:
            import matplotlib.pyplot as plt
            
            fig, axes = plt.subplots(2, 3, figsize=(15, 10))
            
            # pT distribution
            axes[0,0].hist(df['pt'], bins=50, histtype='step', linewidth=2)
            axes[0,0].set_xlabel('pT [GeV]')
            axes[0,0].set_ylabel('Events')
            axes[0,0].set_title('Transverse Momentum')
            axes[0,0].set_yscale('log')
            
            # eta distribution
            axes[0,1].hist(df['eta'], bins=50, histtype='step', linewidth=2)
            axes[0,1].set_xlabel('η')
            axes[0,1].set_ylabel('Events')
            axes[0,1].set_title('Pseudorapidity')
            
            # Energy distribution
            axes[0,2].hist(df['E'], bins=50, histtype='step', linewidth=2)
            axes[0,2].set_xlabel('E [GeV]')
            axes[0,2].set_ylabel('Events')
            axes[0,2].set_title('Energy')
            axes[0,2].set_yscale('log')
            
            # Boost factor
            axes[1,0].hist(df['boost_gamma'], bins=50, histtype='step', linewidth=2)
            axes[1,0].set_xlabel('γ')
            axes[1,0].set_ylabel('Events')
            axes[1,0].set_title('Boost Factor')
            axes[1,0].set_yscale('log')
            
            # 2D: eta vs pT
            h = axes[1,1].hist2d(df['eta'], df['pt'], bins=50, cmap='viridis')
            axes[1,1].set_xlabel('η')
            axes[1,1].set_ylabel('pT [GeV]')
            axes[1,1].set_title('η vs pT')
            plt.colorbar(h[3], ax=axes[1,1])
            
            # Production vertex r vs z
            axes[1,2].hist2d(df['prod_z_mm'], r_prod, bins=50, cmap='viridis',
                            range=[[-10, 10], [0, 5]])
            axes[1,2].set_xlabel('z [mm]')
            axes[1,2].set_ylabel('r [mm]')
            axes[1,2].set_title('Production Vertex')
            
            plt.tight_layout()
            
            outfile = filename.replace('.csv', '_plots.png')
            plt.savefig(outfile, dpi=150)
            print(f"\nPlots saved to: {outfile}")
            
        except ImportError:
            print("\nNote: Install matplotlib for plots (pip install matplotlib)")
    
    print(f"\n{'='*60}\n")

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    
    make_plots = '--plots' in sys.argv
    for arg in sys.argv[1:]:
        if arg.endswith('.csv'):
            analyze_csv(arg, make_plots)
