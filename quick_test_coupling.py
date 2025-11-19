"""
Quick test of coupling limit calculation using manually specified BR vs lifetime data

This demonstrates the conversion without waiting for full analysis
"""
import numpy as np
import matplotlib.pyplot as plt
from hnl_coupling_limit import (
    br_w_to_hnl, hnl_lifetime, compute_coupling_limit,
    create_coupling_mass_plot
)

# Test with a few example points
# These are rough estimates for demonstration
masses = np.array([31, 47, 63])  # GeV

# Manually specify some example BR limits vs lifetime
# (In reality, these come from decayProbPerEvent.py analysis)

# For m=31 GeV: assume BR limit ranges from 1e-4 to 1e-1 over lifetime range
lifetimes_31 = np.logspace(-8, -6, 10)  # 10 ns to 1 μs
br_limits_31 = 1e-3 * (lifetimes_31 / 1e-7)**(-0.3)  # Example scaling

# For m=47 GeV
lifetimes_47 = np.logspace(-8, -6, 10)
br_limits_47 = 2e-3 * (lifetimes_47 / 1e-7)**(-0.3)

# For m=63 GeV
lifetimes_63 = np.logspace(-8, -6, 10)
br_limits_63 = 5e-3 * (lifetimes_63 / 1e-7)**(-0.3)

# Compute coupling limits
coupling_limits = []

for mass, lifetimes, br_limits in [
    (31, lifetimes_31, br_limits_31),
    (47, lifetimes_47, br_limits_47),
    (63, lifetimes_63, br_limits_63)
]:
    coupling_limit = compute_coupling_limit(mass, br_limits, lifetimes, lepton='mu')
    coupling_limits.append(coupling_limit)
    print(f"Mass {mass} GeV: |U_μ|^2 limit = {coupling_limit:.2e}")

coupling_limits = np.array(coupling_limits)

# Create plot
create_coupling_mass_plot(masses, coupling_limits_mu=coupling_limits,
                         save_path='output/images/hnl_coupling_vs_mass_quick_test.png')

print("\nQuick test complete!")
print("This used synthetic BR vs lifetime data for demonstration.")
print("Run the full analysis to get real sensitivity limits.")
