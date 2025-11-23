import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_csv("../output/csv/analysis/HNL_U2_limits_summary.csv")

fig, axes = plt.subplots(1, 3, figsize=(18, 5))

for idx, flavour in enumerate(["electron", "muon", "tau"]):
    ax = axes[idx]
    sel = (df["flavour"] == flavour)
    df_sel = df[sel].sort_values("mass_GeV")

    # Remove NaN values (use eps2_min/eps2_max from PBC analysis)
    df_valid = df_sel[df_sel["eps2_min"].notna() & df_sel["eps2_max"].notna()]
    
    if len(df_valid) > 0:
        mass = df_valid["mass_GeV"].values
        u2_min = df_valid["eps2_min"].values
        u2_max = df_valid["eps2_max"].values
        
        # Plot the excluded island/band
        ax.fill_between(mass, u2_min, u2_max, alpha=0.3, color='red', label='Excluded')
        ax.plot(mass, u2_min, 'r-', linewidth=2, label='Too prompt')
        ax.plot(mass, u2_max, 'b-', linewidth=2, label='Too long-lived')
        
        ax.set_xlabel(r"$m_N$ [GeV]", fontsize=14)
        ax.set_ylabel(r"$|U_{\ell}|^2$", fontsize=14)
        ax.set_title(f"{flavour.capitalize()}-coupled HNL", fontsize=14)
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_xlim([0.2, 50])
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize=10)
    else:
        ax.text(0.5, 0.5, f"No valid limits for {flavour}", 
                ha='center', va='center', transform=ax.transAxes)
        ax.set_xlabel(r"$m_N$ [GeV]", fontsize=14)
        ax.set_ylabel(r"$|U_{\ell}|^2$", fontsize=14)

plt.tight_layout()
plt.savefig("../output/images/HNL_moneyplot_island.png", dpi=150)
print("Saved: ../output/images/HNL_moneyplot_island.png")
