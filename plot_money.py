import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_csv("output/csv/analysis/HNL_U2_limits_summary.csv")

fig, axes = plt.subplots(1, 3, figsize=(18, 5))

for idx, flavour in enumerate(["electron", "muon", "tau"]):
    ax = axes[idx]
    sel = (df["flavour"] == flavour)
    df_sel = df[sel].sort_values("mass_GeV")
    
    # Remove NaN values
    df_valid = df_sel[df_sel["U2_limit"].notna()]
    
    if len(df_valid) > 0:
        ax.loglog(df_valid["mass_GeV"], df_valid["U2_limit"], marker="o", linewidth=2, markersize=8)
        ax.set_xlabel(r"$m_N$ [GeV]", fontsize=14)
        ax.set_ylabel(r"$|U_{\ell}|^2$ limit", fontsize=14)
        ax.set_title(f"{flavour.capitalize()}-coupled HNL", fontsize=14)
        ax.grid(True, which="both", alpha=0.3)
    else:
        ax.text(0.5, 0.5, f"No valid limits for {flavour}", 
                ha='center', va='center', transform=ax.transAxes)
        ax.set_xlabel(r"$m_N$ [GeV]", fontsize=14)
        ax.set_ylabel(r"$|U_{\ell}|^2$ limit", fontsize=14)

plt.tight_layout()
plt.savefig("output/images/HNL_moneyplot_all.png", dpi=150)
print("Saved: output/images/HNL_moneyplot_all.png")
