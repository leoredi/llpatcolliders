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

    # Remove duplicates by keeping best sensitivity for each mass
    # Deduplicate per-mass entries and enforce sensible ordering/positivity before plotting
    df_dedup = df_valid.groupby("mass_GeV", as_index=False).agg({
        "eps2_min": "min",
        "eps2_max": "max",
        "peak_events": "max"
    }).sort_values("mass_GeV").reset_index(drop=True)
    df_dedup = df_dedup[(df_dedup["eps2_min"] > 0) & (df_dedup["eps2_max"] > 0)]
    df_dedup = df_dedup[df_dedup["eps2_min"] <= df_dedup["eps2_max"]]

    if len(df_dedup) > 0:
        # Split at 5 GeV to handle meson→EW transition discontinuity
        df_low = df_dedup[df_dedup["mass_GeV"] < 5.0]
        df_high = df_dedup[df_dedup["mass_GeV"] >= 5.0]

        # Plot low-mass regime (meson production)
        if len(df_low) > 0:
            mass_low = df_low["mass_GeV"].values
            u2_min_low = df_low["eps2_min"].values
            u2_max_low = df_low["eps2_max"].values
            ax.fill_between(mass_low, u2_min_low, u2_max_low, alpha=0.3, color='red')
            ax.plot(mass_low, u2_min_low, 'r-', linewidth=2, marker='o', markersize=3)
            ax.plot(mass_low, u2_max_low, 'b-', linewidth=2, marker='o', markersize=3)

        # Plot high-mass regime (EW production)
        if len(df_high) > 0:
            mass_high = df_high["mass_GeV"].values
            u2_min_high = df_high["eps2_min"].values
            u2_max_high = df_high["eps2_max"].values
            ax.fill_between(mass_high, u2_min_high, u2_max_high, alpha=0.3, color='red')
            ax.plot(mass_high, u2_min_high, 'r-', linewidth=2, marker='s', markersize=4)
            ax.plot(mass_high, u2_max_high, 'b-', linewidth=2, marker='s', markersize=4)

        # Legend
        from matplotlib.patches import Patch
        from matplotlib.lines import Line2D
        legend_elements = [
            Patch(facecolor='red', alpha=0.3, label='Excluded'),
            Line2D([0], [0], color='r', linewidth=2, label='Too long-lived'),
            Line2D([0], [0], color='b', linewidth=2, label='Too prompt')
        ]
        ax.legend(handles=legend_elements, fontsize=10)

        ax.set_xlabel(r"$m_N$ [GeV]", fontsize=14)
        ax.set_ylabel(r"$|U_{\ell}|^2$", fontsize=14)
        ax.set_title(f"{flavour.capitalize()}-coupled HNL", fontsize=14)
        ax.set_yscale('log')
        ax.set_xscale('log')

        if flavour == "tau":
            ax.set_xlim([0.5, 50])
        else:
            ax.set_xlim([0.2, 50])

        ax.grid(True, which="both", alpha=0.3)
    else:
        ax.text(0.5, 0.5, f"No valid limits for {flavour}",
                ha='center', va='center', transform=ax.transAxes)
        ax.set_xlabel(r"$m_N$ [GeV]", fontsize=14)
        ax.set_ylabel(r"$|U_{\ell}|^2$", fontsize=14)

plt.tight_layout()
plt.savefig("../output/images/HNL_moneyplot_island.png", dpi=150, bbox_inches='tight')
print("Saved: ../output/images/HNL_moneyplot_island.png")
