import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import argparse
from pandas.errors import EmptyDataError

def append_tip_point_if_needed(df_all_sorted, df_dedup_valid):
    if len(df_dedup_valid) == 0:
        return df_dedup_valid

    df_all_sorted = df_all_sorted.sort_values("mass_GeV").reset_index(drop=True)
    valid_mask = df_all_sorted["eps2_min"].notna() & df_all_sorted["eps2_max"].notna()
    if not valid_mask.any():
        return df_dedup_valid

    valid_positions = np.flatnonzero(valid_mask.to_numpy())
    last_valid_pos = int(valid_positions[-1])
    if last_valid_pos >= len(df_all_sorted) - 1:
        return df_dedup_valid

    next_pos = last_valid_pos + 1
    if bool(valid_mask.iloc[next_pos]):
        return df_dedup_valid

    peak0 = float(df_all_sorted.loc[last_valid_pos, "peak_events"])
    peak1 = float(df_all_sorted.loc[next_pos, "peak_events"])
    if not (np.isfinite(peak0) and np.isfinite(peak1)):
        return df_dedup_valid

    threshold = 3.0
    if not (peak0 > threshold and peak1 < threshold):
        return df_dedup_valid

    m0 = float(df_all_sorted.loc[last_valid_pos, "mass_GeV"])
    m1 = float(df_all_sorted.loc[next_pos, "mass_GeV"])
    if not (m1 > m0 > 0):
        return df_dedup_valid

    last = df_dedup_valid.iloc[-1]
    width = float(last["eps2_max"]) / float(last["eps2_min"])
    if not np.isfinite(width):
        return df_dedup_valid
    if width <= 1.05:
        return df_dedup_valid

    tip_mass = m0 + (threshold - peak0) * (m1 - m0) / (peak1 - peak0)
    tip_mass = float(np.clip(tip_mass, m0, m1))

    tip_u2 = float(np.sqrt(float(last["eps2_min"]) * float(last["eps2_max"])))
    tip_row = {"mass_GeV": tip_mass, "eps2_min": tip_u2, "eps2_max": tip_u2, "peak_events": threshold}
    tip_df = pd.DataFrame([tip_row])

    if tip_mass / m0 <= 1.06:
        df_dedup_valid = df_dedup_valid.iloc[:-1]
    return pd.concat([df_dedup_valid, tip_df], ignore_index=True)

def main():
    parser = argparse.ArgumentParser(description="Plot exclusion islands for HNL/ALP pipelines")
    parser.add_argument("--particle", choices=["hnl", "alp"], default="hnl", help="Particle hypothesis")
    parser.add_argument("--input", default=None, help="Override input CSV path")
    parser.add_argument("--output", default=None, help="Override output image path")
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parent.parent
    if args.particle == "hnl":
        in_csv = Path(args.input) if args.input else (repo_root / "output/csv/analysis/HNL_U2_limits_summary.csv")
    else:
        in_csv = Path(args.input) if args.input else (repo_root / "output/csv/analysis/ALP_fa_limits_summary.csv")

    try:
        df = pd.read_csv(in_csv)
    except EmptyDataError:
        if args.particle == "hnl":
            df = pd.DataFrame(columns=["mass_GeV", "flavour", "eps2_min", "eps2_max", "peak_events"])
        else:
            df = pd.DataFrame(columns=["mass_GeV", "benchmark", "fa_min", "fa_max", "peak_events"])

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    if args.particle == "hnl":
        panels = [("electron", "100"), ("muon", "010"), ("tau", "001")]
        for idx, (flavour, _) in enumerate(panels):
            ax = axes[idx]
            sel = df["flavour"] == flavour
            df_sel = df[sel].sort_values("mass_GeV")

            df_valid = df_sel[df_sel["eps2_min"].notna() & df_sel["eps2_max"].notna()]
            df_dedup = (
                df_valid.groupby("mass_GeV", as_index=False)
                .agg({"eps2_min": "min", "eps2_max": "max", "peak_events": "max"})
                .sort_values("mass_GeV")
                .reset_index(drop=True)
            )
            df_dedup = df_dedup[(df_dedup["eps2_min"] > 0) & (df_dedup["eps2_max"] > 0)]
            df_dedup = df_dedup[df_dedup["eps2_min"] <= df_dedup["eps2_max"]]
            df_dedup = append_tip_point_if_needed(df_sel, df_dedup)

            if len(df_dedup) > 0:
                mass = df_dedup["mass_GeV"].values
                u2_min = df_dedup["eps2_min"].values
                u2_max = df_dedup["eps2_max"].values

                ax.fill_between(mass, u2_min, u2_max, alpha=0.3, color="red")
                ax.plot(mass, u2_min, "r-", linewidth=2, marker="o", markersize=3)
                ax.plot(mass, u2_max, "b-", linewidth=2, marker="o", markersize=3)

                from matplotlib.patches import Patch
                from matplotlib.lines import Line2D

                legend_elements = [
                    Patch(facecolor="red", alpha=0.3, label="Excluded"),
                    Line2D([0], [0], color="r", linewidth=2, label="Too long-lived"),
                    Line2D([0], [0], color="b", linewidth=2, label="Too prompt"),
                ]
                ax.legend(handles=legend_elements, fontsize=10)

                ax.set_xlabel(r"$m_N$ [GeV]", fontsize=14)
                ax.set_ylabel(r"$|U_{\ell}|^2$", fontsize=14)
                ax.set_title(f"{flavour.capitalize()}-coupled HNL", fontsize=14)
                ax.set_yscale("log")
                ax.set_xscale("log")

                ax.set_xlim([0.5, 50] if flavour == "tau" else [0.2, 50])
                ax.grid(True, which="both", alpha=0.3)
            else:
                ax.text(0.5, 0.5, f"No valid limits for {flavour}", ha="center", va="center", transform=ax.transAxes)
                ax.set_xlabel(r"$m_N$ [GeV]", fontsize=14)
                ax.set_ylabel(r"$|U_{\ell}|^2$", fontsize=14)

        out = Path(args.output) if args.output else (repo_root / "output/images/HNL_moneyplot_island.png")
    else:
        benchmarks = ["BC9", "BC10", "BC11"]
        for idx, bm in enumerate(benchmarks):
            ax = axes[idx]
            if "benchmark" not in df.columns:
                ax.text(0.5, 0.5, "Missing 'benchmark' column", ha="center", va="center", transform=ax.transAxes)
                continue

            df_sel = df[df["benchmark"].astype(str).str.upper() == bm].sort_values("mass_GeV")
            df_valid = df_sel[df_sel["fa_min"].notna() & df_sel["fa_max"].notna()]
            df_dedup = (
                df_valid.groupby("mass_GeV", as_index=False)
                .agg({"fa_min": "min", "fa_max": "max", "peak_events": "max"})
                .sort_values("mass_GeV")
                .reset_index(drop=True)
            )
            df_dedup = df_dedup[(df_dedup["fa_min"] > 0) & (df_dedup["fa_max"] > 0)]
            df_dedup = df_dedup[df_dedup["fa_min"] <= df_dedup["fa_max"]]

            if len(df_dedup) > 0:
                mass = df_dedup["mass_GeV"].values
                inv_fa_prompt = 1.0 / df_dedup["fa_min"].values  # strong coupling (too prompt)
                inv_fa_long = 1.0 / df_dedup["fa_max"].values    # weak coupling (too long-lived)

                y_min = np.minimum(inv_fa_long, inv_fa_prompt)
                y_max = np.maximum(inv_fa_long, inv_fa_prompt)

                ax.fill_between(mass, y_min, y_max, alpha=0.3, color="red")
                ax.plot(mass, y_min, "r-", linewidth=2, marker="o", markersize=3)
                ax.plot(mass, y_max, "b-", linewidth=2, marker="o", markersize=3)

                from matplotlib.patches import Patch
                from matplotlib.lines import Line2D

                legend_elements = [
                    Patch(facecolor="red", alpha=0.3, label="Excluded"),
                    Line2D([0], [0], color="r", linewidth=2, label="Too long-lived"),
                    Line2D([0], [0], color="b", linewidth=2, label="Too prompt"),
                ]
                ax.legend(handles=legend_elements, fontsize=10)

                ax.set_xlabel(r"$m_a$ [GeV]", fontsize=14)
                ax.set_ylabel(r"$1/f_a$ [GeV$^{-1}$]", fontsize=14)
                ax.set_title(f"ALP benchmark {bm}", fontsize=14)
                ax.set_yscale("log")
                ax.set_xscale("log")
                ax.grid(True, which="both", alpha=0.3)
            else:
                ax.text(0.5, 0.5, f"No valid limits for {bm}", ha="center", va="center", transform=ax.transAxes)
                ax.set_xlabel(r"$m_a$ [GeV]", fontsize=14)
                ax.set_ylabel(r"$1/f_a$ [GeV$^{-1}$]", fontsize=14)

        out = Path(args.output) if args.output else (repo_root / "output/images/ALP_moneyplot_island.png")

    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches='tight')
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
