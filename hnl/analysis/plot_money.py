"""E14-16 -- the GRENDEL HNL money plot + A2 metadata + machine-readable bundle.

Built ON `analysis.plot_exclusion` (the proper renderer: excluded-region fill,
open-edge markers, ylim to the scan ceiling 1e-1, segment-aware drawing, and the
FONLL theory ribbon via `band_csv`). On top of that base it layers the two extra
in-scope boundary bands:
  * lower edge: FONLL production band (orange, from `combine_band.py`)  AND
                Bc normalization band (teal, `bc_nuisance.py`);
  * upper edge: decay-model width band (blue, `decay_model_band.py`).
Plus the F17 provenance line, the Umu low-mass `N->l pi` callout, the A2
hypothesis metadata (`run_metadata.json`), and a bundle dir with the tables.

Production + HNL->SM decay physics only; reconstruction / detector / background /
statistics are the partner layer, idealized (background-free, N>=3).

  python -m analysis.plot_money --run central_newgrids_20260623 --out-dir tmp/runs/v1_bundle
"""
from __future__ import annotations

import argparse
import json
import shutil
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from analysis.plot_exclusion import _plot_single_panel, PLOT_U2_MIN, PLOT_U2_MAX  # noqa: E402

HNL_ROOT = Path(__file__).resolve().parent.parent
RUNS = HNL_ROOT / "tmp" / "runs"
LAB = {"Ue": r"$|U_e|^2$", "Umu": r"$|U_\mu|^2$", "Utau": r"$|U_\tau|^2$"}


def _provenance(have_fonll):
    fonll = ("FONLL production (orange)" if have_fonll
             else "FONLL production (band PENDING -- not yet on this figure)")
    return (
        f"In-band: {fonll}; Bc direct-norm ±40% (teal, LHCb arXiv:1910.13404; induced-tau Bc 0.06% negligible); "
        "HNL total-width/lifetime duality δ(m) (blue, cap 20% conservative vs ~10% post-QCD residual).  "
        "NOT banded (named limitations): production form factors, absolute visible-BR norm, kaon flux, FONLL αs.  "
        "Reconstruction / detector / background / statistics IDEALIZED at the partner handoff (background-free, N≥3).  "
        "Island closes at the production/lifetime pinch ~3.7 GeV (last sensitive grid mass 3.6 GeV, gap 0.28 dex; "
        "3.8 GeV grid point insensitive -> physical closure; the plotted tip is a linear log-edge extrapolation)."
    )


def _fonll_dense(band_df, central_df, flavor):
    """Interpolate the (coarse) FONLL band onto the central fine grid in
    log10(U^2) so the orange ribbon is continuous, not sampled at band masses."""
    cm = np.sort(central_df.loc[central_df.flavor == flavor, "mass_GeV"].unique())
    b = band_df[band_df.flavor == flavor].sort_values("mass_GeV")
    out = {"flavor": flavor, "mass_GeV": cm}
    for col in ("u2_min", "u2_max"):
        for side in ("band_lo", "band_hi"):
            y = b[f"{col}_{side}"].to_numpy(float)
            x = b["mass_GeV"].to_numpy(float)
            ok = np.isfinite(y) & (y > 0)
            out[f"{col}_{side}"] = (
                10.0 ** np.interp(cm, x[ok], np.log10(y[ok]), left=np.nan, right=np.nan)
                if ok.sum() >= 2 else np.full(len(cm), np.nan))
    return pd.DataFrame(out)


def _overlay(ax, band_df, flavor, lo_col, hi_col, color, label):
    b = band_df[band_df.flavor == flavor].sort_values("mass_GeV")
    m = b["mass_GeV"].to_numpy(float)
    lo, hi = b[lo_col].to_numpy(float), b[hi_col].to_numpy(float)
    good = np.isfinite(lo) & np.isfinite(hi)
    if good.any():
        ax.fill_between(m, np.where(good, lo, np.nan), np.where(good, hi, np.nan),
                        alpha=0.55, color=color, linewidth=0, zorder=8, label=label)


def _metadata(run, l_int_fb, p_cut_mev, have_fonll):
    in_band = ["HNL total-width/lifetime duality (decay-model band)",
               "direct Bc normalization (B4)"]
    if have_fonll:
        in_band.insert(0, "FONLL heavy-flavor production (scale/PDF/m_Q)")
    lim = ["production form factors (B3)", "absolute visible-BR normalization (B5)",
           "kaon flux/transport (B7)", "FONLL alpha_s (sub-dominant)"]
    if not have_fonll:
        lim.insert(0, "FONLL scale/PDF/m_Q band -- PENDING, not on this figure")
    return {
        "result": "GRENDEL HNL sensitivity projection (single-flavor)",
        "hypothesis": {"mixing": "single-flavor", "flavors": ["Ue", "Umu", "Utau"],
                       "observable": "|U_alpha|^2", "nature": "Majorana",
                       "charge_conjugate_counting": "included (NDecayWidth x2/channel; production both charges)",
                       "production_mixing_equals_decay_mixing": True},
        "beam": {"collision": "pp", "sqrt_s_TeV": 14, "L_int_fb": l_int_fb},
        "limit": {"method": "background-free", "criterion": "N_signal >= 3",
                  "track_momentum_cut_MeV": p_cut_mev},
        "scope": {"in_band": in_band,
                  "idealized_partner_handoff": ["reconstruction", "detector response",
                                                "background", "statistics"]},
        "central_run": run,
        "band_sources": {"fonll_lower_edge": "hnl_band.csv (run_variation_band + combine_band)"
                                             if have_fonll else "PENDING",
                         "decay_model_upper_edge": "decay_model_band_combined.csv (width_band delta(m))",
                         "bc_lower_edge": "bc_nuisance_band.csv (SIGMA_BC_REL_UNCERT=0.40)"},
        "limitations_not_banded": lim,
        "known_features": ["Umu/Utau low-mass step where N->l pi closes (m_mu+m_pi=0.245 GeV) -- "
                           "real physics (REMAINING_WORK item 15)",
                           "island closes at the production/lifetime pinch ~3.68-3.75 GeV (Ue/Umu/Utau): "
                           "last sensitive grid mass 3.6 GeV (0.28 dex gap), 3.8 GeV insensitive -> physical "
                           "closure; the plotted tip is a linear log-edge extrapolation of the two boundaries"],
        "generated": str(date.today()),
        "reproduce": {"env": "HNL_TMP_DIR=tmp HNL_RUN_TAG=<run> HNL_P_CUT=0.600 (hnl conda env)",
                      "steps": ["python -m analysis.run_sensitivity",
                                "run_variation_band.py + python -m analysis.combine_band  # FONLL band -> hnl_band.csv",
                                "python -m analysis.decay_model_band",
                                "python -m analysis.bc_nuisance",
                                "python -m analysis.channel_breakdown",
                                "python -m analysis.plot_money"],
                      "width_band_delta_table": "python -m analysis.width_band  (ROOT-enabled env)"},
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default="central_newgrids_20260623")
    ap.add_argument("--central-csv", default=None)
    ap.add_argument("--fonll-band", default=str(RUNS / "hnl_band.csv"))
    ap.add_argument("--decay-band", default=str(RUNS / "decay_model_band_combined.csv"))
    ap.add_argument("--bc-band", default=str(RUNS / "bc_nuisance_band.csv"))
    ap.add_argument("--breakdown", default=str(RUNS / "channel_breakdown_u2min.csv"))
    ap.add_argument("--l-int-fb", type=float, default=3000.0)
    ap.add_argument("--p-cut-mev", type=int, default=600)
    ap.add_argument("--out-dir", default=str(RUNS / "v1_bundle"))
    a = ap.parse_args(argv)

    central_csv = a.central_csv or str(RUNS / a.run / "analysis_exact_600" / "hnl_sensitivity.csv")
    cen = pd.read_csv(central_csv)
    cen = cen[cen.has_sensitivity == True]                                  # noqa: E712
    have_fonll = Path(a.fonll_band).exists()
    fonll = pd.read_csv(a.fonll_band) if have_fonll else None
    dm = pd.read_csv(a.decay_band)
    bc = pd.read_csv(a.bc_band)

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.6), sharey=True)
    for ax, fl in zip(axes, ["Ue", "Umu", "Utau"]):
        band_df = _fonll_dense(fonll, cen, fl) if have_fonll else None
        _plot_single_panel(ax, cen, fl, is_leftmost=(fl == "Ue"), band_df=band_df,
                           close_island=True)
        _overlay(ax, dm, fl, "u2_max_dm_lo", "u2_max_dm_hi", "steelblue", "decay-model band (upper)")
        _overlay(ax, bc, fl, "u2_min_bc_lo", "u2_min_bc_hi", "teal", "Bc band (lower)")
        ax.set_title(f"HNL {LAB[fl]}", fontsize=13)
        if fl in ("Umu", "Utau"):
            ax.annotate(r"$N\to\ell\pi$ closes" + "\n→ weaker limit", xy=(0.25, 3e-6),
                        xytext=(0.33, 5e-5), fontsize=7, color="dimgray",
                        arrowprops=dict(arrowstyle="->", color="dimgray", lw=0.8))
        ax.legend(fontsize=7, loc="lower right")
    fig.suptitle("GRENDEL HNL sensitivity projection (single-flavor, Majorana, 14 TeV, "
                 f"{a.l_int_fb:.0f} fb$^{{-1}}$, P>{a.p_cut_mev} MeV)", y=1.0, fontsize=13)
    fig.text(0.5, 0.005, _provenance(have_fonll), ha="center", va="bottom",
             fontsize=6.2, style="italic", wrap=True)
    plt.tight_layout(rect=[0, 0.075, 1, 0.97])

    out = Path(a.out_dir); out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / "money_plot.png", dpi=130, bbox_inches="tight")
    fig.savefig(out / "money_plot.pdf", bbox_inches="tight")
    (out / "run_metadata.json").write_text(json.dumps(
        _metadata(a.run, a.l_int_fb, a.p_cut_mev, have_fonll), indent=2))
    for src, name in [(central_csv, "hnl_sensitivity_central.csv"),
                      (a.fonll_band, "hnl_band_fonll.csv"),
                      (a.decay_band, "decay_model_band.csv"),
                      (a.bc_band, "bc_nuisance_band.csv"),
                      (a.breakdown, "channel_breakdown_u2min.csv")]:
        if Path(src).exists():
            shutil.copy2(src, out / name)
    print(f"money plot + bundle -> {out}/  (FONLL band: {'YES' if have_fonll else 'PENDING'})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
