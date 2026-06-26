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

from analysis.plot_exclusion import _plot_single_panel, _closure_vertex, _closure_arc, PLOT_U2_MIN, PLOT_U2_MAX  # noqa: E402

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
        "3.8 GeV grid point insensitive -> physical closure; the plotted tip is a rounded extrapolation to the interpolated pinch)."
    )


def _dex_densify(bm, b_cen, b_lo, b_hi, cm, c_edge, tail=None):
    """Densify one band edge onto the fine central grid.

    The ~18 band anchors are too coarse to plot as absolute edges -- straight
    log-segments between anchors detach from the wiggly dense central line. So we
    interpolate the band *width in dex* (relative to the band-run central) and
    re-apply it to the dense plotted central ``c_edge``, so the ribbon hugs every
    wiggle of the red curve. ``tail = (mass_arc, edge_arc)`` optionally appends the
    rounded closure-nose path so the ribbon collapses onto, and tapers shut with,
    the island's rounded edge instead of poking past it.
    """
    ok = (np.isfinite(b_cen) & (b_cen > 0) & np.isfinite(b_lo) & (b_lo > 0)
          & np.isfinite(b_hi) & (b_hi > 0))
    if ok.sum() >= 2:
        dex_lo = np.log10(b_cen[ok] / b_lo[ok])     # downward half-width [dex]
        dex_hi = np.log10(b_hi[ok] / b_cen[ok])     # upward half-width   [dex]
        lo = c_edge * 10.0 ** (-np.interp(cm, bm[ok], dex_lo))
        hi = c_edge * 10.0 ** (+np.interp(cm, bm[ok], dex_hi))
    else:
        lo = np.full(len(cm), np.nan)
        hi = np.full(len(cm), np.nan)
    m = cm
    if tail is not None and tail[0] is not None:
        tm, ty = tail
        m = np.append(cm, tm)
        lo = np.append(lo, ty)
        hi = np.append(hi, ty)
    return m, lo, hi


def _ribbon(ax, m, lo, hi, color, label, zorder):
    good = np.isfinite(lo) & np.isfinite(hi) & (lo > 0) & (hi > 0)
    if good.any():
        ax.fill_between(m, np.where(good, lo, np.nan), np.where(good, hi, np.nan),
                        alpha=0.5, color=color, linewidth=0, zorder=zorder, label=label)


def _closure_pinch(full_flavor_df):
    """Return the (m_star, u2_star) closure vertex for one flavor's curve, matching
    what `_plot_single_panel(close_island=True)` appends, or None if it does not
    close just past the grid."""
    sens = full_flavor_df[full_flavor_df["has_sensitivity"] == True]   # noqa: E712
    cm = sens["mass_GeV"].to_numpy(float)
    if len(cm) < 2:
        return None
    later = full_flavor_df[full_flavor_df["mass_GeV"] > cm[-1]]
    if later.empty or bool(later.iloc[0]["has_sensitivity"]):
        return None
    cv = _closure_vertex(cm, sens["u2_min"].to_numpy(float), sens["u2_max"].to_numpy(float))
    if cv is not None and cm[-1] < cv[0] <= float(later.iloc[0]["mass_GeV"]):
        return cv
    return None


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
                           "closure; the plotted tip is a rounded extrapolation to the interpolated pinch"],
        "generated": str(date.today()),
        "reproduce": {"env": f"HNL_TMP_DIR=tmp HNL_RUN_TAG=<run> HNL_P_CUT={p_cut_mev/1000:.3f} (hnl conda env)",
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
    # P>100 MeV is the default cut (exact_100); P>600 is legacy (pass
    # --central-csv .../analysis_exact_600/... --p-cut-mev 600 --decay-band ...
    # --bc-band ... to rebuild it).
    ap.add_argument("--run", default="central_newgrids_20260623")
    ap.add_argument("--central-csv", default=None)
    ap.add_argument("--fonll-band", default=str(RUNS / "hnl_band.csv"))
    ap.add_argument("--decay-band", default=str(RUNS / "decay_model_band_100.csv"))
    ap.add_argument("--bc-band", default=str(RUNS / "bc_nuisance_band_100.csv"))
    ap.add_argument("--breakdown", default=str(RUNS / "channel_breakdown_u2min.csv"))
    ap.add_argument("--l-int-fb", type=float, default=3000.0)
    ap.add_argument("--p-cut-mev", type=int, default=100)
    ap.add_argument("--out-dir", default=str(RUNS / "v1_bundle"))
    a = ap.parse_args(argv)

    central_csv = a.central_csv or str(RUNS / a.run / "analysis_exact_100" / "hnl_sensitivity.csv")
    cen = pd.read_csv(central_csv)   # keep insensitive rows so close_island can pinch the dome
    have_fonll = Path(a.fonll_band).exists()
    fonll = pd.read_csv(a.fonll_band) if have_fonll else None
    dm = pd.read_csv(a.decay_band)
    bc = pd.read_csv(a.bc_band)

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.6), sharey=True)
    for ax, fl in zip(axes, ["Ue", "Umu", "Utau"]):
        full = cen[cen.flavor == fl].sort_values("mass_GeV")
        sens = full[full.has_sensitivity == True]                            # noqa: E712
        cm = sens["mass_GeV"].to_numpy(float)
        c_min = sens["u2_min"].to_numpy(float)
        c_max = sens["u2_max"].to_numpy(float)
        # Rounded closure nose (same pinch _plot_single_panel draws); the bands
        # collapse onto its lower/upper arc so nothing pokes past the red island.
        vtx = _closure_pinch(full)
        arc_m = arc_lo = arc_hi = None
        if vtx is not None:
            arc_m, arc_lo, arc_hi = _closure_arc(cm[-1], c_min[-1], c_max[-1], vtx[0], vtx[1])
        tail_lo = (arc_m, arc_lo)   # for lower-edge bands
        tail_hi = (arc_m, arc_hi)   # for upper-edge bands

        # Red central island (+ its rounded closure); bands drawn below by us.
        _plot_single_panel(ax, cen, fl, is_leftmost=(fl == "Ue"), band_df=None,
                           close_island=True)

        # FONLL production band on both edges (hugs the dense central, follows nose).
        if have_fonll:
            bsub = fonll[fonll.flavor == fl].sort_values("mass_GeV")
            bm = bsub["mass_GeV"].to_numpy(float)
            labelled = False
            for col, c_edge, tail in (("u2_min", c_min, tail_lo), ("u2_max", c_max, tail_hi)):
                m, lo, hi = _dex_densify(
                    bm, bsub[f"{col}_central"].to_numpy(float),
                    bsub[f"{col}_band_lo"].to_numpy(float),
                    bsub[f"{col}_band_hi"].to_numpy(float),
                    cm, c_edge, tail=tail)
                _ribbon(ax, m, lo, hi, "orange",
                        None if labelled else "FONLL theory band", zorder=4)
                labelled = True

        # Decay-model width band (upper edge) and Bc normalization band (lower edge).
        dmf = dm[dm.flavor == fl].sort_values("mass_GeV")
        m, lo, hi = _dex_densify(
            dmf["mass_GeV"].to_numpy(float), dmf["u2_max"].to_numpy(float),
            dmf["u2_max_dm_lo"].to_numpy(float), dmf["u2_max_dm_hi"].to_numpy(float),
            cm, c_max, tail=tail_hi)
        _ribbon(ax, m, lo, hi, "steelblue", "decay-model band (upper)", zorder=4)
        bcf = bc[bc.flavor == fl].sort_values("mass_GeV")
        m, lo, hi = _dex_densify(
            bcf["mass_GeV"].to_numpy(float), bcf["u2_min"].to_numpy(float),
            bcf["u2_min_bc_lo"].to_numpy(float), bcf["u2_min_bc_hi"].to_numpy(float),
            cm, c_min, tail=tail_lo)
        _ribbon(ax, m, lo, hi, "teal", "Bc band (lower)", zorder=4)

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
