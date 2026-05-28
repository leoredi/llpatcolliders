#!/usr/bin/env python3
"""
production/decay_engine/generate_induced_tau.py

Driver for induced-tau HNL production:

    pp -> Ds + X         (FONLL charm)
    Ds -> tau nu_tau     (BR = 5.35e-2, 2-body)
    tau -> N + X         (HNLCalc 2-body and 3-body)

    pp -> B+ + X         (FONLL bottom)
    B+ -> tau nu_tau     (BR = 1.09e-4, 2-body)
    tau -> N + X

Note: B0 -> tau nu is helicity-suppressed in SM and excluded.
W -> tau nu is the dominant tau source at LHC but is W-mediated, deferred
to the W/Z production PR.

Weight chain per HNL 4-vector i (from Ds, similar for B+):

    w_i = 2 * sigma_FONLL_charm * f_Ds * BR(Ds->tau nu) * BR(tau->N+X) / N_tau_sampled

The factor 2 is the particle + antiparticle FONLL convention.

Output: output/llp_4vectors/{Ue,Umu,Utau}/tau/mN_{mass}.csv
Format: headerless, 5 columns: weight,E,px,py,pz
"""

import sys
import numpy as np
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "vendored" / "HNLCalc"))

from config_mass_grid import MASS_GRID, format_mass_for_filename
from production.constants import (
    M_TAU, M_DS, M_BPLUS,
    FRAG_B, FRAG_C,
)
from production.fonll.fonll_parser import get_sigma_total
from production.fonll.meson_sampler import sample_meson_4vectors
from production.decay_engine.kinematics import decay_2body, decay_3body_flat

OUTPUT_BASE = PROJECT_ROOT / "output" / "llp_4vectors"

N_POOL = 100_000

# Measured branching ratios for the parent meson -> tau nu_tau (PDG 2024)
BR_DS_TAU_NU = 5.35e-2
BR_BPLUS_TAU_NU = 1.09e-4

# Tau 2-body hadronic channels into N (parent_pdg = 15 (tau-))
# daughter is meson (anti-meson for tau- decay)
TAU_2BODY_MESON_PDGS = [211, 321, 213, 323]   # pi+, K+, rho+, K*+ (HNLCalc convention)


def _init_hnlcalc(flavor):
    """Initialize HNLCalc with unit coupling for the given flavor."""
    from HNLCalc import HNLCalc
    if flavor == "Ue":
        return HNLCalc(ve=1, vmu=0, vtau=0)
    elif flavor == "Umu":
        return HNLCalc(ve=0, vmu=1, vtau=0)
    elif flavor == "Utau":
        return HNLCalc(ve=0, vmu=0, vtau=1)
    else:
        raise ValueError(f"Unknown flavor: {flavor}")


def _eval_tau_2body_br(hnl, meson_pdg, m_N):
    """BR(tau -> meson N) at U^2 = 1 via HNLCalc."""
    try:
        br_expr = hnl.get_2body_br_tau(15, meson_pdg)
        mass = m_N        # noqa: F841 - referenced by eval'd HNLCalc expression
        coupling = 1.0    # noqa: F841
        br_val = eval(br_expr)
        if np.isnan(br_val) or br_val < 0:
            return 0.0
        return float(br_val)
    except Exception:
        return 0.0


def _eval_tau_3body_br(hnl, lep_pid, nu_pid, m_N):
    """BR(tau- -> lep- nu N) at U^2 = 1 via HNLCalc, integrated over phase space."""
    try:
        dbr = hnl.get_3body_dbr_tau(15, -lep_pid, nu_pid)
        m_lep = hnl.masses(lep_pid)
        br_val = hnl.integrate_3body_br(
            dbr, m_N, M_TAU, m_lep, 0.0,
            coupling=1.0, nsample=500, integration="dE",
        )
        if br_val is None or np.isnan(br_val) or br_val < 0:
            return 0.0
        return float(br_val)
    except Exception:
        return 0.0


def compute_tau_production_br_components(hnl, m_N):
    """
    Total BR(tau -> N + X) at U^2 = 1, plus per-channel breakdown for sampling.

    Returns
    -------
    br_2body_channels : list of (meson_mass, br_value)
        2-body hadronic channels with non-zero BR.
    br_3body_channels : list of (lep_mass, br_value)
        3-body leptonic channels (tau- -> lep- nu N) with non-zero BR.
    br_total : float
        Sum of all BRs.
    """
    br_2body_channels = []
    for meson_pdg in TAU_2BODY_MESON_PDGS:
        if m_N >= M_TAU - hnl.masses(meson_pdg):
            continue
        br = _eval_tau_2body_br(hnl, meson_pdg, m_N)
        if br > 0:
            br_2body_channels.append((hnl.masses(meson_pdg), br))

    br_3body_channels = []
    for lep_pid in [11, 13]:
        m_lep = hnl.masses(lep_pid)
        if m_N >= M_TAU - m_lep:
            continue
        # tau- -> lep- nu_tau N   (neutrino is nu_tau, pdg=16)
        br_nt = _eval_tau_3body_br(hnl, lep_pid, 16, m_N)
        if br_nt > 0:
            br_3body_channels.append((m_lep, br_nt))
        # tau- -> lep- nu_lep_bar N   (neutrino is anti-nu_lep, pdg = lep_pid+1)
        br_nl = _eval_tau_3body_br(hnl, lep_pid, lep_pid + 1, m_N)
        if br_nl > 0:
            br_3body_channels.append((m_lep, br_nl))

    br_total = sum(br for _, br in br_2body_channels) + sum(br for _, br in br_3body_channels)
    return br_2body_channels, br_3body_channels, br_total


def _sample_hnl_from_tau(tau_E, tau_px, tau_py, tau_pz, m_N,
                         br_2body_channels, br_3body_channels, br_total, rng):
    """Decay each tau to (N + X) with channel selection weighted by BR."""
    n_events = len(tau_E)
    hnl_4v = np.empty((n_events, 4))

    br_2body_total = sum(br for _, br in br_2body_channels)
    p_2body = br_2body_total / br_total if br_total > 0 else 0.0
    use_2body = rng.random(n_events) < p_2body
    use_3body = ~use_2body

    # 2-body: tau -> meson + N. Pick channel weighted by per-channel BR.
    if use_2body.any() and br_2body_channels:
        br_arr = np.array([br for _, br in br_2body_channels], dtype=float)
        prob = br_arr / br_arr.sum()
        ch_idx = rng.choice(len(br_2body_channels), size=use_2body.sum(), p=prob)
        evt_idx = np.where(use_2body)[0]
        for k in np.unique(ch_idx):
            m_meson, _ = br_2body_channels[int(k)]
            sel = evt_idx[ch_idx == k]
            _, hnl_2b = decay_2body(
                tau_E[sel], tau_px[sel], tau_py[sel], tau_pz[sel],
                M_TAU, m_meson, m_N, rng=rng,
            )
            hnl_4v[sel] = hnl_2b

    # 3-body: tau -> lep + nu + N. Flat phase space.
    if use_3body.any() and br_3body_channels:
        br_arr = np.array([br for _, br in br_3body_channels], dtype=float)
        prob = br_arr / br_arr.sum()
        ch_idx = rng.choice(len(br_3body_channels), size=use_3body.sum(), p=prob)
        evt_idx = np.where(use_3body)[0]
        for k in np.unique(ch_idx):
            m_lep, _ = br_3body_channels[int(k)]
            sel = evt_idx[ch_idx == k]
            # tau -> lep (m_lep), nu (0), N (m_N)
            _, _, hnl_3b = decay_3body_flat(
                tau_E[sel], tau_px[sel], tau_py[sel], tau_pz[sel],
                M_TAU, m_lep, 0.0, m_N, rng=rng,
            )
            hnl_4v[sel] = hnl_3b

    return hnl_4v, use_2body, use_3body


def build_tau_pool(n_pool, rng):
    """
    Build the tau pool by chaining Ds->tau nu and B+ -> tau nu off FONLL mesons.

    Returns
    -------
    tau_4v : ndarray, shape (M, 4)
    tau_w  : ndarray, shape (M,)   per-event weight in pb (sigma * BR / N_per_source)
    """
    # --- Ds source ---
    charm_pool = sample_meson_4vectors(n_pool, "charm", rng=rng)
    sigma_charm = get_sigma_total("charm")
    ds_mask = charm_pool['species_pdg'] == 431
    n_ds = int(ds_mask.sum())
    if n_ds > 0:
        ds_E = charm_pool['E'][ds_mask]
        ds_px = charm_pool['px'][ds_mask]
        ds_py = charm_pool['py'][ds_mask]
        ds_pz = charm_pool['pz'][ds_mask]
        tau_4v_ds, _ = decay_2body(ds_E, ds_px, ds_py, ds_pz,
                                   M_DS, M_TAU, 0.0, rng=rng)
        w_ds = 2.0 * sigma_charm * FRAG_C[431] * BR_DS_TAU_NU / n_ds
        tau_w_ds = np.full(n_ds, w_ds)
    else:
        tau_4v_ds = np.empty((0, 4))
        tau_w_ds = np.empty(0)

    # --- B+ source ---
    bottom_pool = sample_meson_4vectors(n_pool, "bottom", rng=rng)
    sigma_bottom = get_sigma_total("bottom")
    bp_mask = bottom_pool['species_pdg'] == 521
    n_bp = int(bp_mask.sum())
    if n_bp > 0:
        bp_E = bottom_pool['E'][bp_mask]
        bp_px = bottom_pool['px'][bp_mask]
        bp_py = bottom_pool['py'][bp_mask]
        bp_pz = bottom_pool['pz'][bp_mask]
        tau_4v_bp, _ = decay_2body(bp_E, bp_px, bp_py, bp_pz,
                                   M_BPLUS, M_TAU, 0.0, rng=rng)
        w_bp = 2.0 * sigma_bottom * FRAG_B[521] * BR_BPLUS_TAU_NU / n_bp
        tau_w_bp = np.full(n_bp, w_bp)
    else:
        tau_4v_bp = np.empty((0, 4))
        tau_w_bp = np.empty(0)

    tau_4v = np.vstack([tau_4v_ds, tau_4v_bp])
    tau_w = np.concatenate([tau_w_ds, tau_w_bp])

    print(f"  Ds -> tau nu:  {n_ds} tau, w_sum = {tau_w_ds.sum():.3e} pb")
    print(f"  B+ -> tau nu:  {n_bp} tau, w_sum = {tau_w_bp.sum():.3e} pb")
    print(f"  Tau pool:      {len(tau_w)} taus, w_sum = {tau_w.sum():.3e} pb")
    return tau_4v, tau_w


def process_flavor(flavor, tau_4v, tau_w, masses, rng):
    """Decay tau pool into N+X for one flavor across all mass points."""
    hnl = _init_hnlcalc(flavor)
    out_dir = OUTPUT_BASE / flavor / "tau"
    out_dir.mkdir(parents=True, exist_ok=True)

    if len(tau_w) == 0:
        for m_N in masses:
            (out_dir / f"mN_{format_mass_for_filename(m_N)}.csv").write_text("")
        return

    tau_E = tau_4v[:, 0]
    tau_px = tau_4v[:, 1]
    tau_py = tau_4v[:, 2]
    tau_pz = tau_4v[:, 3]

    for m_N in masses:
        mass_label = format_mass_for_filename(m_N)
        csv_path = out_dir / f"mN_{mass_label}.csv"
        if m_N >= M_TAU:
            csv_path.write_text("")
            continue

        br_2body, br_3body, br_total = compute_tau_production_br_components(hnl, m_N)
        if br_total <= 0:
            csv_path.write_text("")
            continue

        hnl_4v, _, _ = _sample_hnl_from_tau(
            tau_E, tau_px, tau_py, tau_pz, m_N,
            br_2body, br_3body, br_total, rng,
        )
        weights = tau_w * br_total
        data = np.column_stack([weights, hnl_4v[:, 0], hnl_4v[:, 1], hnl_4v[:, 2], hnl_4v[:, 3]])
        np.savetxt(csv_path, data, delimiter=",", fmt="%.8e")
        print(f"    {csv_path.name}: {len(weights)} events, "
              f"BR_total={br_total:.3e}, w_sum={weights.sum():.3e}")


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Generate induced-tau HNL CSVs (Ds, B+ -> tau -> N)")
    parser.add_argument("--flavor", choices=["Ue", "Umu", "Utau"], nargs="+",
                        default=["Ue", "Umu", "Utau"])
    parser.add_argument("--n-pool", type=int, default=N_POOL)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--masses", type=float, nargs="+", default=None,
                        help="Custom mass list (default: full grid)")
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    masses = args.masses if args.masses else MASS_GRID
    masses = [m for m in masses if m < M_TAU]  # tau decay closes at m_tau

    print("Building tau pool from Ds and B+ ...")
    tau_4v, tau_w = build_tau_pool(args.n_pool, rng)

    for flavor in args.flavor:
        print(f"\n{'='*60}\nFlavor: {flavor}\n{'='*60}")
        process_flavor(flavor, tau_4v, tau_w, masses, rng)

    print("\nDone.")


if __name__ == "__main__":
    main()
