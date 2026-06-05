#!/usr/bin/env python3
"""Ds/B+ -> tau -> HNL production CSV generation."""

import random
import sys
import numpy as np
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from config_mass_grid import MASS_GRID
from production.constants import (
    M_TAU, M_DS, M_BPLUS,
    FRAG_B, FRAG_C,
)
from production.fonll.fonll_parser import get_sigma_total
from production.fonll.meson_sampler import sample_meson_4vectors
from production.decay_engine.kinematics import decay_2body
from production.decay_engine.tau_decay import (
    init_hnlcalc, compute_tau_production_br_components, sample_hnl_from_tau,
    TAU_2BODY_ASYMMETRY,
)
from production.io import llp_csv_path, write_empty_csv, write_llp_csv
from production.paths import LLP_VECTORS_DIR

OUTPUT_BASE = LLP_VECTORS_DIR

N_POOL = 100_000

BR_DS_TAU_NU = 5.35e-2
BR_BPLUS_TAU_NU = 1.09e-4


def build_tau_pool(n_pool, rng):
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
    hnl = init_hnlcalc(flavor)

    if len(tau_w) == 0:
        for m_N in masses:
            write_empty_csv(llp_csv_path(flavor, "induced_tau", m_N, base=OUTPUT_BASE))
        return

    tau_E = tau_4v[:, 0]
    tau_px = tau_4v[:, 1]
    tau_py = tau_4v[:, 2]
    tau_pz = tau_4v[:, 3]

    for m_N in masses:
        csv_path = llp_csv_path(flavor, "induced_tau", m_N, base=OUTPUT_BASE)
        if m_N >= M_TAU:
            write_empty_csv(csv_path)
            continue

        br_2body, br_3body, br_total = compute_tau_production_br_components(hnl, m_N)
        if br_total <= 0:
            write_empty_csv(csv_path)
            continue

        hnl_4v, _, _ = sample_hnl_from_tau(
            tau_E, tau_px, tau_py, tau_pz, m_N,
            br_2body, br_3body, br_total, rng,
        )
        weights = tau_w * br_total
        write_llp_csv(csv_path, weights, hnl_4v[:, 0], hnl_4v[:, 1], hnl_4v[:, 2], hnl_4v[:, 3])
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
    random.seed(args.seed)
    masses = args.masses if args.masses else MASS_GRID

    print("Building tau pool from Ds and B+ ...")
    tau_4v, tau_w = build_tau_pool(args.n_pool, rng)

    for flavor in args.flavor:
        print(f"\n{'='*60}\nFlavor: {flavor}\n{'='*60}")
        process_flavor(flavor, tau_4v, tau_w, masses, rng)

    print("\nDone.")


if __name__ == "__main__":
    main()
