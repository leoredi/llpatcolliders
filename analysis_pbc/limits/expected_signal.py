"""
limits/expected_signal.py

Reusable physics kernel for the HNL limit pipeline:

  - couplings_from_eps2(): map PBC benchmark -> (Ue², Umu², Utau²)
  - expected_signal_events(): compute N_sig for a geometry dataframe
  - scan_eps2_for_mass(): scan eps² grid and find exclusion interval

This module intentionally contains NO file-discovery, caching, or multiprocessing
logic. Drivers like `limits/run.py` handle I/O and orchestration.
"""

from __future__ import annotations

from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd

from config.production_xsecs import get_parent_sigma_pb
from models.hnl_model_hnlcalc import HNLModel


def couplings_from_eps2(eps2: float, benchmark: str) -> Tuple[float, float, float]:
    """
    Map PBC benchmark (string) to (Ue², Umu², Utau²).
    """
    if benchmark == "100":
        return eps2, 0.0, 0.0
    if benchmark == "010":
        return 0.0, eps2, 0.0
    if benchmark == "001":
        return 0.0, 0.0, eps2
    raise ValueError(f"Unsupported benchmark: {benchmark} (use '100','010','001').")


def _decay_probabilities(
    *,
    beta_gamma: np.ndarray,
    hits_tube: np.ndarray,
    entry_distance_m: np.ndarray,
    path_length_m: np.ndarray,
    ctau0_m: float,
) -> np.ndarray:
    lam = beta_gamma * float(ctau0_m)
    lam = np.where(lam <= 1e-9, 1e-9, lam)  # Prevent divide by zero

    P_decay = np.zeros_like(lam, dtype=float)
    mask_hits = hits_tube & (path_length_m > 0)

    if np.any(mask_hits):
        arg_entry = -entry_distance_m[mask_hits] / lam[mask_hits]
        arg_path = -path_length_m[mask_hits] / lam[mask_hits]
        P_decay[mask_hits] = np.exp(arg_entry) * (-np.expm1(arg_path))

    return P_decay


def expected_signal_events(
    geom_df: pd.DataFrame,
    mass_GeV: float,
    eps2: float,
    benchmark: str,
    lumi_fb: float,
    dirac: bool = False,
) -> float:
    """
    Compute expected signal events N_sig using per-parent counting.

        N_sig = Σ_parents [ L × σ_parent × BR(parent→ℓN) × ε_geom(parent) ]

    Parameters
    ----------
    geom_df : pd.DataFrame
        Geometry dataframe with columns: parent_id, weight, beta_gamma,
        hits_tube, entry_distance, path_length (one row per HNL)
    mass_GeV : float
        HNL mass in GeV
    eps2 : float
        Total coupling squared |U|²
    benchmark : str
        Coupling pattern: "100" (electron), "010" (muon), "001" (tau)
    lumi_fb : float
        Integrated luminosity in fb⁻¹
    dirac : bool
        If True, multiply yield by 2 for Dirac HNL interpretation (N ≠ N̄).
        Default False assumes Majorana HNL.
    """
    # 1) Couplings and HNL model
    Ue2, Umu2, Utau2 = couplings_from_eps2(eps2, benchmark)
    model = HNLModel(mass_GeV=mass_GeV, Ue2=Ue2, Umu2=Umu2, Utau2=Utau2)

    # Proper decay length in metres
    ctau0_m = model.ctau0_m

    # Production BRs per parent
    br_per_parent: Dict[int, float] = model.production_brs()

    # 2) Extract geometry arrays
    required_cols = ["parent_id", "weight", "beta_gamma", "hits_tube", "entry_distance", "path_length"]
    missing_cols = [c for c in required_cols if c not in geom_df.columns]
    if missing_cols:
        print(f"[WARN] geom_df missing columns {missing_cols}. Returning N_sig=0.")
        return 0.0

    if len(geom_df) == 0:
        return 0.0

    parent_id = geom_df["parent_id"].to_numpy()
    weights = geom_df["weight"].to_numpy(dtype=float)
    beta_gamma = geom_df["beta_gamma"].to_numpy(dtype=float)
    hits_tube = geom_df["hits_tube"].to_numpy(dtype=bool)
    entry = geom_df["entry_distance"].to_numpy(dtype=float)
    length = geom_df["path_length"].to_numpy(dtype=float)

    # --- SAFETY CHECK: DROP BAD PARENT IDs ---
    mask_valid = np.isfinite(parent_id)
    if not np.all(mask_valid):
        parent_id = parent_id[mask_valid]
        weights = weights[mask_valid]
        beta_gamma = beta_gamma[mask_valid]
        hits_tube = hits_tube[mask_valid]
        entry = entry[mask_valid]
        length = length[mask_valid]

    if len(parent_id) == 0:
        return 0.0

    # 3) Compute P_decay_i for all HNLs
    P_decay = _decay_probabilities(
        beta_gamma=beta_gamma,
        hits_tube=hits_tube,
        entry_distance_m=entry,
        path_length_m=length,
        ctau0_m=ctau0_m,
    )

    # 4) Group by parent species (per-parent counting)
    unique_parents = np.unique(np.abs(parent_id.astype(int)))
    total_expected = 0.0

    missing_br_pdgs = []
    missing_xsec_pdgs = []

    for pid in unique_parents:
        BR_parent = br_per_parent.get(int(pid), 0.0)
        if BR_parent <= 0.0:
            missing_br_pdgs.append(int(pid))
            continue

        sigma_parent_pb = get_parent_sigma_pb(int(pid))
        if sigma_parent_pb <= 0.0:
            missing_xsec_pdgs.append(int(pid))
            continue

        mask_parent = np.abs(parent_id) == pid
        w = weights[mask_parent]
        P = P_decay[mask_parent]
        w_sum = np.sum(w)
        if w_sum <= 0.0:
            continue

        eff_parent = np.sum(w * P) / w_sum
        # N = L(fb^-1) × σ(pb) × (1e3 fb/pb) × BR × ε
        total_expected += lumi_fb * (sigma_parent_pb * 1e3) * BR_parent * eff_parent

    # Diagnostics: log once per mass point (at first scan point)
    if missing_br_pdgs and eps2 == 1e-12:
        n_lost = int(np.sum(np.isin(parent_id, missing_br_pdgs)))
        print(
            f"[WARN] Mass {mass_GeV:.2f} GeV: {len(missing_br_pdgs)} parent PDG(s) have no HNLCalc BR: {missing_br_pdgs}"
        )
        print(f"       → Discarding {n_lost} events (silent data loss)")

    if missing_xsec_pdgs and eps2 == 1e-12:
        n_lost = int(np.sum(np.isin(parent_id, missing_xsec_pdgs)))
        print(
            f"[WARN] Mass {mass_GeV:.2f} GeV: {len(missing_xsec_pdgs)} parent PDG(s) have no cross-section: {missing_xsec_pdgs}"
        )
        print(f"       → Discarding {n_lost} events (silent data loss)")

    if dirac:
        total_expected *= 2.0

    return float(total_expected)


def expected_signal_events_alp(
    geom_df: pd.DataFrame,
    mass_GeV: float,
    fa_GeV: float,
    benchmark: str,
    lumi_fb: float,
    *,
    production_mode: str,
    visible_br_override: Optional[float] = None,
    # Higgs channels
    br_h_aa: Optional[float] = None,
    lambda_aha_GeV: Optional[float] = None,
    br_h_Za: Optional[float] = None,
    C_aZh: Optional[float] = None,
    Lambda_GeV: float = 1000.0,
    # Z channel
    C_gamma_Z: float = 1.0,
    # FCNC channels
    C_WW: float = 1.0,
) -> float:
    """
    Compute expected signal events for ALPs, reusing the same geometry kernel.

    The geometry dataframe is assumed to represent an ALP sample where the
    generator forced the parent→ALP decay with BR=1.0. This function applies
    the physical production branching ratios analytically, in the same spirit
    as the HNL pipeline.
    """
    from models.alp_model import ALPModel
    from production.alp_production import br_B_to_Ka, br_Z_to_gamma_a, br_higgs_to_Za, br_higgs_to_aa

    alp = ALPModel(mass_GeV=mass_GeV, fa_GeV=fa_GeV, benchmark=benchmark)
    ctau0_m = alp.ctau0_m
    br_vis = float(visible_br_override) if visible_br_override is not None else alp.visible_br

    required_cols = ["parent_id", "weight", "beta_gamma", "hits_tube", "entry_distance", "path_length"]
    missing_cols = [c for c in required_cols if c not in geom_df.columns]
    if missing_cols:
        print(f"[WARN] geom_df missing columns {missing_cols}. Returning N_sig=0.")
        return 0.0
    if len(geom_df) == 0:
        return 0.0

    parent_id = geom_df["parent_id"].to_numpy()
    weights = geom_df["weight"].to_numpy(dtype=float)
    beta_gamma = geom_df["beta_gamma"].to_numpy(dtype=float)
    hits_tube = geom_df["hits_tube"].to_numpy(dtype=bool)
    entry = geom_df["entry_distance"].to_numpy(dtype=float)
    length = geom_df["path_length"].to_numpy(dtype=float)

    mask_valid = np.isfinite(parent_id)
    if not np.all(mask_valid):
        parent_id = parent_id[mask_valid]
        weights = weights[mask_valid]
        beta_gamma = beta_gamma[mask_valid]
        hits_tube = hits_tube[mask_valid]
        entry = entry[mask_valid]
        length = length[mask_valid]
    if len(parent_id) == 0:
        return 0.0

    P_decay = _decay_probabilities(
        beta_gamma=beta_gamma,
        hits_tube=hits_tube,
        entry_distance_m=entry,
        path_length_m=length,
        ctau0_m=ctau0_m,
    )

    unique_parents = np.unique(np.abs(parent_id.astype(int)))
    total_expected = 0.0
    missing_xsec_pdgs = []
    missing_br_pdgs = []

    production_mode = str(production_mode)

    for pid in unique_parents:
        sigma_parent_pb = get_parent_sigma_pb(int(pid), particle="alp")
        if sigma_parent_pb <= 0.0:
            missing_xsec_pdgs.append(int(pid))
            continue

        if production_mode == "Z_to_gamma_a":
            br_eff = br_Z_to_gamma_a(mass_GeV, C_gamma_Z, fa_GeV) if int(pid) == 23 else 0.0
        elif production_mode == "B_to_Ka":
            if int(pid) == 521:
                br_eff = br_B_to_Ka(mass_GeV, fa_GeV, C_WW, meson="B+")
            elif int(pid) == 511:
                br_eff = br_B_to_Ka(mass_GeV, fa_GeV, C_WW, meson="B0")
            elif int(pid) == 531:
                br_eff = br_B_to_Ka(mass_GeV, fa_GeV, C_WW, meson="Bs")
            else:
                br_eff = 0.0
        elif production_mode == "h_to_aa":
            if int(pid) != 25:
                br_eff = 0.0
            elif br_h_aa is not None:
                br_eff = 2.0 * float(br_h_aa)
            elif lambda_aha_GeV is not None:
                br_eff = 2.0 * br_higgs_to_aa(mass_GeV, float(lambda_aha_GeV))
            else:
                br_eff = 0.0
        elif production_mode == "h_to_Za":
            if int(pid) != 25:
                br_eff = 0.0
            elif br_h_Za is not None:
                br_eff = float(br_h_Za)
            elif C_aZh is not None:
                br_eff = br_higgs_to_Za(mass_GeV, float(C_aZh), Lambda_GeV=Lambda_GeV)
            else:
                br_eff = 0.0
        else:
            raise ValueError(
                "Unknown ALP production_mode: "
                f"{production_mode} (use h_to_aa, h_to_Za, Z_to_gamma_a, B_to_Ka)"
            )

        if br_eff <= 0.0:
            missing_br_pdgs.append(int(pid))
            continue

        mask_parent = np.abs(parent_id) == pid
        w = weights[mask_parent]
        P = P_decay[mask_parent]
        w_sum = float(np.sum(w))
        if w_sum <= 0.0:
            continue

        eff_parent = float(np.sum(w * P) / w_sum)
        total_expected += lumi_fb * (sigma_parent_pb * 1e3) * float(br_eff) * eff_parent

    if missing_xsec_pdgs:
        print(f"[WARN] ALP m={mass_GeV:.2f} GeV: missing xsec for PDG(s): {sorted(set(missing_xsec_pdgs))}")
    if missing_br_pdgs and fa_GeV == 1e6:
        print(f"[INFO] ALP m={mass_GeV:.2f} GeV: production BR=0 for PDG(s): {sorted(set(missing_br_pdgs))}")

    return float(total_expected * br_vis)


def scan_eps2_for_mass(
    geom_df: pd.DataFrame,
    mass_GeV: float,
    benchmark: str,
    lumi_fb: float,
    N_limit: float = 2.996,
    dirac: bool = False,
) -> Tuple[np.ndarray, np.ndarray, Optional[float], Optional[float]]:
    eps2_grid = np.logspace(-12, -2, 100)
    Nsig = np.zeros_like(eps2_grid, dtype=float)

    for i, eps2 in enumerate(eps2_grid):
        Nsig[i] = expected_signal_events(
            geom_df=geom_df,
            mass_GeV=mass_GeV,
            eps2=float(eps2),
            benchmark=benchmark,
            lumi_fb=lumi_fb,
            dirac=dirac,
        )

    mask = Nsig >= N_limit
    if not np.any(mask):
        return eps2_grid, Nsig, None, None

    idx_above = np.where(mask)[0]
    eps2_min = float(eps2_grid[idx_above[0]])
    eps2_max = float(eps2_grid[idx_above[-1]])
    return eps2_grid, Nsig, eps2_min, eps2_max
