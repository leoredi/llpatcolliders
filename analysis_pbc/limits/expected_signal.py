"""
limits/expected_signal.py

Reusable physics kernel for the HNL limit pipeline:

  - couplings_from_eps2(): map PBC benchmark -> (Ue², Umu², Utau²)
  - expected_signal_events(): compute N_sig using decay+detector selection
  - scan_eps2_for_mass(): scan eps² grid and find exclusion interval

This module intentionally contains NO file-discovery, caching, or multiprocessing
logic. Drivers like `limits/run.py` handle I/O and orchestration.
"""

from __future__ import annotations

from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd

from config.production_xsecs import get_parent_sigma_pb, get_parent_tau_br
from decay.decay_detector import DecayCache, DecaySelection, compute_decay_acceptance, compute_separation_pass_static
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


def expected_signal_events(
    geom_df: pd.DataFrame,
    mass_GeV: float,
    eps2: float,
    benchmark: str,
    lumi_fb: float,
    dirac: bool = False,
    separation_m: float | None = None,
    decay_seed: int = 12345,
    decay_cache: DecayCache | None = None,
    separation_pass: np.ndarray | None = None,
) -> float:
    """
    Compute expected signal events N_sig using per-parent counting and
    explicit decay + detector selection.

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
    separation_m : float
        Minimum charged-track separation at detector surface (metres).
    decay_seed : int
        RNG seed for decay sampling (for reproducibility).
    decay_cache : DecayCache, optional
        Precomputed decay sampling cache to reuse across eps² points.
    separation_pass : np.ndarray, optional
        Precomputed boolean array of which HNLs pass separation cut.
        If provided, skips expensive per-eps2 recomputation.
    """
    # 1) Couplings and HNL model
    Ue2, Umu2, Utau2 = couplings_from_eps2(eps2, benchmark)
    model = HNLModel(mass_GeV=mass_GeV, Ue2=Ue2, Umu2=Umu2, Utau2=Utau2)

    # Proper decay length in metres
    ctau0_m = model.ctau0_m

    # Production BRs per parent
    br_per_parent: Dict[int, float] = model.production_brs()

    if separation_m is None:
        raise ValueError("Decay-based efficiency is mandatory: provide separation_m.")

    # 2) Extract geometry arrays
    required_cols = [
        "parent_id",
        "weight",
        "beta_gamma",
        "hits_tube",
        "entry_distance",
        "path_length",
        "eta",
        "phi",
    ]
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
    lam = beta_gamma * ctau0_m
    n_clamped = np.sum(lam <= 1e-9)
    if n_clamped > 0:
        print(f"[WARN] {n_clamped}/{len(lam)} HNLs have λ=βγ·cτ₀ ≤ 1e-9 (clamped); check beta_gamma/ctau0 inputs")
    lam = np.where(lam <= 1e-9, 1e-9, lam)

    P_decay = np.zeros_like(lam, dtype=float)
    mask_hits = hits_tube & (length > 0)

    if np.any(mask_hits):
        arg_entry = -entry[mask_hits] / lam[mask_hits]
        arg_path = -length[mask_hits] / lam[mask_hits]
        # Numerically stable: exp(A) * (1 - exp(B)) = exp(A) * (-expm1(B))
        P_decay[mask_hits] = np.exp(arg_entry) * (-np.expm1(arg_path))

    # 4) Apply decay + detector selection (track separation)
    if separation_pass is None:
        # Fallback: compute per eps2 point (slow, for backwards compatibility)
        selection = DecaySelection(separation_m=float(separation_m), seed=decay_seed)
        separation_pass = compute_decay_acceptance(
            geom_df=geom_df,
            mass_GeV=mass_GeV,
            flavour=benchmark_to_flavour(benchmark),
            ctau0_m=ctau0_m,
            mesh=build_mesh_once(),
            selection=selection,
            decay_cache=decay_cache,
        )
    P_decay = P_decay * separation_pass.astype(float)

    # 5) Group by parent species (per-parent counting)
    parent_abs = np.abs(parent_id.astype(int))
    tau_parent_id = None
    if "tau_parent_id" in geom_df.columns:
        tau_parent_id = geom_df["tau_parent_id"].to_numpy(dtype=float)
        # Replace NaN/inf with 0 before cast to avoid RuntimeWarning (EW files lack this column)
        tau_parent_id = np.abs(np.nan_to_num(tau_parent_id, nan=0.0, posinf=0.0, neginf=0.0)).astype(int)

    if tau_parent_id is None:
        mask_from_tau = np.zeros(len(parent_abs), dtype=bool)
    else:
        mask_from_tau = (parent_abs == 15) & (tau_parent_id > 0)

    # Check for malformed data: parent_id == 15 (tau) but tau_parent_id == 0
    # This would use the legacy σ(τ) approximation which ignores B→τ contribution
    mask_tau_parent = (parent_abs == 15)
    if tau_parent_id is None:
        n_bad = int(np.sum(mask_tau_parent))
    else:
        n_bad = int(np.sum(mask_tau_parent & (tau_parent_id == 0)))
    if n_bad > 0:
        raise ValueError(
            f"Found {n_bad} HNLs with parent_id=15 (tau) but tau_parent_id=0 or missing. "
            f"This indicates malformed fromTau data — tau_parent_id must specify the grandfather meson (Ds/B). "
            f"Check production code or regenerate data."
        )

    unique_parents = np.unique(parent_abs[~mask_from_tau])
    total_expected = 0.0

    missing_br_pdgs = []
    missing_xsec_pdgs = []
    missing_tau_br_pdgs = []

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

    # 5b) Tau-decay production: parent -> tau nu, tau -> N X
    if np.any(mask_from_tau):
        assert tau_parent_id is not None  # Guaranteed by mask_from_tau logic above
        br_tau_to_N = br_per_parent.get(15, 0.0)
        if br_tau_to_N <= 0.0:
            missing_br_pdgs.append(15)
        else:
            tau_parents = np.unique(tau_parent_id[mask_from_tau])
            for pid in tau_parents:
                br_parent_to_tau = get_parent_tau_br(int(pid))
                if br_parent_to_tau <= 0.0:
                    missing_tau_br_pdgs.append(int(pid))
                    continue

                sigma_parent_pb = get_parent_sigma_pb(int(pid))
                if sigma_parent_pb <= 0.0:
                    missing_xsec_pdgs.append(int(pid))
                    continue

                mask_parent = mask_from_tau & (tau_parent_id == pid)
                w = weights[mask_parent]
                P = P_decay[mask_parent]
                w_sum = np.sum(w)
                if w_sum <= 0.0:
                    continue

                eff_parent = np.sum(w * P) / w_sum
                total_expected += lumi_fb * (sigma_parent_pb * 1e3) * br_parent_to_tau * br_tau_to_N * eff_parent

    # Diagnostics: log once per mass point (at first scan point)
    if missing_br_pdgs and eps2 == 1e-12:
        n_lost = int(np.sum(np.isin(parent_id, missing_br_pdgs)))
        print(
            f"[WARN] Mass {mass_GeV:.2f} GeV: {len(missing_br_pdgs)} parent PDG(s) have no HNLCalc BR: {missing_br_pdgs}"
        )
        print(f"       → Discarding {n_lost} events (silent data loss)")

    if missing_tau_br_pdgs and eps2 == 1e-12:
        n_lost = int(np.sum(np.isin(tau_parent_id, missing_tau_br_pdgs))) if tau_parent_id is not None else 0
        print(
            f"[WARN] Mass {mass_GeV:.2f} GeV: {len(missing_tau_br_pdgs)} tau-parent PDG(s) have no SM τ BR: {missing_tau_br_pdgs}"
        )
        print(f"       → Discarding {n_lost} events (silent data loss)")

    if missing_xsec_pdgs and eps2 == 1e-12:
        n_lost = int(np.sum(np.isin(parent_id, missing_xsec_pdgs)))
        if tau_parent_id is not None:
            n_lost += int(np.sum(np.isin(tau_parent_id, missing_xsec_pdgs)))
        print(
            f"[WARN] Mass {mass_GeV:.2f} GeV: {len(missing_xsec_pdgs)} parent PDG(s) have no cross-section: {missing_xsec_pdgs}"
        )
        print(f"       → Discarding {n_lost} events (silent data loss)")

    if dirac:
        total_expected *= 2.0

    return float(total_expected)


_MESH_CACHE = None


def build_mesh_once():
    global _MESH_CACHE
    if _MESH_CACHE is None:
        from geometry.per_parent_efficiency import build_drainage_gallery_mesh

        _MESH_CACHE = build_drainage_gallery_mesh()
    return _MESH_CACHE


def benchmark_to_flavour(benchmark: str) -> str:
    if benchmark == "100":
        return "electron"
    if benchmark == "010":
        return "muon"
    if benchmark == "001":
        return "tau"
    raise ValueError(f"Unsupported benchmark: {benchmark} (use '100','010','001').")


def scan_eps2_for_mass(
    geom_df: pd.DataFrame,
    mass_GeV: float,
    benchmark: str,
    lumi_fb: float,
    N_limit: float = 2.996,
    dirac: bool = False,
    separation_m: float | None = None,
    decay_seed: int = 12345,
) -> Tuple[np.ndarray, np.ndarray, Optional[float], Optional[float]]:
    eps2_grid = np.logspace(-12, -2, 100)
    Nsig = np.zeros_like(eps2_grid, dtype=float)

    # Build decay cache and separation_pass ONCE and reuse across all eps2 points
    decay_cache = None
    separation_pass = None
    if separation_m is not None:
        from decay.decay_detector import DecaySelection, build_decay_cache
        flavour = benchmark_to_flavour(benchmark)
        selection = DecaySelection(separation_m=float(separation_m), seed=decay_seed)
        decay_cache = build_decay_cache(geom_df, mass_GeV, flavour, selection)
        # Precompute separation using static method (ctau-independent approximation)
        separation_pass = compute_separation_pass_static(
            geom_df, decay_cache, build_mesh_once(), float(separation_m)
        )

    for i, eps2 in enumerate(eps2_grid):
        Nsig[i] = expected_signal_events(
            geom_df=geom_df,
            mass_GeV=mass_GeV,
            eps2=float(eps2),
            benchmark=benchmark,
            lumi_fb=lumi_fb,
            dirac=dirac,
            separation_m=separation_m,
            decay_seed=decay_seed,
            decay_cache=decay_cache,
            separation_pass=separation_pass,
        )

    mask = Nsig >= N_limit
    if not np.any(mask):
        return eps2_grid, Nsig, None, None

    idx_above = np.where(mask)[0]
    eps2_min = float(eps2_grid[idx_above[0]])
    eps2_max = float(eps2_grid[idx_above[-1]])
    return eps2_grid, Nsig, eps2_min, eps2_max
