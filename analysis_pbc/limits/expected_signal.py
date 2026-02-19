
from __future__ import annotations

from typing import Dict

import numpy as np
import pandas as pd

from config.production_xsecs import get_parent_sigma_pb, get_parent_tau_br
from decay.decay_detector import (
    DecayCache,
    DecaySelection,
    _normalize_separation_policy,
    compute_decay_acceptance,
)
from geometry.per_parent_efficiency import (
    GeometryConfig,
    build_drainage_gallery_mesh,
    geometry_tag,
    normalize_geometry_config,
)
from hnl_models.hnl_model_hnlcalc import HNLModel
from limits.timing_utils import _time_block


def couplings_from_eps2(eps2: float, benchmark: str) -> Tuple[float, float, float]:
    if benchmark == "100":
        return eps2, 0.0, 0.0
    if benchmark == "010":
        return 0.0, eps2, 0.0
    if benchmark == "001":
        return 0.0, 0.0, eps2
    raise ValueError(f"Unsupported benchmark: {benchmark} (use '100','010','001').")


def _normalize_qcd_mode(series: pd.Series | None, n: int) -> np.ndarray:
    if series is None:
        return np.full(n, "auto", dtype=object)
    mode = series.astype(str).str.strip()
    mode = mode.replace({"": "auto", "nan": "auto", "None": "auto"})
    return mode.to_numpy(dtype=object)


def _normalize_sigma_gen_pb(series: pd.Series | None, n: int) -> np.ndarray:
    if series is None:
        return np.full(n, np.nan, dtype=float)
    return pd.to_numeric(series, errors="coerce").to_numpy(dtype=float)


def _context_groups(qcd_mode: np.ndarray, sigma_gen_pb: np.ndarray) -> dict[tuple[str, str], np.ndarray]:
    grouped: dict[tuple[str, str], list[int]] = {}
    for idx, (mode_raw, sigma_raw) in enumerate(zip(qcd_mode, sigma_gen_pb)):
        mode = str(mode_raw) if mode_raw is not None else "auto"
        sigma_key = f"{float(sigma_raw):.12g}" if np.isfinite(sigma_raw) else "nan"
        key = (mode, sigma_key)
        if key not in grouped:
            grouped[key] = []
        grouped[key].append(idx)
    return {k: np.array(v, dtype=int) for k, v in grouped.items()}


def _resolve_parent_sigma_pb(
    pid: int,
    qcd_mode: str,
    sigma_gen_pb: float,
) -> float:
    qcd_mode = qcd_mode or "auto"

    if qcd_mode == "hardccbar":
        if np.isfinite(sigma_gen_pb) and sigma_gen_pb > 0.0:
            return get_parent_sigma_pb(pid, sigma_ccbar_pb=float(sigma_gen_pb))
        raise ValueError(
            f"Hard-QCD mode '{qcd_mode}' for parent {pid} has no valid sigma_gen_pb "
            f"({sigma_gen_pb}). Falling back to inclusive cross-section can strongly "
            f"bias yields for a pTHat-sliced sample. "
            f"Provide sigma_gen_pb in .meta.json or use qcd_mode='auto'."
        )

    if qcd_mode in {"hardbbbar", "hardBc"}:
        if np.isfinite(sigma_gen_pb) and sigma_gen_pb > 0.0:
            return get_parent_sigma_pb(pid, sigma_bbbar_pb=float(sigma_gen_pb))
        raise ValueError(
            f"Hard-QCD mode '{qcd_mode}' for parent {pid} has no valid sigma_gen_pb "
            f"({sigma_gen_pb}). Falling back to inclusive cross-section can strongly "
            f"bias yields for a pTHat-sliced sample. "
            f"Provide sigma_gen_pb in .meta.json or use qcd_mode='auto'."
        )

    return get_parent_sigma_pb(pid)


def expected_signal_events(
    geom_df: pd.DataFrame,
    mass_GeV: float,
    eps2: float,
    benchmark: str,
    lumi_fb: float,
    dirac: bool = False,
    separation_m: float | None = None,
    max_separation_m: float | None = None,
    separation_policy: str = "all-pairs-min",
    decay_seed: int = 12345,
    p_min_GeV: float = 0.6,
    reco_efficiency: float = 1.0,
    decay_cache: DecayCache | None = None,
    separation_pass: np.ndarray | None = None,
    ctau0_m: float | None = None,
    br_per_parent: Dict[int, float] | None = None,
    br_scale: float | None = None,
    timing: dict | None = None,
    geometry_config: GeometryConfig | None = None,
) -> float:
    if timing is not None:
        timing["count_eps2_calls"] = timing.get("count_eps2_calls", 0) + 1
    if ctau0_m is None or br_per_parent is None:
        Ue2, Umu2, Utau2 = couplings_from_eps2(eps2, benchmark)
        with _time_block(timing, "time_model_s"):
            model = HNLModel(mass_GeV=mass_GeV, Ue2=Ue2, Umu2=Umu2, Utau2=Utau2)
        with _time_block(timing, "time_ctau_br_s"):
            ctau0_m = model.ctau0_m
            br_per_parent = model.production_brs()

    if ctau0_m is None or br_per_parent is None:
        raise ValueError("ctau0_m and br_per_parent must be available (either computed or provided).")

    if separation_m is None:
        raise ValueError("Decay-based efficiency is mandatory: provide separation_m.")

    separation_policy_norm = _normalize_separation_policy(separation_policy)
    if max_separation_m is not None and float(max_separation_m) <= float(separation_m):
        raise ValueError(
            f"max_separation_m must be > separation_m (got {max_separation_m} <= {separation_m})."
        )

    geometry_cfg = normalize_geometry_config(geometry_config)

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
    qcd_mode = _normalize_qcd_mode(geom_df["qcd_mode"] if "qcd_mode" in geom_df.columns else None, len(geom_df))
    sigma_gen_pb = _normalize_sigma_gen_pb(
        geom_df["sigma_gen_pb"] if "sigma_gen_pb" in geom_df.columns else None,
        len(geom_df),
    )

    mask_valid = np.isfinite(parent_id)
    if not np.all(mask_valid):
        parent_id = parent_id[mask_valid]
        weights = weights[mask_valid]
        beta_gamma = beta_gamma[mask_valid]
        hits_tube = hits_tube[mask_valid]
        entry = entry[mask_valid]
        length = length[mask_valid]
        qcd_mode = qcd_mode[mask_valid]
        sigma_gen_pb = sigma_gen_pb[mask_valid]

    if len(parent_id) == 0:
        return 0.0

    with _time_block(timing, "time_pdecay_s"):
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
            P_decay[mask_hits] = np.exp(arg_entry) * (-np.expm1(arg_path))

    if separation_pass is None:
        with _time_block(timing, "time_separation_s"):
            selection = DecaySelection(
                separation_m=float(separation_m),
                seed=decay_seed,
                p_min_GeV=p_min_GeV,
                max_separation_m=max_separation_m,
                separation_policy=separation_policy_norm,
            )
            separation_pass = compute_decay_acceptance(
                geom_df=geom_df,
                mass_GeV=mass_GeV,
                flavour=benchmark_to_flavour(benchmark),
                ctau0_m=ctau0_m,
                mesh=build_mesh_once(geometry_cfg),
                selection=selection,
                decay_cache=decay_cache,
            )
    if not np.all(mask_valid):
        separation_pass = separation_pass[mask_valid]
    P_decay = P_decay * separation_pass.astype(float)

    parent_abs = np.abs(parent_id.astype(int))
    tau_parent_id = None
    if "tau_parent_id" in geom_df.columns:
        tau_parent_id = geom_df["tau_parent_id"].to_numpy(dtype=float)
        if not np.all(mask_valid):
            tau_parent_id = tau_parent_id[mask_valid]
        tau_parent_id = np.abs(np.nan_to_num(tau_parent_id, nan=0.0, posinf=0.0, neginf=0.0)).astype(int)

    if tau_parent_id is None:
        mask_from_tau = np.zeros(len(parent_abs), dtype=bool)
    else:
        mask_from_tau = (parent_abs == 15) & (tau_parent_id > 0)

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


    with _time_block(timing, "time_parent_sum_s"):
        for pid in unique_parents:
            BR_parent = br_per_parent.get(int(pid), 0.0)
            if br_scale is not None:
                BR_parent *= br_scale
            if BR_parent <= 0.0:
                missing_br_pdgs.append(int(pid))
                continue

            mask_parent = np.abs(parent_id) == pid
            w_parent = weights[mask_parent]
            P_parent = P_decay[mask_parent]
            qcd_parent = qcd_mode[mask_parent]
            sigma_parent_arr = sigma_gen_pb[mask_parent]

            if len(w_parent) == 0:
                continue

            grouped = _context_groups(qcd_parent, sigma_parent_arr)
            for (ctx_mode, sigma_key), local_idx in grouped.items():
                sigma_val = float(sigma_key) if sigma_key != "nan" else np.nan
                sigma_parent_pb = _resolve_parent_sigma_pb(
                    int(pid),
                    ctx_mode,
                    sigma_val,
                )
                if sigma_parent_pb <= 0.0:
                    missing_xsec_pdgs.append(int(pid))
                    continue

                w = w_parent[local_idx]
                P = P_parent[local_idx]
                w_sum = np.sum(w)
                if w_sum <= 0.0:
                    continue

                eff_parent = np.sum(w * P) / w_sum
                total_expected += lumi_fb * (sigma_parent_pb * 1e3) * BR_parent * eff_parent

    if np.any(mask_from_tau):
        with _time_block(timing, "time_tau_parent_sum_s"):
            assert tau_parent_id is not None
            br_tau_to_N = br_per_parent.get(15, 0.0)
            if br_scale is not None:
                br_tau_to_N *= br_scale
            if br_tau_to_N <= 0.0:
                missing_br_pdgs.append(15)
            else:
                tau_parents = np.unique(tau_parent_id[mask_from_tau])
                for pid in tau_parents:
                    br_parent_to_tau = get_parent_tau_br(int(pid))
                    if br_parent_to_tau <= 0.0:
                        missing_tau_br_pdgs.append(int(pid))
                        continue

                    mask_parent = mask_from_tau & (tau_parent_id == pid)
                    w_parent = weights[mask_parent]
                    P_parent = P_decay[mask_parent]
                    qcd_parent = qcd_mode[mask_parent]
                    sigma_parent_arr = sigma_gen_pb[mask_parent]

                    if len(w_parent) == 0:
                        continue

                    grouped = _context_groups(qcd_parent, sigma_parent_arr)
                    for (ctx_mode, sigma_key), local_idx in grouped.items():
                        sigma_val = float(sigma_key) if sigma_key != "nan" else np.nan
                        sigma_parent_pb = _resolve_parent_sigma_pb(
                            int(pid),
                            ctx_mode,
                            sigma_val,
                        )
                        if sigma_parent_pb <= 0.0:
                            missing_xsec_pdgs.append(int(pid))
                            continue

                        w = w_parent[local_idx]
                        P = P_parent[local_idx]
                        w_sum = np.sum(w)
                        if w_sum <= 0.0:
                            continue

                        eff_parent = np.sum(w * P) / w_sum
                        total_expected += (
                            lumi_fb * (sigma_parent_pb * 1e3) * br_parent_to_tau * br_tau_to_N * eff_parent
                        )

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

    total_expected *= reco_efficiency

    return float(total_expected)


_MESH_CACHE: dict[str, object] = {}


def build_mesh_once(geometry_config: GeometryConfig | None = None):
    cfg = normalize_geometry_config(geometry_config)
    tag = geometry_tag(cfg)
    mesh = _MESH_CACHE.get(tag)
    if mesh is None:
        mesh = build_drainage_gallery_mesh(cfg)
        _MESH_CACHE[tag] = mesh
    return mesh


def benchmark_to_flavour(benchmark: str) -> str:
    if benchmark == "100":
        return "electron"
    if benchmark == "010":
        return "muon"
    if benchmark == "001":
        return "tau"
    raise ValueError(f"Unsupported benchmark: {benchmark} (use '100','010','001').")


