"""
limits/u2_limit_calculator.py

End-to-end HNL reinterpretation for the CMS drainage gallery LLP detector.
"""

from __future__ import annotations

import sys
import re
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Tuple, Optional, Dict

import numpy as np
import pandas as pd

# ----------------------------------------------------------------------
# 0. Helpers: repository paths & I/O
# ----------------------------------------------------------------------

THIS_FILE = Path(__file__).resolve()
ANALYSIS_DIR = THIS_FILE.parent           # .../analysis_pbc_test/limits
REPO_ROOT = ANALYSIS_DIR.parents[1]       # .../llpatcolliders
OUTPUT_DIR = REPO_ROOT / "output" / "csv"
SIM_DIR = OUTPUT_DIR / "simulation_new"
GEOM_CACHE_DIR = OUTPUT_DIR / "geometry"
ANALYSIS_OUT_DIR = OUTPUT_DIR / "analysis"

# Ensure the analysis_pbc_test root is on sys.path
ANALYSIS_ROOT = ANALYSIS_DIR.parent
if str(ANALYSIS_ROOT) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_ROOT))

# Import project modules
try:
    from geometry.per_parent_efficiency import (
        build_drainage_gallery_mesh,
        preprocess_hnl_csv,
    )
    from models.hnl_model_hnlcalc import HNLModel
    from config.production_xsecs import get_parent_sigma_pb
except ImportError as e:
    print(f"CRITICAL: Could not import project modules. Check python path.\nError: {e}")
    sys.exit(1)


# ----------------------------------------------------------------------
# 1. Coupling benchmark mapping
# ----------------------------------------------------------------------

def couplings_from_eps2(eps2: float, benchmark: str) -> Tuple[float, float, float]:
    """
    Map PBC benchmark (string) to (Ue^2, Umu^2, Utau^2).
    """
    if benchmark == "100":
        return eps2, 0.0, 0.0
    elif benchmark == "010":
        return 0.0, eps2, 0.0
    elif benchmark == "001":
        return 0.0, 0.0, eps2
    else:
        raise ValueError(f"Unsupported benchmark: {benchmark} (use '100','010','001').")


# ----------------------------------------------------------------------
# 2. Expected signal yield N_sig(m, eps^2)
# ----------------------------------------------------------------------

def expected_signal_events(
    geom_df: pd.DataFrame,
    mass_GeV: float,
    eps2: float,
    benchmark: str,
    lumi_fb: float,
) -> float:
    """
    Compute expected signal events N_sig using per-parent counting.

    Multi-HNL Event Handling
    -------------------------
    Pythia events can contain MULTIPLE HNLs from different parent mesons
    (e.g., one event might produce D0→N, D+→N, B0→N simultaneously).

    We treat each HNL independently because each comes from a distinct
    production process with its own cross-section:

        N_sig = Σ_parents [ L × σ_parent × BR(parent→ℓN) × ε_geom(parent) ]

    This is **per-parent** counting, not per-event counting. We do NOT use
    event-level logic like P_event = 1 - ∏(1 - P_i), because:

    1. Different parents have different cross-sections (σ_D ≠ σ_B ≠ σ_K)
    2. Each parent meson represents an independent production channel
    3. This matches MATHUSLA/ANUBIS/CODEX-b/PBC methodology

    Example: If event #44 produces [B0→N, Bs→N, D+→N, Ds→N], we count
    4 independent signal contributions with their respective σ and BR values.

    Parameters
    ----------
    geom_df : pd.DataFrame
        Geometry dataframe with columns: parent_id, weight, beta_gamma,
        hits_tube, entry_distance, path_length (one row per HNL)
    mass_GeV : float
        HNL mass in GeV
    eps2 : float
        Total coupling squared |U|² = |Ue|² + |Umu|² + |Utau|²
    benchmark : str
        Coupling pattern: "100" (electron), "010" (muon), "001" (tau)
    lumi_fb : float
        Integrated luminosity in fb⁻¹

    Returns
    -------
    float
        Expected number of signal events N_sig
    """
    # 1) Couplings and HNL model
    Ue2, Umu2, Utau2 = couplings_from_eps2(eps2, benchmark)
    model = HNLModel(mass_GeV=mass_GeV, Ue2=Ue2, Umu2=Umu2, Utau2=Utau2)

    # Proper decay length in METRES
    ctau0_m = model.ctau0_m

    # Production BRs per parent
    br_per_parent: Dict[int, float] = model.production_brs()

    # 2) Extract geometry arrays
    df = geom_df  # Reference copy
    
    required_cols = ["parent_id", "weight", "beta_gamma", "hits_tube", "entry_distance", "path_length"]
    missing_cols = [c for c in required_cols if c not in df.columns]
    if missing_cols:
        print(f"[WARN] geom_df missing columns {missing_cols}. Returning N_sig=0.")
        return 0.0

    if len(df) == 0:
        return 0.0

    # Extract raw arrays
    parent_id = df["parent_id"].to_numpy()
    weights = df["weight"].to_numpy(dtype=float)
    beta_gamma = df["beta_gamma"].to_numpy(dtype=float)
    hits_tube = df["hits_tube"].to_numpy(dtype=bool)
    entry = df["entry_distance"].to_numpy(dtype=float)
    length = df["path_length"].to_numpy(dtype=float)

    # --- SAFETY CHECK: DROP BAD PARENT IDs ---
    # This prevents IntCastingNaNError if geometry has NaNs
    mask_valid = np.isfinite(parent_id)
    if not np.all(mask_valid):
        # We silently drop them here to keep the scan fast; 
        # logging happens in the worker function usually.
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
    lam = np.where(lam <= 1e-9, 1e-9, lam) # Prevent divide by zero

    P_decay = np.zeros_like(lam, dtype=float)
    
    mask_hits = hits_tube & (length > 0)
    
    if np.any(mask_hits):
        arg_entry = -entry[mask_hits] / lam[mask_hits]
        arg_path  = -length[mask_hits] / lam[mask_hits]
        
        # Numerically stable: exp(A) * (1 - exp(B)) = exp(A) * (-expm1(B))
        prob_in_tube = np.exp(arg_entry) * (-np.expm1(arg_path))
        P_decay[mask_hits] = prob_in_tube

    # 4) Group by parent species
    # Safe to cast to int now because we filtered non-finites
    unique_parents = np.unique(np.abs(parent_id.astype(int)))
    total_expected = 0.0

    # Track diagnostics for missing PDG codes
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

        # Select geometry rows for this parent
        # We use the array from step 2 (already filtered)
        mask_parent = np.abs(parent_id) == pid
        
        w = weights[mask_parent]
        P = P_decay[mask_parent]
        w_sum = np.sum(w)

        if w_sum <= 0.0:
            continue
            
        eff_parent = np.sum(w * P) / w_sum

        # N = L * sigma * BR * eff (1 pb = 1000 fb)
        total_expected += lumi_fb * (sigma_parent_pb * 1e3) * BR_parent * eff_parent

    # Diagnostic warnings for missing PDG codes (only warn once per mass point)
    # We suppress detailed logging during the 100-point eps2 scan to avoid spam
    # The worker function (_scan_single_mass) should log these once per mass
    if missing_br_pdgs and eps2 == 1e-12:  # Only log at first eps2 point
        n_lost = np.sum(np.isin(parent_id, missing_br_pdgs))
        print(f"[WARN] Mass {mass_GeV:.2f} GeV: {len(missing_br_pdgs)} parent PDG(s) have no HNLCalc BR: {missing_br_pdgs}")
        print(f"       → Discarding {n_lost} events (silent data loss)")

    if missing_xsec_pdgs and eps2 == 1e-12:  # Only log at first eps2 point
        n_lost = np.sum(np.isin(parent_id, missing_xsec_pdgs))
        print(f"[WARN] Mass {mass_GeV:.2f} GeV: {len(missing_xsec_pdgs)} parent PDG(s) have no cross-section: {missing_xsec_pdgs}")
        print(f"       → Discarding {n_lost} events (silent data loss)")

    return float(total_expected)


# ----------------------------------------------------------------------
# 3. Scan eps^2 for a fixed mass
# ----------------------------------------------------------------------

def scan_eps2_for_mass(
    geom_df: pd.DataFrame,
    mass_GeV: float,
    benchmark: str,
    lumi_fb: float,
    N_limit: float = 3.0,
) -> Tuple[np.ndarray, np.ndarray, Optional[float], Optional[float]]:
    
    # Log grid from 1e-12 to 1e-2
    eps2_grid = np.logspace(-12, -2, 100)
    Nsig = np.zeros_like(eps2_grid, dtype=float)

    for i, eps2 in enumerate(eps2_grid):
        Nsig[i] = expected_signal_events(
            geom_df=geom_df,
            mass_GeV=mass_GeV,
            eps2=float(eps2),
            benchmark=benchmark,
            lumi_fb=lumi_fb,
        )

    mask = Nsig >= N_limit
    
    if not np.any(mask):
        return eps2_grid, Nsig, None, None
    
    idx_above = np.where(mask)[0]
    eps2_min = float(eps2_grid[idx_above[0]])
    eps2_max = float(eps2_grid[idx_above[-1]])

    return eps2_grid, Nsig, eps2_min, eps2_max


# ----------------------------------------------------------------------
# 4. Worker for parallel mass scan
# ----------------------------------------------------------------------

def _scan_single_mass(
    mass_val: float,
    mass_str: str,
    flavour: str,
    benchmark: str,
    lumi_fb: float,
    csv_path: Path,
    geom_path: Path,
) -> Optional[dict]:
    
    # Build mesh locally
    try:
        mesh = build_drainage_gallery_mesh()
    except Exception as e:
        print(f"[Worker Error] Mesh build failed: {e}")
        return None
        
    origin = (0.0, 0.0, 0.0)

    # Load or Create Geometry DataFrame
    if geom_path.exists():
        try:
            geom_df = pd.read_csv(geom_path)
            if "entry_distance" not in geom_df.columns:
                print(f"[WARN] Cache corrupt for {mass_str}, rebuilding...")
                raise ValueError("Corrupt geometry cache")
        except Exception:
            geom_df = preprocess_hnl_csv(str(csv_path), mesh, origin=origin)
            geom_df.to_csv(geom_path, index=False)
    else:
        if not csv_path.exists():
            return None
        geom_df = preprocess_hnl_csv(str(csv_path), mesh, origin=origin)
        geom_df.to_csv(geom_path, index=False)

    # --- Data Cleaning ---
    # Check all critical columns for NaN/inf to prevent propagation into physics calculations
    cols_to_check = ["parent_id", "weight", "beta_gamma", "entry_distance", "path_length"]
    mask_valid = np.ones(len(geom_df), dtype=bool)

    for col in cols_to_check:
        if col in geom_df.columns:
            mask_valid &= geom_df[col].notna().to_numpy() & np.isfinite(geom_df[col].to_numpy())
        else:
            print(f"[WARN] m={mass_str} {flavour}: Missing required column '{col}'.")
            return None

    n_total = len(geom_df)
    n_valid = mask_valid.sum()

    if n_valid < n_total:
        n_dropped = n_total - n_valid
        print(f"[INFO] m={mass_str} {flavour}: Dropping {n_dropped} rows with NaNs in geometry/weights.")
        geom_df = geom_df[mask_valid].copy()

    if len(geom_df) == 0:
        print(f"[WARN] m={mass_str} {flavour}: No valid events left after cleaning.")
        return None

    # --- Weight Sanity Check ---
    # Weights should be RELATIVE MC weights (typically 0.1-10), not absolute cross-sections.
    # If weights are suspiciously large (>> 1000), they might be absolute σ values,
    # which would cause double-counting since we normalize to get_parent_sigma_pb().
    if "weight" in geom_df.columns:
        w_mean = geom_df["weight"].mean()
        w_max = geom_df["weight"].max()

        if w_max > 1e6:
            print(f"[ERROR] m={mass_str} {flavour}: Weights suspiciously large (max={w_max:.2e})!")
            print(f"        This looks like absolute cross-section (pb), not relative MC weight.")
            print(f"        Will cause DOUBLE-COUNTING of cross-section. Check CSV generation!")
            print(f"        Expected: weight = pythia.info.weight() (relative, ~1.0)")
            print(f"        NOT: weight = pythia.info.sigmaGen() (absolute σ in pb)")
            return None
        elif w_max > 1000:
            print(f"[WARN] m={mass_str} {flavour}: Weights unusually large (max={w_max:.2e}, mean={w_mean:.2e})")
            print(f"       Verify these are relative MC weights, not absolute cross-sections.")
    # ---------------------

    # Run Physics Scan
    eps2_grid, Nsig, eps2_min, eps2_max = scan_eps2_for_mass(
        geom_df=geom_df,
        mass_GeV=float(mass_val),
        benchmark=benchmark,
        lumi_fb=lumi_fb,
        N_limit=3.0,
    )

    peak_evts = Nsig.max()
    log_msg = f"m={mass_str:<4} {flavour:<4}: Peak Events={peak_evts:.2e}"
    if eps2_min is not None:
        log_msg += f" | Limit: [{eps2_min:.1e}, {eps2_max:.1e}]"
    else:
        log_msg += " | No Sensitivity"
        
    print(log_msg)

    return dict(
        mass_GeV=float(mass_val),
        flavour=flavour,
        benchmark=benchmark,
        eps2_min=eps2_min,
        eps2_max=eps2_max,
        peak_events=peak_evts
    )


# ----------------------------------------------------------------------
# 5. High-level driver
# ----------------------------------------------------------------------

def run_reach_scan(
    flavour: str,
    benchmark: str,
    lumi_fb: float,
    csv_dir: Path = SIM_DIR,
    geom_cache_dir: Path = GEOM_CACHE_DIR,
    results_out: Path = None,
    n_jobs: int = 4,
) -> pd.DataFrame:
    
    if results_out is None:
        results_out = ANALYSIS_OUT_DIR / "HNL_U2_limits_summary.csv"

    geom_cache_dir.mkdir(parents=True, exist_ok=True)
    ANALYSIS_OUT_DIR.mkdir(parents=True, exist_ok=True)

    # 1. Auto-detect files and capture exact mass string
    # New format: HNL_XpYGeV_flavour_regime.csv (e.g., HNL_2p60GeV_muon_beauty.csv)
    # Old format: HNL_mass_X_flavour_Meson.csv (for backwards compatibility)
    # Regimes: kaon, charm, beauty (meson production), ew (electroweak), fromTau/direct (tau)

    # New format patterns
    pattern_new = re.compile(rf"HNL_([0-9]+p[0-9]+)GeV_{flavour}_(kaon|charm|beauty|ew)(?:_direct|_fromTau)?\.csv")
    # Old format patterns
    pattern_old_meson = re.compile(rf"HNL_mass_([0-9p\.]+)_{flavour}_Meson\.csv")
    pattern_old_ew = re.compile(rf"HNL_mass_([0-9p\.]+)_{flavour}_EW\.csv")

    available_files = []
    for f in csv_dir.glob(f"*{flavour}*.csv"):
        # Skip empty files (failed EW simulations are 93B with only CSV header)
        file_size = f.stat().st_size
        if file_size < 1000:  # Less than 1KB = empty file
            print(f"[SKIP] Empty file ({file_size}B): {f.name}")
            continue

        match_new = pattern_new.search(f.name)
        match_old_meson = pattern_old_meson.search(f.name)
        match_old_ew = pattern_old_ew.search(f.name)

        if match_new:
            mass_str = match_new.group(1)   # e.g., "2p60"
            regime = match_new.group(2)     # e.g., "beauty"
            mass_val = float(mass_str.replace('p', '.'))  # Convert 2p60 → 2.60
            # For tau, skip "fromTau" files to avoid double counting (use only direct)
            if "_fromTau" not in f.name:
                available_files.append((mass_val, mass_str, f, regime))
        elif match_old_meson:
            mass_str = match_old_meson.group(1)
            mass_val = float(mass_str.replace('p', '.'))
            available_files.append((mass_val, mass_str, f, "Meson"))
        elif match_old_ew:
            mass_str = match_old_ew.group(1)
            mass_val = float(mass_str.replace('p', '.'))
            available_files.append((mass_val, mass_str, f, "EW"))
    
    # Sort by mass value
    available_files.sort(key=lambda x: x[0])

    if not available_files:
        print(f"[WARN] No Simulation CSVs found for flavour {flavour} in {csv_dir}")
        return pd.DataFrame()

    print(f"--- Starting Scan: {flavour} (Benchmark {benchmark}) ---")
    print(f"Found {len(available_files)} mass points.")

    tasks = []
    for mass_val, mass_str, csv_path, regime in available_files:
        # Use consistent geometry cache naming
        # New format: HNL_2p60GeV_muon_geom.csv (regime-independent)
        # Old format: HNL_mass_2p6_muon_Meson_geom.csv (kept for backwards compatibility)
        if "HNL_mass_" in csv_path.name:
            # Old format
            geom_path = geom_cache_dir / f"HNL_mass_{mass_str}_{flavour}_{regime}_geom.csv"
        else:
            # New format (regime-independent cache)
            geom_path = geom_cache_dir / f"HNL_{mass_str}GeV_{flavour}_geom.csv"
        tasks.append((mass_val, mass_str, flavour, benchmark, lumi_fb, csv_path, geom_path))

    results = []
    
    # 2. Execute
    with ProcessPoolExecutor(max_workers=n_jobs) as executor:
        future_to_mass = {
            executor.submit(_scan_single_mass, *t): t[0] for t in tasks
        }
        for fut in as_completed(future_to_mass):
            try:
                res = fut.result()
                if res:
                    results.append(res)
            except Exception as e:
                print(f"Worker Exception: {e}")

    df_res = pd.DataFrame(results).sort_values("mass_GeV")

    # 3. Merge (Overwrite old results for this benchmark)
    if results_out.exists():
        existing_df = pd.read_csv(results_out)
        mask_keep = ~(
            (existing_df["flavour"] == flavour) & 
            (existing_df["benchmark"] == benchmark)
        )
        existing_df = existing_df[mask_keep]
        final_df = pd.concat([existing_df, df_res], ignore_index=True)
    else:
        final_df = df_res

    final_df.to_csv(results_out, index=False)
    print(f"Saved {len(df_res)} points to {results_out}\n")
    return df_res


# ----------------------------------------------------------------------
# 6. Main execution
# ----------------------------------------------------------------------

if __name__ == "__main__":
    L_HL_LHC_FB = 3000.0
    N_CORES = 4 

    # 1. Electrons (Ve)
    run_reach_scan("electron", "100", L_HL_LHC_FB, n_jobs=N_CORES)

    # 2. Muons (Vmu)
    run_reach_scan("muon", "010", L_HL_LHC_FB, n_jobs=N_CORES)

    # 3. Taus (Vtau)
    run_reach_scan("tau", "001", L_HL_LHC_FB, n_jobs=N_CORES)