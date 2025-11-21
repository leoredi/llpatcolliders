import os
import glob
import math
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

SPEED_OF_LIGHT = 299792458.0  # m/s

# ---------- MODEL-SPECIFIC PIECES (TO FILL) ----------
# These are implemented using standard Atre/PBC-style approximations:
#   - Total width  Γ_N ∝ G_F^2 m_N^5 |U_ℓ|^2  (neglecting hadronic threshold wiggles)
#   - Production from W:   BR(W→ℓN) / |U_ℓ|^2 = BR(W→ℓν) × ρ_W(m_N)
#   - Production from B:   BR(B→XℓN) / |U_ℓ|^2 ≈ BR(B→Xℓν) × ρ_B(m_N)  (very crude)
#
# For BC6, BC7, BC8 (one non-zero flavour at a time), this captures the correct
# scaling of cτ and production with |U_ℓ|^2 and m_N. If you later want full
# numerical precision, you can replace these with a table or call out to HNLCalc.


def _normalize_flavour(flavour):
    """Map various flavour labels to a canonical key."""
    f = flavour.lower()
    if f in {"e", "ele", "electron"}:
        return "e"
    if f in {"mu", "muon", "m"}:
        return "mu"
    if f in {"tau", "ta"}:
        return "tau"
    # Fallback: just return lowercased string
    return f


def ctau0_m_for_U2_eq_1(mass_GeV, flavour):
    """
    Proper decay length cτ0(m) in meters for |U_ℓ|^2 = 1.

    We use the standard approximate scaling for a Majorana HNL
    that mixes with a *single* active flavour (BC6/7/8 pattern):

        Γ_N(m, |U_ℓ|^2=1) ≈ G_F^2 m^5 / (96 π^3)

    This is the usual Atre/PBC-style expression (up to O(1) factors from
    detailed hadronic thresholds and NC vs CC structure). See e.g.
    Atre et al. (JHEP 05 (2009) 030) and the PBC / Drewes benchmarks.

    We then compute:

        c τ_0(m) = (ħ c) / Γ_N(m, |U_ℓ|^2=1)

    and return it in meters.

    For BC6 (e), BC7 (μ), BC8 (τ) with only one non-zero mixing,
    the lifetime is essentially flavour-independent at this level,
    so we ignore tiny flavour differences here.
    """
    if mass_GeV <= 0:
        raise ValueError("mass_GeV must be positive")

    # Constants in natural units
    G_F = 1.1663787e-5       # GeV^-2
    hbar_GeV_s = 6.582119569e-25  # ħ in GeV·s

    # Leading scaling: Γ ∝ G_F^2 m^5
    gamma_GeV = (G_F ** 2) * (mass_GeV ** 5) / (96.0 * math.pi ** 3)

    # Lifetime and cτ
    tau_s = hbar_GeV_s / gamma_GeV
    ctau_m = SPEED_OF_LIGHT * tau_s
    return ctau_m


def f_prod(mass_GeV, flavour, production_mode):
    """
    Production factor f_prod(m): BR(parent→ℓN) / |U_ℓ|^2 at |U_ℓ|^2 = 1.

    For BC6/BC7/BC8 we assume only one non-zero mixing |U_ℓ|^2, so:
        BR(parent→ℓN) = |U_ℓ|^2 × f_prod(m)

    Implementations:

    1) production_mode in {"W", "WZ"}:
       Use the standard relation
           Γ(W→ℓN) = |U_ℓ|^2 Γ(W→ℓν) ρ_W(x),
           ρ_W(x) = (1 - x)^2 (1 + x/2),  x = m_N^2 / m_W^2
       ⇒ BR(W→ℓN) / |U_ℓ|^2 = BR(W→ℓν) × ρ_W(x).

       We take BR(W→ℓν) ≈ 0.108 for each lepton flavour.

    2) production_mode == "B":
       Very crude effective ansatz for inclusive B→XℓN, based on
       measured inclusive semileptonic BRs and the same phase-space
       factor ρ_B(x) = (1 - x)^2 (1 + x/2), x = m_N^2 / m_B^2:

         BR(B→X e ν) ≈ BR(B→X μ ν) ≈ 0.105
         BR(B→X τ ν) ≈ 0.025

       Then:
         BR(B→XℓN) / |U_ℓ|^2 ≈ BR(B→Xℓν) × ρ_B(x)

       This is *not* a precision model for BC6/7/8; it just gives the
       correct order of magnitude and m_N threshold behaviour. Replace
       with a proper meson-production model if you need exact PBC lines.
    """
    mode = production_mode.upper()
    flav = _normalize_flavour(flavour)

    # --- W production (HL-LHC / ATLAS/CMS-style BC6/7) ---
    if mode in {"W", "WZ"}:
        m_W = 80.379  # GeV
        if mass_GeV >= m_W:
            return 0.0

        BR_W_lnu = 0.108  # per flavour, approx PDG
        x = (mass_GeV / m_W) ** 2
        rho_W = (1.0 - x) ** 2 * (1.0 + 0.5 * x)
        return BR_W_lnu * rho_W

    # --- B production (LHCb-style) ---
    if mode == "B":
        # Use an effective B mass and inclusive semileptonic BRs
        m_B = 5.279  # GeV, representative B meson mass

        if mass_GeV >= m_B:
            return 0.0

        # Inclusive semileptonic branching fractions (approx PDG)
        if flav in {"e", "mu"}:
            BR_B_lnu = 0.105  # B→X e ν or μ ν
        elif flav == "tau":
            BR_B_lnu = 0.025  # B→X τ ν
        else:
            # default to e/μ-like if flavour is weird
            BR_B_lnu = 0.105

        x = (mass_GeV / m_B) ** 2
        rho_B = (1.0 - x) ** 2 * (1.0 + 0.5 * x)

        return BR_B_lnu * rho_B

    # If you introduce more production modes later (Z, Higgs, etc.),
    # add them here.
    raise ValueError(f"Unknown production_mode '{production_mode}'")


# ---------- GENERIC MAPPER: BR_limit(cτ) → U2_limit(m) ----------

def load_br_vs_ctau(csv_file):
    df = pd.read_csv(csv_file)
    # Sometimes ctau is strictly monotonic; we assume it is.
    # We interpolate BR_limit as a function of log10(ctau_m) for stability.
    ctau = df['ctau_m'].values
    br_lim = df['BR_limit'].values

    # Remove infinities / zeros if any
    mask = np.isfinite(br_lim) & (br_lim > 0) & np.isfinite(ctau) & (ctau > 0)
    ctau = ctau[mask]
    br_lim = br_lim[mask]

    # If no valid points remain, return None (no sensitivity)
    if len(ctau) < 2:
        return None

    log_ctau = np.log10(ctau)
    log_br_lim = np.log10(br_lim)

    f = interp1d(
        log_ctau,
        log_br_lim,
        kind='linear',
        bounds_error=False,
        fill_value=np.inf  # outside range → enormous BR_limit
    )
    return f


def find_U2_exclusion(mass_GeV, flavour, production_mode,
                      br_limit_interp,
                      U2_grid=None):
    """
    For a fixed mass, find the excluded |U|^2 band (island).

    Returns (U2_min, U2_max) where:
    - U2_min: too prompt (decays before detector)
    - U2_max: too long-lived (passes through detector)

    The excluded region is U2_min < |U|^2 < U2_max.
    """
    if U2_grid is None:
        # e.g. scan |U|^2 from 1e-15 to 1
        U2_grid = np.logspace(-15, 0, 300)

    ctau0 = ctau0_m_for_U2_eq_1(mass_GeV, flavour)
    f_p = f_prod(mass_GeV, flavour, production_mode)

    # If production factor is zero (e.g. above kinematic threshold), nothing is excluded
    if f_p <= 0.0:
        return None, None

    # For each U2:
    #   ctau(U2) = ctau0 / U2
    #   BR_phys(U2) = U2 * f_p
    #   BR_lim  = BR_limit(ctau(U2))
    # We find where BR_phys crosses BR_lim.
    ratios = []

    for U2 in U2_grid:
        ctau = ctau0 / U2
        log_ctau = np.log10(ctau)
        log_br_lim = br_limit_interp(log_ctau)
        br_lim = 10 ** log_br_lim

        br_phys = U2 * f_p
        ratios.append(br_phys / br_lim)

    ratios = np.array(ratios)

    # Excluded region is where br_phys > br_lim ⇒ ratio > 1.
    mask = ratios > 1.0
    if not np.any(mask):
        return None, None  # not excluded anywhere in this U2 range

    # Find contiguous excluded region
    excluded_indices = np.where(mask)[0]

    # Lower boundary (too prompt)
    U2_min = U2_grid[excluded_indices[0]]

    # Upper boundary (too long-lived)
    U2_max = U2_grid[excluded_indices[-1]]

    return U2_min, U2_max


def main():
    in_pattern = "../output/csv/analysis/HNL_mass_*_BR_vs_ctau.csv"
    rows = []

    for fname in glob.glob(in_pattern):
        base = os.path.basename(fname).replace(".csv", "")
        parts = base.split("_")
        mass = float(parts[2])
        # For low-mass "B-mode", you may have electron/muon/tau; adjust parse if needed.
        flavour = parts[3]

        # Choose production mode by mass:
        #   - below ~5 GeV: B decays (LHCb-style)
        #   - above ~5 GeV: W decays (HL-LHC-style)
        production_mode = "B" if mass < 5.0 else "W"

        print(f"Processing {base} ... (mass={mass} GeV, flavour={flavour}, mode={production_mode})")

        br_interp = load_br_vs_ctau(fname)

        # Skip if no valid data (no sensitivity)
        if br_interp is None:
            print(f"  -> No valid BR limits (no detector acceptance)")
            U2_min, U2_max = None, None
        else:
            U2_min, U2_max = find_U2_exclusion(
                mass, flavour, production_mode,
                br_limit_interp=br_interp
            )

        rows.append({
            "mass_GeV": mass,
            "flavour": flavour,
            "production_mode": production_mode,
            "U2_min": U2_min if U2_min is not None else np.nan,
            "U2_max": U2_max if U2_max is not None else np.nan,
        })

    df_out = pd.DataFrame(rows).sort_values(["flavour", "mass_GeV"])
    out_file = "../output/csv/analysis/HNL_U2_limits_summary.csv"
    df_out.to_csv(out_file, index=False)
    print(f"\nSaved |U|^2 limits to: {out_file}")
    print(df_out)


if __name__ == "__main__":
    main()