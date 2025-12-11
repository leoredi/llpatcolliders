"""
config/production_xsecs.py

Standard LHC Production Cross-Sections for Physics Beyond Colliders (PBC) Analysis.

Reference: CERN-PBC-REPORT-2018-007 / MATHUSLA Physics Case
Energy: 14 TeV (HL-LHC)

These cross-sections represent the TOTAL production of parent mesons at the LHC,
accounting for both the quark-pair production cross-section (ccbar, bbbar) and
the fragmentation fractions (probability that a c/b quark hadronizes into a
specific meson species).

Units: All cross-sections in picobarns (pb)
Conversion: 1 mb = 10^9 pb, 1 μb = 10^6 pb

Per-Parent Counting Methodology
---------------------------------
These cross-sections are used for PER-PARENT counting, not per-event:

    N_sig = Σ_parents [ L × σ_parent × BR(parent→ℓN) × ε_geom(parent) ]

Each parent meson (D, B, K, etc.) represents an independent production channel.
A single pp collision can produce multiple parent mesons → multiple HNLs.
We count each parent's contribution separately because they have different
cross-sections and kinematics.

This matches MATHUSLA/ANUBIS/CODEX-b/AL3X methodology.
"""

# ==========================================
# BASE CROSS-SECTIONS (in picobarns - pb)
# ==========================================

# Sigma(pp -> ccbar) ~ 24 mb at 14 TeV
# This is the QCD production of charm quark pairs
SIGMA_CCBAR_PB = 24.0 * 1e9  # 24 mb = 2.4 × 10^10 pb

# Sigma(pp -> bbbar) ~ 500 μb at 14 TeV
# This is the QCD production of beauty quark pairs
SIGMA_BBBAR_PB = 500.0 * 1e6  # 500 μb = 5.0 × 10^8 pb

# ==========================================
# FRAGMENTATION FRACTIONS
# (Approximate probability quark -> specific meson)
# ==========================================

# Charm Fragmentation (f_c)
# Probability that a c-quark hadronizes into each meson species
FRAG_C_D0     = 0.59  # D0 / D0bar
FRAG_C_DPLUS  = 0.24  # D+ / D-
FRAG_C_DS     = 0.10  # Ds+ / Ds-
FRAG_C_LAMBDA = 0.06  # Λc+ / Λc-

# Beauty Fragmentation (f_b)
# Probability that a b-quark hadronizes into each meson species
FRAG_B_B0     = 0.40  # B0 / B0bar
FRAG_B_BPLUS  = 0.40  # B+ / B-
FRAG_B_BS     = 0.10  # Bs0 / Bs0bar
FRAG_B_LAMBDA = 0.10  # Λb0 / Λb0bar

# ==========================================
# KAON PRODUCTION (Light QCD)
# ==========================================
# Kaon production is dominated by soft QCD processes
# Approximate cross-section for K+/K- production
SIGMA_KAON_PB = 5.0 * 1e10  # ~50 mb (very approximate, soft QCD)

# ==========================================
# ELECTROWEAK PRODUCTION (W/Z BOSONS)
# ==========================================
# W and Z boson production at 14 TeV LHC
# Reference: ATLAS/CMS measurements + NLO calculations
SIGMA_W_PB = 2.0 * 1e8  # σ(pp→W) ~ 200 nb (W+ + W- combined)
SIGMA_Z_PB = 6.0 * 1e7  # σ(pp→Z) ~ 60 nb

# ==========================================
# PARENT PRODUCTION LOOKUP
# Returns: Production Cross Section in pb
# ==========================================

def get_parent_sigma_pb(parent_pdg: int) -> float:
    """
    Returns the total production cross-section for a specific parent meson species.

    This is used to normalize the event yield:
        N_parent = L × σ_parent

    where L is the integrated luminosity in fb^-1.

    Parameters:
    -----------
    parent_pdg : int
        PDG ID of the parent meson (sign doesn't matter - we use abs())

    Returns:
    --------
    float
        Production cross-section in picobarns (pb)

    Notes:
    ------
    - The factor of 2 accounts for particle + antiparticle production
    - Fragmentation fractions are normalized such that Σ_i f_i ≈ 1
    - These are approximate values; real analysis should use NLO QCD calculations
    """
    pid = abs(int(parent_pdg))

    # --- KAON SECTOR (Light QCD) ---
    # NOTE: Only K± (321) is configured in Pythia for HNL production.
    # K_S (310) and K_L (130) are NOT configured, so we don't include them here
    # to avoid inconsistency (having xsec but no events).
    if pid == 321:  # K+ / K-
        return SIGMA_KAON_PB

    # --- CHARM SECTOR ---
    if pid == 421:  # D0 / D0bar
        return SIGMA_CCBAR_PB * FRAG_C_D0 * 2  # *2 for particle+antiparticle
    if pid == 411:  # D+ / D-
        return SIGMA_CCBAR_PB * FRAG_C_DPLUS * 2
    if pid == 431:  # Ds+ / Ds-
        return SIGMA_CCBAR_PB * FRAG_C_DS * 2
    if pid == 4122:  # Λc+ / Λc-
        return SIGMA_CCBAR_PB * FRAG_C_LAMBDA * 2

    # --- BEAUTY SECTOR ---
    if pid == 511:  # B0 / B0bar
        return SIGMA_BBBAR_PB * FRAG_B_B0 * 2
    if pid == 521:  # B+ / B-
        return SIGMA_BBBAR_PB * FRAG_B_BPLUS * 2
    if pid == 531:  # Bs0 / Bs0bar
        return SIGMA_BBBAR_PB * FRAG_B_BS * 2
    if pid == 541:  # Bc+ / Bc- (very rare)
        return SIGMA_BBBAR_PB * 0.001 * 2  # ~0.1% fragmentation
    if pid == 5122:  # Λb0 / Λb0bar
        return SIGMA_BBBAR_PB * FRAG_B_LAMBDA * 2

    # Other beauty baryons (very rough estimates)
    if pid == 5232:  # Σb
        return SIGMA_BBBAR_PB * 0.03 * 2
    if pid == 5332:  # Ωb
        return SIGMA_BBBAR_PB * 0.01 * 2

    # --- TAU (Approximation) ---
    # If Taus are primary parents in your list, they likely come from Ds
    # We approximate σ_tau ~ σ_Ds × BR(Ds→τν)
    # BR(Ds → τ ν) ~ 5.5%
    if pid == 15:
        return (SIGMA_CCBAR_PB * FRAG_C_DS * 2) * 0.055

    # --- ELECTROWEAK BOSONS ---
    if pid == 24:  # W± bosons
        return SIGMA_W_PB
    if pid == 23:  # Z boson
        return SIGMA_Z_PB

    # Default fallback (shouldn't happen if inputs are clean)
    print(f"[WARNING] Unknown parent PDG {pid} in cross-section lookup. Returning 0.")
    return 0.0


def get_sigma_summary():
    """
    Print a summary of production cross-sections for common parents.
    Useful for debugging and verification.
    """
    common_parents = [
        (321, "K+"),
        (421, "D0"),
        (411, "D+"),
        (431, "Ds+"),
        (511, "B0"),
        (521, "B+"),
        (531, "Bs0"),
        (541, "Bc+"),
        (4122, "Λc+"),
        (5122, "Λb0"),
        (15, "τ"),
        (24, "W±"),
        (23, "Z"),
    ]

    print("=" * 70)
    print("LHC Production Cross-Sections (14 TeV, PBC Standard)")
    print("=" * 70)
    print(f"Base rates:")
    print(f"  σ(ccbar) = {SIGMA_CCBAR_PB:.2e} pb = {SIGMA_CCBAR_PB/1e9:.1f} mb")
    print(f"  σ(bbbar) = {SIGMA_BBBAR_PB:.2e} pb = {SIGMA_BBBAR_PB/1e6:.1f} μb")
    print()
    print("Parent meson cross-sections:")
    print(f"{'Parent':<10} {'PDG':>6} {'σ (pb)':>15} {'σ (nb)':>15}")
    print("-" * 70)

    for pdg, name in common_parents:
        sigma_pb = get_parent_sigma_pb(pdg)
        sigma_nb = sigma_pb / 1e3
        print(f"{name:<10} {pdg:>6} {sigma_pb:>15.3e} {sigma_nb:>15.3e}")

    print("=" * 70)


if __name__ == "__main__":
    # Print summary when run as script
    get_sigma_summary()
