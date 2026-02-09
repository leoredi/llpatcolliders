"""
config/production_xsecs.py

Standard LHC Production Cross-Sections for Physics Beyond Colliders (PBC) Analysis.

Primary references (baseline values / methodology):
- CERN-PBC-REPORT-2018-007 (Physics Beyond Colliders LLP sensitivity inputs)
- MATHUSLA physics case: arXiv:1811.00927 (uses the same per-parent counting structure)
- FONLL/LHCb 13/14 TeV reference for SOTA cross-section normalization:
  - σ(ccbar): FONLL NLO+NLL, ~23.6 mb at 14 TeV (Cacciari et al.)
  - σ(bbbar): FONLL/LHCb measurement, ~495 μb at 14 TeV
  - σ(Bc):   CMS/LHCb measurements + BCVEGPY/FONLL, ~0.9 μb at 14 TeV

Notes:
- These are **inclusive, order-of-magnitude** cross-sections suitable for fast
  sensitivity projections. They are not a replacement for a full experimental
  systematics model.
- `SIGMA_W_PB`/`SIGMA_Z_PB` are **inclusive boson production** cross-sections at
  14 TeV (i.e. not multiplied by SM leptonic BRs). Leptonic+HNL branching is
  applied later via `HNLModel.production_brs()`.
- `SIGMA_KAON_PB` is especially uncertain (soft QCD dominated) and is the main
  normalization systematic for the kaon regime.
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

# Sigma(pp -> ccbar) ~ 23.6 mb at 14 TeV
# FONLL NLO+NLL calculation (Cacciari, Greco, Nason)
# Reference: FONLL/LHCb 13/14 TeV reference; consistent with LHCb measurements
# extrapolated to 14 TeV. Previous PBC value was 24 mb (LO-inspired).
SIGMA_CCBAR_PB = 23.6 * 1e9  # 23.6 mb = 2.36 × 10^10 pb

# Sigma(pp -> bbbar) ~ 495 μb at 14 TeV
# FONLL NLO+NLL calculation, validated against LHCb measurements at 7/8/13 TeV
# Reference: FONLL/LHCb 13/14 TeV reference; LHCb-PAPER-2016-031 (13 TeV: 495 ± 52 μb)
# extrapolated to 14 TeV with ~1% increase → ~495 μb central value.
# Previous PBC value was 500 μb (rounded).
SIGMA_BBBAR_PB = 495.0 * 1e6  # 495 μb = 4.95 × 10^8 pb

# Sigma(pp -> Bc) ~ 0.9 μb at 14 TeV (TOTAL Bc+ + Bc- production)
# Bc production is NOT well-modeled by bbbar × fragmentation fraction alone.
# The Bc requires both a b and c quark, making it a special case.
# Reference: FONLL/LHCb 13/14 TeV reference; BCVEGPY (Chang et al.) + CMS/LHCb
# measurements. Literature value ~0.9 μb at LHC energies.
# Note: This is an independent measurement, not derived from SIGMA_BBBAR_PB.
SIGMA_BC_PB = 0.9 * 1e6  # 0.9 μb = 9.0 × 10^5 pb

# ==========================================
# FRAGMENTATION FRACTIONS
# (Approximate probability quark -> specific meson)
# ==========================================

# Charm Fragmentation (f_c)
# Probability that a c-quark hadronizes into each meson species
# These are approximate ground-state fractions; the remaining O(1%) goes to other
# charm hadrons / excited states and is neglected at this projection level.
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
# Kaon production is dominated by soft QCD processes.
# K± cross-section (very approximate, soft QCD dominated):
SIGMA_KAON_PB = 5.0 * 1e10  # ~50 mb for K+ + K-

# K_L cross-section: approximate as ~½ of K± (isospin symmetry).
# This is a stopgap; for precision, measure N(K_L)/N(K±) ratio from Pythia minbias.
# Note: K_S is not included — its contribution is suppressed by τ_S/τ_L ≈ 1/570.
SIGMA_KL_PB = SIGMA_KAON_PB * 0.5  # ~25 mb

# ==========================================
# ELECTROWEAK PRODUCTION (W/Z BOSONS)
# ==========================================
# W and Z boson production at 14 TeV LHC
# Reference: standard inclusive 14 TeV projections used in PBC-style studies
# (e.g. CERN-PBC-REPORT-2018-007; LHC SM cross-section summaries / HXSWG-style inputs).
# These are inclusive production cross-sections (not multiplied by leptonic BRs).
# NLO K-factor for electroweak HNL production (W/Z → ℓN).
# MadGraph LO → NLO correction factor; defined in run_hnl_scan.py:75 as "used in analysis".
# Reference: standard NLO/LO ratio for W/Z + heavy neutral lepton production.
K_FACTOR_EW = 1.3

SIGMA_W_PB = 2.0 * 1e8  # σ(pp→W) ~ 200 nb (W+ + W- combined)
SIGMA_Z_PB = 6.0 * 1e7  # σ(pp→Z) ~ 60 nb

# ==========================================
# PARENT PRODUCTION LOOKUP
# Returns: Production Cross Section in pb
# ==========================================

def get_parent_sigma_pb(
    parent_pdg: int,
    sigma_ccbar_pb: float | None = None,
    sigma_bbbar_pb: float | None = None,
    sigma_bc_pb: float | None = None,
) -> float:
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
    - Base σ(ccbar) and σ(bbbar) use FONLL NLO+NLL values (SOTA)
    - σ(Bc) uses dedicated BCVEGPY/FONLL measurement (not derived from bbbar)
    """
    pid = abs(int(parent_pdg))

    sigma_ccbar = float(SIGMA_CCBAR_PB if sigma_ccbar_pb is None else sigma_ccbar_pb)
    sigma_bbbar = float(SIGMA_BBBAR_PB if sigma_bbbar_pb is None else sigma_bbbar_pb)
    if sigma_bc_pb is None:
        if sigma_bbbar_pb is None:
            sigma_bc = float(SIGMA_BC_PB)
        else:
            # When only a bb̄ slice cross-section is provided, scale Bc by the
            # same Bc/bb̄ ratio as the inclusive defaults.
            sigma_bc = float(sigma_bbbar * (SIGMA_BC_PB / SIGMA_BBBAR_PB))
    else:
        sigma_bc = float(sigma_bc_pb)

    # --- KAON SECTOR (Light QCD) ---
    if pid == 321:  # K+ / K-
        return SIGMA_KAON_PB
    if pid == 130:  # K_L (long-lived neutral kaon)
        return SIGMA_KL_PB
    # Note: K_S (310) is not included — suppressed by τ_S/τ_L ≈ 1/570

    # --- CHARM SECTOR ---
    if pid == 421:  # D0 / D0bar
        return sigma_ccbar * FRAG_C_D0 * 2  # *2 for particle+antiparticle
    if pid == 411:  # D+ / D-
        return sigma_ccbar * FRAG_C_DPLUS * 2
    if pid == 431:  # Ds+ / Ds-
        return sigma_ccbar * FRAG_C_DS * 2
    if pid == 4122:  # Λc+ / Λc-
        return sigma_ccbar * FRAG_C_LAMBDA * 2

    # --- BEAUTY SECTOR ---
    if pid == 511:  # B0 / B0bar
        return sigma_bbbar * FRAG_B_B0 * 2
    if pid == 521:  # B+ / B-
        return sigma_bbbar * FRAG_B_BPLUS * 2
    if pid == 531:  # Bs0 / Bs0bar
        return sigma_bbbar * FRAG_B_BS * 2
    if pid == 541:  # Bc+ / Bc- (FONLL/LHCb 13/14 TeV reference)
        return sigma_bc  # σ(Bc) ~ 0.9 μb (Bc+ + Bc- inclusive)
    if pid == 5122:  # Λb0 / Λb0bar
        return sigma_bbbar * FRAG_B_LAMBDA * 2

    # Other beauty baryons (very rough estimates)
    if pid == 5232:  # Σb
        return sigma_bbbar * 0.03 * 2
    if pid == 5332:  # Ωb
        return sigma_bbbar * 0.01 * 2

    # --- TAU (Approximation) ---
    # If Taus are primary parents in your list, they likely come from Ds
    # We approximate σ_tau ~ σ_Ds × BR(Ds→τν)
    # BR(Ds → τ ν) = 5.48 ± 0.23% (PDG 2024)
    if pid == 15:
        return (sigma_ccbar * FRAG_C_DS * 2) * 0.0548

    # --- ELECTROWEAK BOSONS (with NLO K-factor) ---
    if pid == 24:  # W± bosons
        return SIGMA_W_PB * K_FACTOR_EW
    if pid == 23:  # Z boson
        return SIGMA_Z_PB * K_FACTOR_EW

    # Default fallback (shouldn't happen if inputs are clean)
    print(f"[WARNING] Unknown parent PDG {pid} in cross-section lookup. Returning 0.")
    return 0.0


def get_parent_tau_br(parent_pdg: int) -> float:
    """
    Return SM BR(parent -> tau + nu) for use in fromTau weighting.

    These are the SM branching ratios for meson → τν decays, used to weight
    the tau-decay production chain: parent → τν, τ → N X.

    References:
        - PDG 2024 for Ds → τν
        - World average / R(D*) measurements for B semitauonic
    """
    pid = abs(int(parent_pdg))

    # Ds+ -> tau+ nu_tau (dominant tau source in charm regime)
    # PDG 2024: (5.48 ± 0.23)%
    if pid == 431:
        return 0.0548

    # B0 -> D(*) tau nu (semitauonic, combined D + D*)
    # BR(B0 → D− τ+ ντ) ≈ 0.86% (older B-factory)
    # BR(B0 → D*− τ+ ντ) ≈ 1.40% (world average from R(D*))
    # Combined: ~2.3%
    if pid == 511:
        return 0.023

    # B+ -> D(*) tau nu (semitauonic, combined D0 + D*0)
    # Similar to B0 by isospin
    if pid == 521:
        return 0.023

    # Bs -> Ds(*) tau nu (semitauonic)
    # Less precisely measured, assume similar to B0/B+
    if pid == 531:
        return 0.023

    # Bc+ -> tau+ nu_tau (purely leptonic)
    # BR(Bc→τντ) ≈ 2.4%  (lattice QCD: HPQCD 2020, arXiv:2007.06956;
    # consistent with PDG 2024 indirect constraints).
    # This is an analysis-level approximation — the exact value depends on
    # f_Bc and |V_cb|, with ~10-15% theoretical uncertainty.
    if pid == 541:
        return 0.024

    return 0.0


def get_sigma_summary():
    """
    Print a summary of production cross-sections for common parents.
    Useful for debugging and verification.
    """
    common_parents = [
        (321, "K±"),
        (130, "K_L"),
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
    print(f"Base rates (FONLL/LHCb 13/14 TeV reference):")
    print(f"  σ(ccbar) = {SIGMA_CCBAR_PB:.2e} pb = {SIGMA_CCBAR_PB/1e9:.1f} mb")
    print(f"  σ(bbbar) = {SIGMA_BBBAR_PB:.2e} pb = {SIGMA_BBBAR_PB/1e6:.1f} μb")
    print(f"  σ(Bc)    = {SIGMA_BC_PB:.2e} pb = {SIGMA_BC_PB/1e6:.1f} μb")
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
