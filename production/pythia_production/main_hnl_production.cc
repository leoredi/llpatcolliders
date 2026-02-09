// ==========================================================================
// main_hnl_production.cc
//
// Publication-quality HNL production simulation for far-detector studies.
// Follows methodology of MATHUSLA, ANUBIS, and Physics Beyond Colliders.
//
// Usage: ./main_hnl_production <mass_GeV> <flavor> [nEvents] [mode] [qcdMode] [pTHatMin]
//   flavor: electron, muon, tau (PBC benchmarks BC6/BC7/BC8)
//   mode: 'direct' (default) or 'fromTau' (tau coupling only)
//   qcdMode: 'auto' (default), 'hardBc', 'hardccbar', 'hardbbbar'
//   pTHatMin: override pTHat minimum in GeV (default: mode-dependent)
//
// Production modes (for maximum tau coupling reach):
//   MODE A ("direct"):  B/Ds/W → τ N     (mixing at meson/W vertex)
//   MODE B ("fromTau"): B/Ds/W → τ ν, τ → N X  (mixing at tau decay)
//   → Both modes are O(U_tau²), combine in analysis for maximum sensitivity
//   → Electron and muon use 'direct' mode only
//
// SOTA QCD modes (for transverse detector searches):
//   qcdMode = "auto":      Standard regime-based card selection (default)
//   qcdMode = "hardBc":    Bc production via gg→bb̄/qq̄→bb̄, pTHatMin=15 GeV
//   qcdMode = "hardccbar": Hard cc̄ with pTHatMin cut for high-pT D mesons
//   qcdMode = "hardbbbar": Hard bb̄ with pTHatMin cut for high-pT B mesons
//   → These modes enhance statistics in the kinematic region relevant for
//     transverse detectors (MATHUSLA, CODEX-b)
//
// Output: CSV file with HNL 4-vectors and parent information
//
// ==========================================================================
// CRITICAL: Normalization Strategy
// ==========================================================================
//
// This code uses Pythia as a KINEMATIC GENERATOR ONLY. All physical
// cross-sections and branching ratios are applied externally in Stage 2.
//
// DIVISION OF LABOR:
//
// Stage 1 (This Code - Pythia):
//   → Generates HNL 4-vectors with proper kinematic correlations
//   → Tracks parent species (PDG codes) for each HNL
//   → Records production vertices and boost factors
//   → Internal decay BRs (e.g., "BR=1.0") control RELATIVE sampling
//     of different topologies (2-body vs 3-body, etc.) for kinematics
//   → These internal BRs are NOT physical and do NOT enter final signal
//
// Stage 2 (Analysis Pipeline - HNLCalc + Geometry):
//   → Applies σ(pp → parent) from experimental measurements
//   → Applies BR_inclusive(parent → ℓN) from HNLCalc theory
//   → Computes geometric acceptance ε_geom via ray-tracing
//   → Calculates decay probability P_decay from HNL lifetime
//
// SIGNAL CALCULATION:
//   N_sig = Σ_parents [ L × σ_parent × BR_inclusive × ε_geom × P_decay ]
//
// NO DOUBLE-COUNTING:
//   Even if a parent (e.g., D⁺) has multiple Pythia channels (2-body + 3-body),
//   ALL events from that parent represent the INCLUSIVE parent→ℓN process.
//   The channel mixture approximates inclusive kinematics. Physical
//   normalizations come entirely from HNLCalc, not Pythia.
//
// MAJORANA vs DIRAC:
//   → Simulation generates one Majorana HNL state (N = N̄)
//   → For Dirac interpretation (N ≠ N̄), multiply final yields by factor 2
//   → This factor is NOT included in output CSVs or intermediate results
//   → Apply scaling in final plots with clear labeling
//
// ==========================================================================
//
// References:
//   - arXiv:1805.08567 (HNL phenomenology)
//   - arXiv:1901.09966 (PBC benchmarks)
//   - arXiv:2103.11494 (Pythia validation for HNL)
//   - arXiv:2405.07330 (HNLCalc package)
//
// ==========================================================================

#include "Pythia8/Pythia.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <map>
#include <iomanip>

using namespace Pythia8;

// ==========================================================================
// Physical constants and PDG IDs
// ==========================================================================

// HNL PDG ID: 9900012 (matches MadGraph convention for consistency across production modes)
// Note: Some experiments use 9900015 (SHiP/MATHUSLA), but we use 9900012 to match
// our MadGraph model (SM_HeavyN_CKM_AllMasses_LO) and ensure unified analysis
const int HNL_ID = 9900012;

// Charged mesons that can produce HNL via M -> ℓ N (2-body leptonic)
const std::vector<int> CHARGED_MESONS_2BODY = {
    321,   // K+
    411,   // D+
    431,   // Ds+
    521,   // B+
    541    // Bc+
};

// Neutral mesons/baryons that can produce HNL via semileptonic (3-body)
const std::vector<int> NEUTRAL_MESONS_3BODY = {
    421,   // D0 -> K ℓ N
    511,   // B0 -> D ℓ N  
    531    // Bs -> Ds ℓ N
};

// Baryons that can produce HNL
const std::vector<int> BARYONS_3BODY = {
    4122,  // Λc -> Λ ℓ N or p K ℓ N
    5122   // Λb -> Λc ℓ N
};

// Lepton masses (GeV)
const double M_ELECTRON = 0.000511;
const double M_MUON = 0.10566;
const double M_TAU = 1.777;

// Meson masses (GeV) - for kinematic checks
const std::map<int, double> MESON_MASSES = {
    {130, 0.498},   // K_L (neutral kaon, long-lived)
    {321, 0.494},   // K+
    {411, 1.870},   // D+
    {421, 1.865},   // D0
    {431, 1.968},   // Ds+
    {511, 5.280},   // B0
    {521, 5.279},   // B+
    {531, 5.367},   // Bs
    {541, 6.275}    // Bc+
};

// ==========================================================================
// Helper functions
// ==========================================================================

// Get lepton ID and mass from flavor string
void getLeptonInfo(const std::string& flavor, int& leptonID, int& neutrinoID, 
                   double& leptonMass, std::string& flavorLabel) {
    if (flavor == "electron" || flavor == "e") {
        leptonID = 11;
        neutrinoID = 12;
        leptonMass = M_ELECTRON;
        flavorLabel = "electron";
    } else if (flavor == "muon" || flavor == "mu" || flavor == "μ") {
        leptonID = 13;
        neutrinoID = 14;
        leptonMass = M_MUON;
        flavorLabel = "muon";
    } else if (flavor == "tau" || flavor == "τ") {
        leptonID = 15;
        neutrinoID = 16;
        leptonMass = M_TAU;
        flavorLabel = "tau";
    } else {
        std::cerr << "Unknown flavor: " << flavor << std::endl;
        std::cerr << "Valid options: electron, muon, tau" << std::endl;
        exit(1);
    }
}

// Check if decay is kinematically allowed (2-body: M -> ℓ N)
bool isKinematicallyAllowed2Body(double mParent, double mLepton, double mHNL) {
    return (mHNL < mParent - mLepton);
}

// Determine production regime based on HNL mass and flavor
// For tau coupling: kaons cannot produce taus (m_K < m_tau), so use charm/beauty
std::string getProductionRegime(double mHNL, const std::string& flavor = "") {
    // Tau coupling requires heavy meson parents (Ds, B) - kaons cannot produce taus
    if (flavor == "tau") {
        if (mHNL < 2.0) return "charm";   // Ds-dominated for tau
        return "beauty";                   // B-dominated for tau
    }
    // Electron/muon coupling: standard mass-based regime
    if (mHNL < 0.5) return "kaon";        // Kaon-dominated regime
    if (mHNL < 2.0) return "charm";       // Charm-dominated regime
    return "beauty";                       // Beauty regime (2.0-10.0 GeV)
}

// Convert mass to filename-safe label
std::string massToLabel(double mass) {
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(2) << mass;
    std::string label = ss.str();
    for (char& c : label) {
        if (c == '.') c = 'p';
    }
    return label;
}

// Find physical parent (skip copies, find original decaying particle)
int findPhysicalParent(const Event& event, int iParticle, int forbiddenId = 0) {
    if (iParticle < 0 || iParticle >= event.size()) return 0;
    
    // Get the top copy of the HNL (earliest in decay chain)
    int iTop = event[iParticle].iTopCopy();
    if (iTop < 0 || iTop >= event.size()) iTop = iParticle;
    
    // Get mother of the top copy
    int iMother = event[iTop].mother1();
    if (iMother <= 0 || iMother >= event.size()) return 0;
    
    // Optionally guard against pathological self-parenting (HNL -> HNL)
    if (forbiddenId != 0 && std::abs(event[iMother].id()) == std::abs(forbiddenId)) {
        return 0;
    }
    
    return event[iMother].id();
}

// ==========================================================================
// Configure forced decays for meson production
// ==========================================================================
//
// We implement:
//   - 2-body leptonic: M+ → ℓ+ N (K+, D+, Ds+, B+, Bc+)
//   - 3-body semileptonic: M → M' ℓ N (representative channels)
//
// For semileptonic decays, we use ONE representative exclusive channel
// per parent meson (e.g., D0 → K ℓ N, B0 → D ℓ N) with phase-space
// kinematics (meMode=0). This is validated by arXiv:2103.11494 as
// adequate for sensitivity estimates.
//
// IMPORTANT: The branching ratios here are artificially set to ~100%.
// Actual inclusive BRs must be applied as weights from external
// calculations (e.g., HNLCalc using formulas from arXiv:1805.08567).
// ==========================================================================

void configureMesonDecays(Pythia& pythia, int leptonID,
                          double mHNL, double mLepton, bool verbose = true) {

    std::string hnl = std::to_string(HNL_ID);
    std::string lep = std::to_string(leptonID);
    std::string lepBar = std::to_string(-leptonID);

    int nChannelsConfigured = 0;

    if (verbose) {
        std::cout << "\n=== Configuring HNL decay channels ===" << std::endl;
        std::cout << "HNL mass: " << mHNL << " GeV" << std::endl;
        std::cout << "Lepton ID: " << leptonID << " (mass " << mLepton << " GeV)" << std::endl;
    }

    // -----------------------------------------------------------------------
    // 2-body leptonic decays: M+ -> ℓ+ N
    // -----------------------------------------------------------------------
    
    // K+ -> ℓ+ N (and K- -> ℓ- Nbar)
    // NOTE: Kaons have mayDecay=off by default in Pythia (long-lived for detector sim).
    // We must explicitly enable decays for HNL production.
    if (isKinematicallyAllowed2Body(MESON_MASSES.at(321), mLepton, mHNL)) {
        pythia.readString("321:mayDecay = on");
        pythia.readString("-321:mayDecay = on");
        pythia.readString("321:onMode = off");
        pythia.readString("321:addChannel = 1 1.0 0 " + lepBar + " " + hnl);
        pythia.readString("-321:onMode = off");
        pythia.readString("-321:addChannel = 1 1.0 0 " + lep + " " + hnl);
        if (verbose) std::cout << "  K± -> ℓ N : ENABLED (mayDecay forced on)" << std::endl;
        nChannelsConfigured++;
    } else if (verbose) {
        std::cout << "  K± -> ℓ N : DISABLED (kinematically forbidden)" << std::endl;
    }
    
    // D+ -> ℓ+ N
    if (isKinematicallyAllowed2Body(MESON_MASSES.at(411), mLepton, mHNL)) {
        pythia.readString("411:onMode = off");
        pythia.readString("411:addChannel = 1 1.0 0 " + lepBar + " " + hnl);
        pythia.readString("-411:onMode = off");
        pythia.readString("-411:addChannel = 1 1.0 0 " + lep + " " + hnl);
        if (verbose) std::cout << "  D± -> ℓ N : ENABLED" << std::endl;
        nChannelsConfigured++;
    } else if (verbose) {
        std::cout << "  D± -> ℓ N : DISABLED (kinematically forbidden)" << std::endl;
    }
    
    // Ds+ -> ℓ+ N
    if (isKinematicallyAllowed2Body(MESON_MASSES.at(431), mLepton, mHNL)) {
        pythia.readString("431:onMode = off");
        pythia.readString("431:addChannel = 1 1.0 0 " + lepBar + " " + hnl);
        pythia.readString("-431:onMode = off");
        pythia.readString("-431:addChannel = 1 1.0 0 " + lep + " " + hnl);
        if (verbose) std::cout << "  Ds± -> ℓ N : ENABLED" << std::endl;
        nChannelsConfigured++;
    } else if (verbose) {
        std::cout << "  Ds± -> ℓ N : DISABLED (kinematically forbidden)" << std::endl;
    }
    
    // B+ -> ℓ+ N
    if (isKinematicallyAllowed2Body(MESON_MASSES.at(521), mLepton, mHNL)) {
        pythia.readString("521:onMode = off");
        pythia.readString("521:addChannel = 1 1.0 0 " + lepBar + " " + hnl);
        pythia.readString("-521:onMode = off");
        pythia.readString("-521:addChannel = 1 1.0 0 " + lep + " " + hnl);
        if (verbose) std::cout << "  B± -> ℓ N : ENABLED" << std::endl;
        nChannelsConfigured++;
    } else if (verbose) {
        std::cout << "  B± -> ℓ N : DISABLED (kinematically forbidden)" << std::endl;
    }
    
    // Bc+ -> ℓ+ N (rare, but physically correct)
    if (isKinematicallyAllowed2Body(MESON_MASSES.at(541), mLepton, mHNL)) {
        pythia.readString("541:onMode = off");
        pythia.readString("541:addChannel = 1 1.0 0 " + lepBar + " " + hnl);
        pythia.readString("-541:onMode = off");
        pythia.readString("-541:addChannel = 1 1.0 0 " + lep + " " + hnl);
        if (verbose) std::cout << "  Bc± -> ℓ N : ENABLED" << std::endl;
        nChannelsConfigured++;
    } else if (verbose) {
        std::cout << "  Bc± -> ℓ N : DISABLED (kinematically forbidden)" << std::endl;
    }
    
    // -----------------------------------------------------------------------
    // 3-body semileptonic decays: M -> M' ℓ N
    // Note: Using meMode=0 (phase space) for simplicity
    //       For proper matrix elements, use external decay tools
    // -----------------------------------------------------------------------

    // K_L -> π± ℓ∓ N (semileptonic, 3-body)
    // K_L is self-conjugate: both π⁺ℓ⁻ and π⁻ℓ⁺ are allowed with equal weight.
    // Note: K_S is omitted — its contribution is suppressed by τ_S/τ_L ≈ 1/570
    // relative to K_L (HNLCalc handles this via lifetime in BR calculation).
    // NOTE: K_L has mayDecay=off by default in Pythia (long-lived for detector sim).
    double mKL = MESON_MASSES.at(130);
    double mPiCharged = 0.140;  // π± mass
    if (mHNL + mLepton + mPiCharged < mKL) {
        pythia.readString("130:mayDecay = on");
        pythia.readString("130:onMode = off");
        // K_L → π⁻ ℓ⁺ N
        pythia.readString("130:addChannel = 1 0.5 0 -211 " + lepBar + " " + hnl);
        // K_L → π⁺ ℓ⁻ N
        pythia.readString("130:addChannel = 1 0.5 0 211 " + lep + " " + hnl);
        if (verbose) std::cout << "  K_L -> π ℓ N : ENABLED (3-body, mayDecay forced on)" << std::endl;
        nChannelsConfigured++;
    } else if (verbose) {
        std::cout << "  K_L -> π ℓ N : DISABLED (kinematically forbidden)" << std::endl;
    }

    // D0 -> K- ℓ+ N (semileptonic)
    double mD0 = MESON_MASSES.at(421);
    double mK = MESON_MASSES.at(321);
    if (mHNL + mLepton + mK < mD0) {
        pythia.readString("421:onMode = off");
        pythia.readString("421:addChannel = 1 1.0 0 -321 " + lepBar + " " + hnl);
        pythia.readString("-421:onMode = off");
        pythia.readString("-421:addChannel = 1 1.0 0 321 " + lep + " " + hnl);
        if (verbose) std::cout << "  D0 -> K ℓ N : ENABLED (3-body)" << std::endl;
        nChannelsConfigured++;
    } else if (verbose) {
        std::cout << "  D0 -> K ℓ N : DISABLED (kinematically forbidden)" << std::endl;
    }

    // D+ -> K0bar ℓ+ N (semileptonic, K0bar = -311)
    double mDplus = MESON_MASSES.at(411);
    double mK0 = 0.498;  // K0 mass
    if (mHNL + mLepton + mK0 < mDplus) {
        // Note: D+ 2-body leptonic D+ → ℓ N is already added above
        // This adds the semileptonic channel
        pythia.readString("411:addChannel = 1 0.5 0 -311 " + lepBar + " " + hnl);
        pythia.readString("-411:addChannel = 1 0.5 0 311 " + lep + " " + hnl);
        if (verbose) std::cout << "  D± -> K0 ℓ N : ENABLED (3-body)" << std::endl;
        nChannelsConfigured++;
    } else if (verbose) {
        std::cout << "  D± -> K0 ℓ N : DISABLED (kinematically forbidden)" << std::endl;
    }

    // B0 -> D- ℓ+ N (semileptonic)
    double mB0 = MESON_MASSES.at(511);
    double mDminus = MESON_MASSES.at(411);
    if (mHNL + mLepton + mDminus < mB0) {
        pythia.readString("511:onMode = off");
        pythia.readString("511:addChannel = 1 1.0 0 -411 " + lepBar + " " + hnl);
        pythia.readString("-511:onMode = off");
        pythia.readString("-511:addChannel = 1 1.0 0 411 " + lep + " " + hnl);
        if (verbose) std::cout << "  B0 -> D ℓ N : ENABLED (3-body)" << std::endl;
        nChannelsConfigured++;
    } else if (verbose) {
        std::cout << "  B0 -> D ℓ N : DISABLED (kinematically forbidden)" << std::endl;
    }

    // B+ -> D0bar ℓ+ N (semileptonic, D0bar = -421)
    double mBplus = MESON_MASSES.at(521);
    double mD0mass = MESON_MASSES.at(421);
    if (mHNL + mLepton + mD0mass < mBplus) {
        // Note: B+ 2-body leptonic B+ → ℓ N is already added above
        // This adds the semileptonic channel
        pythia.readString("521:addChannel = 1 0.5 0 -421 " + lepBar + " " + hnl);
        pythia.readString("-521:addChannel = 1 0.5 0 421 " + lep + " " + hnl);
        if (verbose) std::cout << "  B± -> D0 ℓ N : ENABLED (3-body)" << std::endl;
        nChannelsConfigured++;
    } else if (verbose) {
        std::cout << "  B± -> D0 ℓ N : DISABLED (kinematically forbidden)" << std::endl;
    }

    // Bs -> Ds- ℓ+ N (semileptonic)
    double mBs = MESON_MASSES.at(531);
    double mDs = MESON_MASSES.at(431);
    if (mHNL + mLepton + mDs < mBs) {
        pythia.readString("531:onMode = off");
        pythia.readString("531:addChannel = 1 1.0 0 -431 " + lepBar + " " + hnl);
        pythia.readString("-531:onMode = off");
        pythia.readString("-531:addChannel = 1 1.0 0 431 " + lep + " " + hnl);
        if (verbose) std::cout << "  Bs -> Ds ℓ N : ENABLED (3-body)" << std::endl;
        nChannelsConfigured++;
    } else if (verbose) {
        std::cout << "  Bs -> Ds ℓ N : DISABLED (kinematically forbidden)" << std::endl;
    }
    
    // Λb -> Λc ℓ- N (baryon semileptonic)
    double mLb = 5.620;  // Lambda_b mass
    double mLc = 2.286;  // Lambda_c mass
    if (mHNL + mLepton + mLc < mLb) {
        pythia.readString("5122:onMode = off");
        pythia.readString("5122:addChannel = 1 1.0 0 4122 " + lep + " " + hnl);
        pythia.readString("-5122:onMode = off");
        pythia.readString("-5122:addChannel = 1 1.0 0 -4122 " + lepBar + " " + hnl);
        if (verbose) std::cout << "  Λb -> Λc ℓ N : ENABLED (3-body)" << std::endl;
        nChannelsConfigured++;
    } else if (verbose) {
        std::cout << "  Λb -> Λc ℓ N : DISABLED (kinematically forbidden)" << std::endl;
    }

    // Λc -> Λ ℓ+ N (baryon semileptonic)
    double mLambda0 = 1.115;  // Lambda^0 mass
    if (mHNL + mLepton + mLambda0 < mLc) {
        pythia.readString("4122:onMode = off");
        pythia.readString("4122:addChannel = 1 1.0 0 3122 " + lep + " " + hnl);
        pythia.readString("-4122:onMode = off");
        pythia.readString("-4122:addChannel = 1 1.0 0 -3122 " + lepBar + " " + hnl);
        if (verbose) std::cout << "  Λc -> Λ ℓ N : ENABLED (3-body)" << std::endl;
        nChannelsConfigured++;
    } else if (verbose) {
        std::cout << "  Λc -> Λ ℓ N : DISABLED (kinematically forbidden)" << std::endl;
    }
    
    if (verbose) {
        std::cout << "Total channels configured: " << nChannelsConfigured << std::endl;
        std::cout << "==========================================\n" << std::endl;
    }
}

// ==========================================================================
// Configure meson decays to τν for "fromTau" production mode
// ==========================================================================
//
// Forces parent mesons to decay to τν, ensuring every event produces a tau.
// This avoids wasting CPU on events where mesons decay to other channels.
//
// SM branching fractions (for reference - NOT used here, applied via HNLCalc):
//   Ds → τν:       ~5.3%
//   B → D τν:      ~0.9%   (world average)
//   B → D* τν:     ~1.4%   (world average from R(D*) measurements)
//
// By forcing 100% here, we get ~20-50x speedup. Physical BRs are applied
// externally in the analysis pipeline (consistent with direct mode methodology).
//
// The tau_parent_id column in output identifies the grandfather meson for
// correct BR weighting in HNLCalc.
// ==========================================================================

void configureMesonDecaysToTauNu(Pythia& pythia, double mHNL, bool verbose = true) {

    if (verbose) {
        std::cout << "\n=== Configuring meson → τν decays (for fromTau mode) ===" << std::endl;
        std::cout << "HNL mass: " << mHNL << " GeV" << std::endl;
    }

    int nChannelsConfigured = 0;

    // -----------------------------------------------------------------------
    // Ds± → τ± ντ (dominant tau source in charm regime)
    // SM BR ~ 5.3%
    // -----------------------------------------------------------------------
    double mDs = MESON_MASSES.at(431);
    if (M_TAU < mDs) {
        pythia.readString("431:onMode = off");
        pythia.readString("431:addChannel = 1 1.0 0 -15 16");   // Ds+ → τ+ ντ
        pythia.readString("-431:onMode = off");
        pythia.readString("-431:addChannel = 1 1.0 0 15 -16");  // Ds- → τ- ν̄τ
        if (verbose) std::cout << "  Ds± → τ ν : ENABLED" << std::endl;
        nChannelsConfigured++;
    }

    // -----------------------------------------------------------------------
    // B+ → D̄0 τ+ ντ and B+ → D̄*0 τ+ ντ (semileptonic)
    // SM BR(B→Dτν) ~ 0.9%, BR(B→D*τν) ~ 1.4%
    // Weight ratio ~0.4:0.6 to approximate relative BRs
    // -----------------------------------------------------------------------
    double mBplus = MESON_MASSES.at(521);
    double mD0 = MESON_MASSES.at(421);
    double mDstar0 = 2.007;  // D*0 mass

    pythia.readString("521:onMode = off");
    pythia.readString("-521:onMode = off");

    bool bplus_d_ok = (M_TAU + mD0 < mBplus);
    bool bplus_dstar_ok = (M_TAU + mDstar0 < mBplus);

    if (bplus_d_ok && bplus_dstar_ok) {
        // Both channels open - weight by approximate BR ratio
        pythia.readString("521:addChannel = 1 0.4 0 -421 -15 16");   // B+ → D̄0 τ+ ντ
        pythia.readString("521:addChannel = 1 0.6 0 -423 -15 16");   // B+ → D̄*0 τ+ ντ
        pythia.readString("-521:addChannel = 1 0.4 0 421 15 -16");
        pythia.readString("-521:addChannel = 1 0.6 0 423 15 -16");
        if (verbose) std::cout << "  B± → D(*)0 τ ν : ENABLED (D:D* = 0.4:0.6)" << std::endl;
        nChannelsConfigured++;
    } else if (bplus_d_ok) {
        pythia.readString("521:addChannel = 1 1.0 0 -421 -15 16");
        pythia.readString("-521:addChannel = 1 1.0 0 421 15 -16");
        if (verbose) std::cout << "  B± → D0 τ ν : ENABLED (D* closed)" << std::endl;
        nChannelsConfigured++;
    } else if (verbose) {
        std::cout << "  B± → D(*)0 τ ν : DISABLED (kinematically forbidden)" << std::endl;
    }

    // -----------------------------------------------------------------------
    // B0 → D− τ+ ντ and B0 → D*− τ+ ντ (semileptonic)
    // -----------------------------------------------------------------------
    double mB0 = MESON_MASSES.at(511);
    double mDminus = MESON_MASSES.at(411);
    double mDstarMinus = 2.010;  // D*- mass

    pythia.readString("511:onMode = off");
    pythia.readString("-511:onMode = off");

    bool b0_d_ok = (M_TAU + mDminus < mB0);
    bool b0_dstar_ok = (M_TAU + mDstarMinus < mB0);

    if (b0_d_ok && b0_dstar_ok) {
        pythia.readString("511:addChannel = 1 0.4 0 -411 -15 16");   // B0 → D− τ+ ντ
        pythia.readString("511:addChannel = 1 0.6 0 -413 -15 16");   // B0 → D*− τ+ ντ
        pythia.readString("-511:addChannel = 1 0.4 0 411 15 -16");
        pythia.readString("-511:addChannel = 1 0.6 0 413 15 -16");
        if (verbose) std::cout << "  B0 → D(*)± τ ν : ENABLED (D:D* = 0.4:0.6)" << std::endl;
        nChannelsConfigured++;
    } else if (b0_d_ok) {
        pythia.readString("511:addChannel = 1 1.0 0 -411 -15 16");
        pythia.readString("-511:addChannel = 1 1.0 0 411 15 -16");
        if (verbose) std::cout << "  B0 → D± τ ν : ENABLED (D* closed)" << std::endl;
        nChannelsConfigured++;
    } else if (verbose) {
        std::cout << "  B0 → D(*)± τ ν : DISABLED (kinematically forbidden)" << std::endl;
    }

    // -----------------------------------------------------------------------
    // Bs → Ds− τ+ ντ and Bs → Ds*− τ+ ντ (semileptonic)
    // -----------------------------------------------------------------------
    double mBs = MESON_MASSES.at(531);
    double mDsPlus = MESON_MASSES.at(431);
    double mDsstar = 2.112;  // Ds*- mass

    pythia.readString("531:onMode = off");
    pythia.readString("-531:onMode = off");

    bool bs_ds_ok = (M_TAU + mDsPlus < mBs);
    bool bs_dsstar_ok = (M_TAU + mDsstar < mBs);

    if (bs_ds_ok && bs_dsstar_ok) {
        pythia.readString("531:addChannel = 1 0.4 0 -431 -15 16");   // Bs → Ds− τ+ ντ
        pythia.readString("531:addChannel = 1 0.6 0 -433 -15 16");   // Bs → Ds*− τ+ ντ
        pythia.readString("-531:addChannel = 1 0.4 0 431 15 -16");
        pythia.readString("-531:addChannel = 1 0.6 0 433 15 -16");
        if (verbose) std::cout << "  Bs → Ds(*) τ ν : ENABLED (Ds:Ds* = 0.4:0.6)" << std::endl;
        nChannelsConfigured++;
    } else if (bs_ds_ok) {
        pythia.readString("531:addChannel = 1 1.0 0 -431 -15 16");
        pythia.readString("-531:addChannel = 1 1.0 0 431 15 -16");
        if (verbose) std::cout << "  Bs → Ds τ ν : ENABLED (Ds* closed)" << std::endl;
        nChannelsConfigured++;
    } else if (verbose) {
        std::cout << "  Bs → Ds(*) τ ν : DISABLED (kinematically forbidden)" << std::endl;
    }

    // -----------------------------------------------------------------------
    // Bc+ → τ+ ντ (purely leptonic)
    // SM BR(Bc→τντ) ≈ 2.4% (lattice QCD: HPQCD 2020, arXiv:2007.06956)
    // m(Bc) = 6.275 GeV >> m(τ) = 1.777 GeV — always kinematically open.
    // Without forcing this channel, hardBc + fromTau produces near-zero
    // statistics because Bc→τν is rare in SM branching table.
    // -----------------------------------------------------------------------
    double mBc = MESON_MASSES.at(541);
    if (M_TAU < mBc) {
        pythia.readString("541:onMode = off");
        pythia.readString("541:addChannel = 1 1.0 0 -15 16");   // Bc+ → τ+ ντ
        pythia.readString("-541:onMode = off");
        pythia.readString("-541:addChannel = 1 1.0 0 15 -16");  // Bc- → τ- ν̄τ
        if (verbose) std::cout << "  Bc± → τ ν : ENABLED" << std::endl;
        nChannelsConfigured++;
    } else if (verbose) {
        std::cout << "  Bc± → τ ν : DISABLED (kinematically forbidden)" << std::endl;
    }

    if (verbose) {
        std::cout << "Total τν channels configured: " << nChannelsConfigured << std::endl;
        std::cout << "==========================================\n" << std::endl;
    }
}

// ==========================================================================
// Configure tau decays for "fromTau" production mode
// ==========================================================================
//
// PHYSICS: For BC8 (tau coupling), there are TWO independent O(U_tau²) sources:
//
//   MODE A ("direct"):  B/Ds/W → τ N  (mixing at meson/W vertex)
//   MODE B ("fromTau"): B/Ds/W → τ ν → N X  (mixing at tau decay)
//
// To avoid O(U⁴) contamination, we generate these as SEPARATE samples:
//   - "direct" mode: Mesons/W forced to τN, taus decay SM
//   - "fromTau" mode: Mesons forced to τν, taus forced to NX
//
// The two samples are combined in the analysis pipeline.
//
// This function configures the tau → NX part of MODE B.
// ==========================================================================

void configureTauDecays(Pythia& pythia, double mHNL, bool verbose = true) {
    std::string hnl = std::to_string(HNL_ID);
    double mTau = M_TAU;

    if (verbose) {
        std::cout << "\n=== Configuring tau → N X decays (MODE B: fromTau) ===" << std::endl;
        std::cout << "HNL mass: " << mHNL << " GeV" << std::endl;
    }

    // Turn off all SM tau decays
    pythia.readString("15:onMode = off");
    pythia.readString("-15:onMode = off");

    // Kinematic thresholds for each channel
    double mPi = 0.140;   // charged pion
    double mRho = 0.775;  // rho mass
    double m3Pi = 3 * mPi;
    double mMu = M_MUON;  // 0.106 GeV
    double mE = M_ELECTRON;  // 0.0005 GeV

    // 2-body hadronic channels
    const bool allow_pi = (mHNL + mPi < mTau);      // < 1.637 GeV
    const bool allow_rho = (mHNL + mRho < mTau);    // < 1.002 GeV
    const bool allow_tripi = (mHNL + m3Pi < mTau);  // < 1.357 GeV

    // 3-body leptonic channels: τ → ℓ ν N (effective limit is m_N < m_τ - m_ℓ)
    const bool allow_mu = (mHNL + mMu < mTau);      // < 1.671 GeV
    const bool allow_e = (mHNL + mE < mTau);        // < 1.777 GeV

    // Representative kinematics mixture weights (NOT physical BRs).
    // Physical τ→NX branching is applied later by HNLCalc in the analysis.
    // Weights prioritize: hadronic when available, leptonic near endpoint.
    double pi_weight = 0.0;
    double rho_weight = 0.0;
    double tripi_weight = 0.0;
    double mu_weight = 0.0;
    double e_weight = 0.0;

    if (allow_pi) {
        // Hadronic channels available - use them primarily
        if (allow_rho && allow_tripi) {
            rho_weight = 0.45;
            tripi_weight = 0.25;
            pi_weight = 0.20;
            mu_weight = 0.05;
            e_weight = 0.05;
        } else if (allow_rho && !allow_tripi) {
            rho_weight = 0.50;
            pi_weight = 0.35;
            mu_weight = 0.08;
            e_weight = 0.07;
        } else if (!allow_rho && allow_tripi) {
            tripi_weight = 0.25;
            pi_weight = 0.55;
            mu_weight = 0.10;
            e_weight = 0.10;
        } else {
            // Only π and leptonic
            pi_weight = 0.70;
            mu_weight = 0.15;
            e_weight = 0.15;
        }
    } else if (allow_mu) {
        // Hadronic closed, but μ channel still open (1.637 < m_N < 1.671 GeV)
        mu_weight = 0.50;
        e_weight = 0.50;
    } else if (allow_e) {
        // Only e channel open (1.671 < m_N < 1.777 GeV)
        e_weight = 1.00;
    }

    if (verbose) {
        std::cout << "  τ→NX channels: "
                  << "π=" << pi_weight
                  << ", ρ=" << rho_weight
                  << ", 3π=" << tripi_weight
                  << ", μν=" << mu_weight
                  << ", eν=" << e_weight
                  << std::endl;
    }

    int nChannels = 0;

    if (rho_weight > 0.0) {
        pythia.readString("15:addChannel = 1 " + std::to_string(rho_weight) + " 0 -213 " + hnl);
        pythia.readString("-15:addChannel = 1 " + std::to_string(rho_weight) + " 0 213 " + hnl);
        if (verbose) std::cout << "  τ → ρ N : ENABLED" << std::endl;
        nChannels++;
    }

    if (tripi_weight > 0.0) {
        pythia.readString("15:addChannel = 1 " + std::to_string(tripi_weight) + " 0 -211 -211 211 " + hnl);
        pythia.readString("-15:addChannel = 1 " + std::to_string(tripi_weight) + " 0 211 211 -211 " + hnl);
        if (verbose) std::cout << "  τ → 3π N : ENABLED" << std::endl;
        nChannels++;
    }

    if (pi_weight > 0.0) {
        pythia.readString("15:addChannel = 1 " + std::to_string(pi_weight) + " 0 -211 " + hnl);
        pythia.readString("-15:addChannel = 1 " + std::to_string(pi_weight) + " 0 211 " + hnl);
        if (verbose) std::cout << "  τ → π N : ENABLED" << std::endl;
        nChannels++;
    }

    // Leptonic channels: τ- → ℓ- ν̄_ℓ N (3-body, use phase space meMode=0)
    // PDG codes: e=11, νe=12, μ=13, νμ=14
    if (mu_weight > 0.0) {
        pythia.readString("15:addChannel = 1 " + std::to_string(mu_weight) + " 0 13 -14 " + hnl);
        pythia.readString("-15:addChannel = 1 " + std::to_string(mu_weight) + " 0 -13 14 " + hnl);
        if (verbose) std::cout << "  τ → μ ν̄ N : ENABLED" << std::endl;
        nChannels++;
    }

    if (e_weight > 0.0) {
        pythia.readString("15:addChannel = 1 " + std::to_string(e_weight) + " 0 11 -12 " + hnl);
        pythia.readString("-15:addChannel = 1 " + std::to_string(e_weight) + " 0 -11 12 " + hnl);
        if (verbose) std::cout << "  τ → e ν̄ N : ENABLED" << std::endl;
        nChannels++;
    }

    if (nChannels == 0 && verbose) {
        std::cout << "  WARNING: No tau decay channels available at this mass!" << std::endl;
    }

    // Note: In reality, τ → N + X has many channels (π, ρ, ℓνν, etc.)
    // We use representative modes (π and, when allowed, ρ) for geometric
    // acceptance. Physical branching ratios are applied via HNLCalc in analysis.

    if (verbose) {
        std::cout << "==========================================\n" << std::endl;
    }
}

// ==========================================================================
// Main function
// ==========================================================================

int main(int argc, char* argv[]) {

    // -----------------------------------------------------------------------
    // Parse command line arguments
    // -----------------------------------------------------------------------

    if (argc < 3) {
        std::cout << "Usage: " << argv[0] << " <mass_GeV> <flavor> [nEvents] [mode] [qcdMode] [pTHatMin]" << std::endl;
        std::cout << "  mass_GeV:  HNL mass in GeV" << std::endl;
        std::cout << "  flavor:    electron, muon, tau (PBC benchmark BC6/7/8)" << std::endl;
        std::cout << "  nEvents:   optional, default 100000" << std::endl;
        std::cout << "  mode:      optional, 'direct' (default) or 'fromTau' (tau only)" << std::endl;
        std::cout << "  qcdMode:   optional QCD production mode (default: auto)" << std::endl;
        std::cout << "  pTHatMin:  optional pTHat minimum in GeV (default: mode-dependent)" << std::endl;
        std::cout << "\nProduction modes (tau coupling only):" << std::endl;
        std::cout << "  direct:  B/Ds/W → τ N  (mixing at meson/W vertex)" << std::endl;
        std::cout << "  fromTau: B/Ds/W → τ ν, then τ → N X  (mixing at tau decay)" << std::endl;
        std::cout << "  → Both modes are O(U_tau²), combine in analysis for maximum reach" << std::endl;
        std::cout << "\nQCD modes (SOTA for transverse detectors):" << std::endl;
        std::cout << "  auto:      Standard regime-based card selection (default)" << std::endl;
        std::cout << "  hardBc:    Bc production via gg→bb̄/qq̄→bb̄, pTHatMin=15 GeV" << std::endl;
        std::cout << "  hardccbar: Hard cc̄ with pTHatMin (default 10 GeV)" << std::endl;
        std::cout << "  hardbbbar: Hard bb̄ with pTHatMin (default 10 GeV)" << std::endl;
        std::cout << "\nExamples:" << std::endl;
        std::cout << "  " << argv[0] << " 0.3 muon                            # 300 MeV muon-coupled" << std::endl;
        std::cout << "  " << argv[0] << " 2.0 electron                        # 2 GeV electron-coupled" << std::endl;
        std::cout << "  " << argv[0] << " 3.0 tau 100000 direct               # 3 GeV tau, direct" << std::endl;
        std::cout << "  " << argv[0] << " 3.0 tau 100000 fromTau              # 3 GeV tau, from tau" << std::endl;
        std::cout << "  " << argv[0] << " 4.0 muon 500000 direct hardBc       # Bc production mode" << std::endl;
        std::cout << "  " << argv[0] << " 2.0 muon 100000 direct hardccbar 10 # Hard cc̄, pTHat>10" << std::endl;
        std::cout << "  " << argv[0] << " 3.0 muon 100000 direct hardbbbar 15 # Hard bb̄, pTHat>15" << std::endl;
        return 1;
    }

    double mHNL = std::stod(argv[1]);
    std::string flavor = argv[2];
    int nEvents = (argc >= 4) ? std::stoi(argv[3]) : 100000;
    std::string productionMode = (argc >= 5) ? argv[4] : "direct";
    std::string qcdMode = (argc >= 6) ? argv[5] : "auto";
    double pTHatMinUser = (argc >= 7) ? std::stod(argv[6]) : -1.0;

    // Validate production mode
    if (productionMode != "direct" && productionMode != "fromTau") {
        std::cerr << "Error: Invalid production mode '" << productionMode << "'" << std::endl;
        std::cerr << "Must be 'direct' or 'fromTau'" << std::endl;
        return 1;
    }

    // Validate QCD mode
    if (qcdMode != "auto" && qcdMode != "hardBc" &&
        qcdMode != "hardccbar" && qcdMode != "hardbbbar") {
        std::cerr << "Error: Invalid QCD mode '" << qcdMode << "'" << std::endl;
        std::cerr << "Must be 'auto', 'hardBc', 'hardccbar', or 'hardbbbar'" << std::endl;
        return 1;
    }

    // Validate mode-flavor combination
    if (productionMode == "fromTau" && flavor != "tau") {
        std::cerr << "Error: 'fromTau' mode only valid for tau coupling" << std::endl;
        std::cerr << "For electron/muon, use 'direct' mode only" << std::endl;
        return 1;
    }

    // Check kinematic limit for fromTau mode: τ → N X requires mN < mτ
    // Channels close at different thresholds:
    //   τ → ρ N:  m_N < 1.00 GeV
    //   τ → 3π N: m_N < 1.36 GeV
    //   τ → π N:  m_N < 1.64 GeV
    //   τ → μ ν N: m_N < 1.67 GeV
    //   τ → e ν N: m_N < 1.78 GeV (practical limit)
    double mTauMinusE = M_TAU - M_ELECTRON;  // ~1.777 GeV
    if (productionMode == "fromTau" && mHNL > mTauMinusE) {
        std::cerr << "Error: 'fromTau' mode kinematically forbidden for mHNL = "
                  << mHNL << " GeV" << std::endl;
        std::cerr << "τ → N X requires mHNL < mτ ≈ " << mTauMinusE << " GeV" << std::endl;
        std::cerr << "Use 'direct' mode instead for this mass range." << std::endl;
        return 1;
    }

    int leptonID, neutrinoID;
    double mLepton;
    std::string flavorLabel;
    getLeptonInfo(flavor, leptonID, neutrinoID, mLepton, flavorLabel);

    // Determine production regime (tau uses charm/beauty, not kaon)
    std::string regime = getProductionRegime(mHNL, flavorLabel);

    // Override regime for specific QCD modes
    if (qcdMode == "hardBc") {
        regime = "Bc";
    } else if (qcdMode == "hardccbar") {
        regime = "charm";
    } else if (qcdMode == "hardbbbar") {
        regime = "beauty";
    }

    // Determine effective pTHatMin
    double effectivePTHatMin = 0.0;
    if (qcdMode == "hardBc") {
        effectivePTHatMin = (pTHatMinUser > 0.0) ? pTHatMinUser : 15.0;
    } else if (qcdMode == "hardccbar" || qcdMode == "hardbbbar") {
        effectivePTHatMin = (pTHatMinUser > 0.0) ? pTHatMinUser : 10.0;
    }

    std::cout << "============================================" << std::endl;
    std::cout << "HNL Production Simulation" << std::endl;
    std::cout << "============================================" << std::endl;
    std::cout << "HNL mass:        " << mHNL << " GeV" << std::endl;
    std::cout << "Coupling:        " << flavorLabel << " (BC" << (leptonID == 11 ? "6" : leptonID == 13 ? "7" : "8") << ")" << std::endl;
    std::cout << "Production mode: " << regime << std::endl;
    if (flavorLabel == "tau") {
        std::cout << "Tau mode:        " << productionMode << std::endl;
    }
    if (qcdMode != "auto") {
        std::cout << "QCD mode:        " << qcdMode << std::endl;
        std::cout << "pTHatMin:        " << effectivePTHatMin << " GeV" << std::endl;
    }
    std::cout << "Events:          " << nEvents << std::endl;
    std::cout << "============================================\n" << std::endl;
    
    // -----------------------------------------------------------------------
    // Initialize Pythia
    // -----------------------------------------------------------------------
    
    Pythia pythia;

    // Choose appropriate card based on regime
    std::string cardName;
    if (regime == "Bc") {
        cardName = "hnl_Bc.cmnd";
    } else if (regime == "kaon") {
        cardName = "hnl_Kaon.cmnd";
    } else if (regime == "charm") {
        cardName = "hnl_Dmeson.cmnd";
    } else if (regime == "beauty") {
        cardName = "hnl_Bmeson.cmnd";
    } else {
        std::cerr << "ERROR: Unknown regime '" << regime << "'. Cannot select card file." << std::endl;
        return 1;
    }

    // Try to read card from current directory, then from parent/cards
    std::string cardFile = "cards/" + cardName;
    bool cardOK = pythia.readFile(cardFile);

    if (!cardOK) {
        std::cerr << "Warning: Could not read " << cardFile
                  << ", trying ../cards/..." << std::endl;
        cardFile = "../cards/" + cardName;
        cardOK = pythia.readFile(cardFile);
    }

    if (!cardOK) {
        std::cerr << "ERROR: Could not read card file '" << cardName << "'" << std::endl;
        std::cerr << "Searched in: cards/" << cardName << " and ../cards/" << cardName << std::endl;
        std::cerr << "Card files are required for reproducible physics settings." << std::endl;
        std::cerr << "Please ensure the cards/ directory is present." << std::endl;
        return 1;
    }
    std::cout << "Using card file: " << cardFile << std::endl;

    // -----------------------------------------------------------------------
    // Apply QCD mode overrides (pTHatMin for hard-QCD slicing)
    // -----------------------------------------------------------------------
    // For transverse detectors (MATHUSLA, CODEX-b), the relevant kinematic
    // region is high-pT. Standard MinBias/HardQCD with pTHatMin=0 wastes
    // CPU on low-pT events that never reach the detector. Applying a
    // pTHatMin cut enhances statistics in the relevant phase space.
    // -----------------------------------------------------------------------

    if (qcdMode == "hardccbar") {
        // Override: force HardQCD:hardccbar with pTHatMin cut
        pythia.readString("SoftQCD:all = off");
        pythia.readString("HardQCD:all = off");
        pythia.readString("HardQCD:hardccbar = on");
        pythia.readString("PhaseSpace:pTHatMin = " + std::to_string(effectivePTHatMin));
        std::cout << "QCD override: HardQCD:hardccbar with pTHatMin = "
                  << effectivePTHatMin << " GeV" << std::endl;
    } else if (qcdMode == "hardbbbar") {
        // Override: force HardQCD:hardbbbar with pTHatMin cut
        pythia.readString("SoftQCD:all = off");
        pythia.readString("HardQCD:all = off");
        pythia.readString("HardQCD:hardbbbar = on");
        pythia.readString("PhaseSpace:pTHatMin = " + std::to_string(effectivePTHatMin));
        std::cout << "QCD override: HardQCD:hardbbbar with pTHatMin = "
                  << effectivePTHatMin << " GeV" << std::endl;
    } else if (qcdMode == "hardBc") {
        // Bc card already sets gg2bbbar + qqbar2bbbar + pTHatMin=15
        // Apply user override if provided
        if (pTHatMinUser > 0.0) {
            pythia.readString("PhaseSpace:pTHatMin = " + std::to_string(pTHatMinUser));
            std::cout << "Bc mode: pTHatMin overridden to " << pTHatMinUser << " GeV" << std::endl;
        }
    }

    // -----------------------------------------------------------------------
    // Define HNL particle
    // -----------------------------------------------------------------------
    // We define a single HNL state (PDG 9900012). This is appropriate for
    // Majorana HNL where N = N̄. For Dirac HNL interpretation, multiply
    // final yields by factor 2 (or generate both ±9900012).
    // PBC benchmarks BC6/BC7/BC8 assume Majorana.
    //
    // spinType=2 for spin-1/2 fermion (HNL is a sterile neutrino)
    std::ostringstream hnlDef;
    hnlDef << HNL_ID << ":new = N Nbar 2 0 0";
    pythia.readString(hnlDef.str());
    
    hnlDef.str("");
    hnlDef << HNL_ID << ":m0 = " << mHNL;
    pythia.readString(hnlDef.str());
    
    hnlDef.str("");
    hnlDef << HNL_ID << ":tau0 = 1.0e12";  // Stable (very long ctau in mm)
    pythia.readString(hnlDef.str());
    
    hnlDef.str("");
    hnlDef << HNL_ID << ":mayDecay = off";
    pythia.readString(hnlDef.str());
    
    // -----------------------------------------------------------------------
    // Configure decay channels based on production mode
    // -----------------------------------------------------------------------

    if (flavorLabel == "tau" && productionMode == "fromTau") {
        // MODE B: Tau-decay production (tau coupling only)
        // Parents (B/Ds) forced to τν, then τ → N X
        // → Force meson decays to τν (avoids ~95-98% CPU waste)
        // → Force tau decay to N X
        configureMesonDecaysToTauNu(pythia, mHNL);
        configureTauDecays(pythia, mHNL);

    } else {
        // MODE A: Direct production (default for all flavors)
        // For e/μ: K/D/B → ℓ N  (only mode available)
        // For τ:   B/Ds → τ N  (mixing at meson vertex)
        // → Force meson decays to ℓN
        // → Keep tau decays at SM defaults (no τ → N X)

        configureMesonDecays(pythia, leptonID, mHNL, mLepton);
    }
    
    // Set number of events
    pythia.settings.mode("Main:numberOfEvents", nEvents);
    
    // Reduce output
    pythia.readString("Init:showChangedSettings = on");
    pythia.readString("Init:showChangedParticleData = on");
    pythia.readString("Next:numberCount = 10000");
    pythia.readString("Next:numberShowEvent = 0");
    
    // Initialize
    if (!pythia.init()) {
        std::cerr << "Pythia initialization failed!" << std::endl;
        return 1;
    }
    
    // -----------------------------------------------------------------------
    // Open output file
    // -----------------------------------------------------------------------

    std::ostringstream outFileName;
    outFileName << "HNL_" << massToLabel(mHNL) << "GeV_" << flavorLabel << "_" << regime;

    // For tau coupling, distinguish direct vs fromTau production
    if (flavorLabel == "tau") {
        outFileName << "_" << productionMode;
    }

    // For non-auto QCD modes, include mode and pTHatMin in filename
    if (qcdMode != "auto") {
        outFileName << "_" << qcdMode;
        if (effectivePTHatMin > 0.0) {
            outFileName << "_pTHat" << massToLabel(effectivePTHatMin);
        }
    }

    outFileName << ".csv";
    
    std::ofstream outFile(outFileName.str());
    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open output file: " << outFileName.str() << std::endl;
        return 1;
    }
    
    // CSV header
    outFile << "event,weight,hnl_id,parent_pdg,tau_parent_id,pt,eta,phi,p,E,mass,"
            << "prod_x_mm,prod_y_mm,prod_z_mm,beta_gamma" << std::endl;
    
    // -----------------------------------------------------------------------
    // Event loop
    // -----------------------------------------------------------------------
    
    int nHNLfound = 0;
    int nEventsProcessed = 0;
    int nBcFiltered = 0;  // Count HNLs rejected by Bc parent filter

    // In hardBc mode, only accept HNLs from Bc± (541) parents
    const bool filterBcParent = (qcdMode == "hardBc");

    for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
        if (!pythia.next()) continue;
        nEventsProcessed++;

        double weight = pythia.info.weight();

        // Search for HNL in the event
        for (int i = 0; i < pythia.event.size(); ++i) {
            const Particle& p = pythia.event[i];

            if (std::abs(p.id()) != HNL_ID) continue;

            // Find parent
            int parentPdg = findPhysicalParent(pythia.event, i, HNL_ID);

            // Extract tau grandfather (must happen BEFORE Bc filter, since
            // fromTau events have parentPdg==15 and we need tauParentId
            // to decide whether the tau came from a Bc).
            int tauParentId = 0;
            if (std::abs(parentPdg) == 15) {
                int iTop = pythia.event[i].iTopCopy();
                if (iTop < 0 || iTop >= (int)pythia.event.size()) iTop = i;
                int iTau = pythia.event[iTop].mother1();
                if (iTau > 0 && iTau < (int)pythia.event.size()) {
                    int tauParentPdg = findPhysicalParent(pythia.event, iTau, 0);
                    if (tauParentPdg != 0) {
                        tauParentId = std::abs(tauParentPdg);
                    }
                }
            }

            // In Bc mode, only keep HNLs originating from Bc± (PDG 541).
            // Accept both direct production (Bc→ℓN, parentPdg==541) and
            // tau-chain production (Bc→τν, τ→NX, parentPdg==15 with
            // tauParentId==541).  Other B-meson parents (B+, B0, Bs) are
            // handled by the standard beauty mode with their own normalization.
            if (filterBcParent) {
                bool isFromBc = (std::abs(parentPdg) == 541)
                              || (std::abs(parentPdg) == 15 && tauParentId == 541);
                if (!isFromBc) {
                    nBcFiltered++;
                    continue;
                }
            }
            
            // Get production vertex (in mm)
            double prodX = p.xProd();
            double prodY = p.yProd();
            double prodZ = p.zProd();

            // Calculate boost factor (use actual particle mass for robustness)
            double mass = p.m();
            if (mass <= 0.0 || !std::isfinite(mass)) {
                std::cerr << "WARNING: Invalid HNL mass " << mass
                          << " in event " << iEvent
                          << ", using input mass " << mHNL << std::endl;
                mass = mHNL;
            }
            // Sanity check: mass must be positive (guaranteed by construction)
            if (mass <= 0.0) {
                std::cerr << "FATAL: Both p.m() and mHNL are non-positive!" << std::endl;
                return 1;
            }
            // Compute β γ = p / m (NOT the Lorentz factor γ = E / m)
            // This is the quantity needed for decay length calculations: λ = βγ cτ₀
            double betaGamma = p.pAbs() / mass;
            
            // Write to CSV
            outFile << iEvent << ","
                    << weight << ","
                    << p.id() << ","
                    << parentPdg << ","
                    << tauParentId << ","
                    << p.pT() << ","
                    << p.eta() << ","
                    << p.phi() << ","
                    << p.pAbs() << ","
                    << p.e() << ","
                    << p.m() << ","
                    << prodX << ","
                    << prodY << ","
                    << prodZ << ","
                    << betaGamma << std::endl;
            
            nHNLfound++;
        }
    }
    
    outFile.close();

    // -----------------------------------------------------------------------
    // Final statistics
    // -----------------------------------------------------------------------

    pythia.stat();

    const double sigmaGenMb = pythia.info.sigmaGen();
    const double sigmaGenPb = sigmaGenMb * 1e9;

    // Write run metadata for downstream normalization (especially hard-QCD slices).
    // Keeping this in a sidecar avoids bloating per-event CSV rows.
    const std::string metaFileName = outFileName.str() + ".meta.json";
    {
        std::ofstream metaFile(metaFileName);
        if (metaFile.is_open()) {
            metaFile << "{\n";
            metaFile << "  \"mass_GeV\": " << mHNL << ",\n";
            metaFile << "  \"flavor\": \"" << flavorLabel << "\",\n";
            metaFile << "  \"production_mode\": \"" << productionMode << "\",\n";
            metaFile << "  \"regime\": \"" << regime << "\",\n";
            metaFile << "  \"qcd_mode\": \"" << qcdMode << "\",\n";
            metaFile << "  \"pthat_min_gev\": ";
            if (effectivePTHatMin > 0.0) {
                metaFile << effectivePTHatMin;
            } else {
                metaFile << "null";
            }
            metaFile << ",\n";
            metaFile << "  \"sigma_gen_mb\": " << sigmaGenMb << ",\n";
            metaFile << "  \"sigma_gen_pb\": " << sigmaGenPb << ",\n";
            metaFile << "  \"events_requested\": " << nEvents << ",\n";
            metaFile << "  \"events_processed\": " << nEventsProcessed << ",\n";
            metaFile << "  \"hnls_found\": " << nHNLfound << ",\n";
            metaFile << "  \"bc_parent_filter\": " << (filterBcParent ? "true" : "false") << "\n";
            metaFile << "}\n";
        } else {
            std::cerr << "[WARN] Could not write metadata sidecar: " << metaFileName << std::endl;
        }
    }

    std::cout << "\n============================================" << std::endl;
    std::cout << "Summary" << std::endl;
    std::cout << "============================================" << std::endl;
    std::cout << "Events generated:  " << nEventsProcessed << std::endl;
    std::cout << "HNLs found:        " << nHNLfound << std::endl;
    if (filterBcParent) {
        std::cout << "Bc-filtered out:   " << nBcFiltered << " (non-Bc parents rejected)" << std::endl;
    }
    const double effPercent = (nEventsProcessed > 0)
        ? (100.0 * nHNLfound / nEventsProcessed)
        : 0.0;
    std::cout << "Efficiency:        " << effPercent << "%" << std::endl;
    if (qcdMode != "auto") {
        std::cout << "QCD mode:          " << qcdMode << std::endl;
        std::cout << "pTHatMin:          " << effectivePTHatMin << " GeV" << std::endl;
        std::cout << "sigmaGen:          " << sigmaGenMb << " mb (" << sigmaGenPb << " pb)" << std::endl;
    }
    std::cout << "Output file:       " << outFileName.str() << std::endl;
    std::cout << "Metadata file:     " << metaFileName << std::endl;
    std::cout << "============================================" << std::endl;
    
    return 0;
}
