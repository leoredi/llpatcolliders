// ==========================================================================
// main_hnl_production.cc
//
// Publication-quality HNL production simulation for far-detector studies.
// Follows methodology of MATHUSLA, ANUBIS, and Physics Beyond Colliders.
//
// Usage: ./main_hnl_production <mass_GeV> <flavor> [nEvents] [mode]
//   flavor: electron, muon, tau (PBC benchmarks BC6/BC7/BC8)
//   mode: 'direct' (default) or 'fromTau' (tau coupling only)
//
// Production modes (for maximum tau coupling reach):
//   MODE A ("direct"):  B/Ds/W → τ N     (mixing at meson/W vertex)
//   MODE B ("fromTau"): B/Ds/W → τ ν, τ → N X  (mixing at tau decay)
//   → Both modes are O(U_tau²), combine in analysis for maximum sensitivity
//   → Electron and muon use 'direct' mode only
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

const int HNL_ID = 9900015;

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

// Determine production regime based on HNL mass
std::string getProductionRegime(double mHNL) {
    if (mHNL < 0.5) return "kaon";       // Kaon-dominated regime
    if (mHNL < 2.0) return "charm";      // Charm-dominated regime
    return "beauty";                      // Beauty regime (2.0-10.0 GeV)
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
int findPhysicalParent(const Event& event, int iHNL) {
    if (iHNL < 0 || iHNL >= event.size()) return 0;
    
    // Get the top copy of the HNL (earliest in decay chain)
    int iTop = event[iHNL].iTopCopy();
    if (iTop < 0 || iTop >= event.size()) iTop = iHNL;
    
    // Get mother of the top copy
    int iMother = event[iTop].mother1();
    if (iMother <= 0 || iMother >= event.size()) return 0;
    
    // If mother is also HNL, something went wrong
    if (std::abs(event[iMother].id()) == HNL_ID) {
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
    if (isKinematicallyAllowed2Body(MESON_MASSES.at(321), mLepton, mHNL)) {
        pythia.readString("321:onMode = off");
        pythia.readString("321:addChannel = 1 1.0 0 " + lepBar + " " + hnl);
        pythia.readString("-321:onMode = off");
        pythia.readString("-321:addChannel = 1 1.0 0 " + lep + " " + hnl);
        if (verbose) std::cout << "  K± -> ℓ N : ENABLED" << std::endl;
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
    
    if (verbose) {
        std::cout << "Total channels configured: " << nChannelsConfigured << std::endl;
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
//   - "fromTau" mode: Mesons/W decay SM to τν, taus forced to NX
//
// The two samples are combined in the analysis pipeline.
//
// This function configures MODE B (tau-decay production).
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

    // τ- → π- N (2-body, representative mode for acceptance)
    double mPi = 0.140;  // charged pion
    if (mHNL + mPi < mTau) {
        pythia.readString("15:addChannel = 1 1.0 0 -211 " + hnl);
        pythia.readString("-15:addChannel = 1 1.0 0 211 " + hnl);
        if (verbose) std::cout << "  τ → π N : ENABLED" << std::endl;
    } else if (verbose) {
        std::cout << "  τ → π N : DISABLED (kinematically forbidden)" << std::endl;
        std::cout << "  WARNING: No tau decay channels available at this mass!" << std::endl;
    }

    // Note: In reality, τ → N + X has many channels (π, ρ, ℓνν, etc.)
    // We use one representative mode (τ → π N) for geometric acceptance.
    // Physical branching ratios will be applied via HNLCalc in analysis.

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
        std::cout << "Usage: " << argv[0] << " <mass_GeV> <flavor> [nEvents] [mode]" << std::endl;
        std::cout << "  mass_GeV: HNL mass in GeV" << std::endl;
        std::cout << "  flavor: electron, muon, tau (PBC benchmark BC6/7/8)" << std::endl;
        std::cout << "  nEvents: optional, default 100000" << std::endl;
        std::cout << "  mode: optional, 'direct' (default) or 'fromTau' (tau only)" << std::endl;
        std::cout << "\nProduction modes (tau coupling only):" << std::endl;
        std::cout << "  direct:  B/Ds/W → τ N  (mixing at meson/W vertex)" << std::endl;
        std::cout << "  fromTau: B/Ds/W → τ ν, then τ → N X  (mixing at tau decay)" << std::endl;
        std::cout << "  → Both modes are O(U_tau²), combine in analysis for maximum reach" << std::endl;
        std::cout << "\nExamples:" << std::endl;
        std::cout << "  " << argv[0] << " 0.3 muon              # 300 MeV muon-coupled" << std::endl;
        std::cout << "  " << argv[0] << " 2.0 electron          # 2 GeV electron-coupled" << std::endl;
        std::cout << "  " << argv[0] << " 3.0 tau 100000 direct # 3 GeV tau, direct production" << std::endl;
        std::cout << "  " << argv[0] << " 3.0 tau 100000 fromTau # 3 GeV tau, from tau decay" << std::endl;
        return 1;
    }

    double mHNL = std::stod(argv[1]);
    std::string flavor = argv[2];
    int nEvents = (argc >= 4) ? std::stoi(argv[3]) : 100000;
    std::string productionMode = (argc >= 5) ? argv[4] : "direct";

    // Validate production mode
    if (productionMode != "direct" && productionMode != "fromTau") {
        std::cerr << "Error: Invalid production mode '" << productionMode << "'" << std::endl;
        std::cerr << "Must be 'direct' or 'fromTau'" << std::endl;
        return 1;
    }

    // Validate mode-flavor combination
    if (productionMode == "fromTau" && flavor != "tau") {
        std::cerr << "Error: 'fromTau' mode only valid for tau coupling" << std::endl;
        std::cerr << "For electron/muon, use 'direct' mode only" << std::endl;
        return 1;
    }

    int leptonID, neutrinoID;
    double mLepton;
    std::string flavorLabel;
    getLeptonInfo(flavor, leptonID, neutrinoID, mLepton, flavorLabel);

    // Determine production regime
    std::string regime = getProductionRegime(mHNL);

    std::cout << "============================================" << std::endl;
    std::cout << "HNL Production Simulation" << std::endl;
    std::cout << "============================================" << std::endl;
    std::cout << "HNL mass:        " << mHNL << " GeV" << std::endl;
    std::cout << "Coupling:        " << flavorLabel << " (BC" << (leptonID == 11 ? "6" : leptonID == 13 ? "7" : "8") << ")" << std::endl;
    std::cout << "Production mode: " << regime << std::endl;
    if (flavorLabel == "tau") {
        std::cout << "Tau mode:        " << productionMode << std::endl;
    }
    std::cout << "Events:          " << nEvents << std::endl;
    std::cout << "============================================\n" << std::endl;
    
    // -----------------------------------------------------------------------
    // Initialize Pythia
    // -----------------------------------------------------------------------
    
    Pythia pythia;

    // Choose appropriate card based on regime
    std::string cardName;
    if (regime == "kaon") {
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
        std::cerr << "Warning: Could not read " << cardFile << std::endl;
        std::cerr << "Using default settings..." << std::endl;

        // Fallback: set basic parameters directly
        pythia.readString("Beams:idA = 2212");
        pythia.readString("Beams:idB = 2212");
        pythia.readString("Beams:eCM = 14000.");
        pythia.readString("Tune:pp = 14");

        if (regime == "kaon") {
            pythia.readString("SoftQCD:nonDiffractive = on");
        } else if (regime == "charm") {
            pythia.readString("HardQCD:hardccbar = on");
        } else if (regime == "beauty") {
            pythia.readString("HardQCD:hardbbbar = on");
        }
    } else {
        std::cout << "Using card file: " << cardFile << std::endl;
    }
    
    // -----------------------------------------------------------------------
    // Define HNL particle
    // -----------------------------------------------------------------------
    // We define a single HNL state (PDG 9900015). This is appropriate for
    // Majorana HNL where N = N̄. For Dirac HNL interpretation, multiply
    // final yields by factor 2 (or generate both ±9900015).
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
        // Parents (B/Ds/W) decay SM to τν, then τ → N X
        // → Keep meson/W decays at SM defaults
        // → Force tau decay to N X
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

    outFileName << ".csv";
    
    std::ofstream outFile(outFileName.str());
    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open output file: " << outFileName.str() << std::endl;
        return 1;
    }
    
    // CSV header
    outFile << "event,weight,hnl_id,parent_pdg,pt,eta,phi,p,E,mass,"
            << "prod_x_mm,prod_y_mm,prod_z_mm,boost_gamma" << std::endl;
    
    // -----------------------------------------------------------------------
    // Event loop
    // -----------------------------------------------------------------------
    
    int nHNLfound = 0;
    int nEventsProcessed = 0;
    
    for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
        if (!pythia.next()) continue;
        nEventsProcessed++;
        
        double weight = pythia.info.weight();
        
        // Search for HNL in the event
        for (int i = 0; i < pythia.event.size(); ++i) {
            const Particle& p = pythia.event[i];
            
            if (std::abs(p.id()) != HNL_ID) continue;
            
            // Find parent
            int parentPdg = findPhysicalParent(pythia.event, i);
            
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
            double boostGamma = p.e() / mass;
            
            // Write to CSV
            outFile << iEvent << "," 
                    << weight << ","
                    << p.id() << ","
                    << parentPdg << ","
                    << p.pT() << ","
                    << p.eta() << ","
                    << p.phi() << ","
                    << p.pAbs() << ","
                    << p.e() << ","
                    << p.m() << ","
                    << prodX << ","
                    << prodY << ","
                    << prodZ << ","
                    << boostGamma << std::endl;
            
            nHNLfound++;
        }
    }
    
    outFile.close();
    
    // -----------------------------------------------------------------------
    // Final statistics
    // -----------------------------------------------------------------------
    
    pythia.stat();
    
    std::cout << "\n============================================" << std::endl;
    std::cout << "Summary" << std::endl;
    std::cout << "============================================" << std::endl;
    std::cout << "Events generated:  " << nEventsProcessed << std::endl;
    std::cout << "HNLs found:        " << nHNLfound << std::endl;
    std::cout << "Efficiency:        " << (100.0 * nHNLfound / nEventsProcessed) << "%" << std::endl;
    std::cout << "Output file:       " << outFileName.str() << std::endl;
    std::cout << "============================================" << std::endl;
    
    return 0;
}
