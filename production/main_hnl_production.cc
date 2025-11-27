// ==========================================================================
// main_hnl_production.cc
//
// Publication-quality HNL production simulation for far-detector studies.
// Follows methodology of MATHUSLA, ANUBIS, and Physics Beyond Colliders.
//
// Usage: ./main_hnl_production <mass_GeV> <flavor> [nEvents]
//   flavor: electron, muon, tau
//
// Output: CSV file with HNL 4-vectors and parent information
//
// References:
//   - arXiv:1805.08567 (HNL phenomenology)
//   - arXiv:1901.09966 (PBC benchmarks)
//   - arXiv:2103.11494 (Pythia validation for HNL)
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

// W/Z for EW production
const int W_PLUS = 24;
const int Z_BOSON = 23;

// Tau lepton (for BC8)
const int TAU = 15;

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
    if (mHNL < 5.0) return "beauty";     // Beauty-dominated regime
    return "ew";                          // Electroweak regime
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

void configureMesonDecays(Pythia& pythia, int leptonID, int neutrinoID, 
                          double mHNL, double mLepton, bool verbose = true) {
    
    std::string hnl = std::to_string(HNL_ID);
    std::string lep = std::to_string(leptonID);
    std::string lepBar = std::to_string(-leptonID);
    std::string nu = std::to_string(neutrinoID);
    
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
// Configure W/Z decays for EW production
// ==========================================================================

void configureEWDecays(Pythia& pythia, int leptonID, int neutrinoID, 
                       double mHNL, double mLepton, bool verbose = true) {
    
    std::string hnl = std::to_string(HNL_ID);
    std::string lep = std::to_string(leptonID);
    std::string lepBar = std::to_string(-leptonID);
    std::string nu = std::to_string(neutrinoID);
    std::string nuBar = std::to_string(-neutrinoID);
    
    if (verbose) {
        std::cout << "\n=== Configuring EW decay channels ===" << std::endl;
    }
    
    // W+ -> ℓ+ N
    double mW = 80.4;
    if (mHNL + mLepton < mW) {
        pythia.readString("24:onMode = off");
        pythia.readString("24:addChannel = 1 1.0 0 " + lepBar + " " + hnl);
        pythia.readString("-24:onMode = off");
        pythia.readString("-24:addChannel = 1 1.0 0 " + lep + " " + hnl);
        if (verbose) std::cout << "  W± -> ℓ N : ENABLED" << std::endl;
    }
    
    // Z -> ν N (for neutral current production)
    // Note: This is subdominant but physically correct
    double mZ = 91.2;
    if (mHNL < mZ) {
        pythia.readString("23:onMode = off");
        pythia.readString("23:addChannel = 1 1.0 0 " + nu + " " + hnl);
        pythia.readString("23:addChannel = 1 1.0 0 " + nuBar + " " + hnl);
        if (verbose) std::cout << "  Z -> ν N : ENABLED" << std::endl;
    }
    
    // Also configure top -> W b -> ℓ N b
    // The W from top decay will use the W decay channel above
    if (verbose) {
        std::cout << "  t -> W b -> ℓ N b : Uses W channel above" << std::endl;
        std::cout << "==========================================\n" << std::endl;
    }
}

// ==========================================================================
// Configure tau decays for BC8 (tau-coupled HNL)
// ==========================================================================

void configureTauDecays(Pythia& pythia, double mHNL, bool verbose = true) {
    
    std::string hnl = std::to_string(HNL_ID);
    double mTau = M_TAU;
    
    if (verbose) {
        std::cout << "\n=== Configuring tau decay channels ===" << std::endl;
    }
    
    // τ -> π N (2-body)
    double mPi = 0.140;
    if (mHNL + mPi < mTau) {
        pythia.readString("15:onMode = off");
        pythia.readString("15:addChannel = 1 1.0 0 -211 " + hnl);  // τ- -> π- N
        pythia.readString("-15:onMode = off");
        pythia.readString("-15:addChannel = 1 1.0 0 211 " + hnl);  // τ+ -> π+ N
        if (verbose) std::cout << "  τ -> π N : ENABLED" << std::endl;
    }
    
    // τ -> ρ N (2-body, ρ -> ππ)
    double mRho = 0.775;
    if (mHNL + mRho < mTau) {
        // Note: Pythia will decay ρ -> ππ automatically
        pythia.readString("15:addChannel = 1 0.5 0 -213 " + hnl);  // τ- -> ρ- N
        pythia.readString("-15:addChannel = 1 0.5 0 213 " + hnl);  // τ+ -> ρ+ N
        if (verbose) std::cout << "  τ -> ρ N : ENABLED" << std::endl;
    }
    
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
        std::cout << "Usage: " << argv[0] << " <mass_GeV> <flavor> [nEvents]" << std::endl;
        std::cout << "  mass_GeV: HNL mass in GeV" << std::endl;
        std::cout << "  flavor: electron, muon, tau (PBC benchmark BC6/7/8)" << std::endl;
        std::cout << "  nEvents: optional, default 100000" << std::endl;
        std::cout << "\nExamples:" << std::endl;
        std::cout << "  " << argv[0] << " 0.3 muon          # 300 MeV muon-coupled" << std::endl;
        std::cout << "  " << argv[0] << " 2.0 electron      # 2 GeV electron-coupled" << std::endl;
        std::cout << "  " << argv[0] << " 10.0 muon 500000  # 10 GeV muon-coupled, 500k events" << std::endl;
        return 1;
    }
    
    double mHNL = std::stod(argv[1]);
    std::string flavor = argv[2];
    int nEvents = (argc >= 4) ? std::stoi(argv[3]) : 100000;
    
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
        cardName = "hnl_EW.cmnd";
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
        } else {
            pythia.readString("WeakSingleBoson:ffbar2W = on");
            pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
            pythia.readString("Top:gg2ttbar = on");
            pythia.readString("Top:qqbar2ttbar = on");
        }
    } else {
        std::cout << "Using card file: " << cardFile << std::endl;
    }
    
    // Define HNL particle (spinType=1 for fermion, not 2 for vector boson)
    std::ostringstream hnlDef;
    hnlDef << HNL_ID << ":new = N Nbar 1 0 0";
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
    
    // Configure decay channels
    if (regime == "kaon" || regime == "charm" || regime == "beauty") {
        configureMesonDecays(pythia, leptonID, neutrinoID, mHNL, mLepton);
    } else {
        configureEWDecays(pythia, leptonID, neutrinoID, mHNL, mLepton);
    }
    
    // For tau-coupled HNL (BC8), also configure tau decays
    if (leptonID == 15) {
        configureTauDecays(pythia, mHNL);
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
    outFileName << "HNL_" << massToLabel(mHNL) << "GeV_" << flavorLabel << "_" << regime << ".csv";
    
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
            if (mass <= 0.0) mass = mHNL;          // fallback to input mass
            if (mass <= 0.0) mass = 1e-6;          // hard floor to avoid div by zero
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
