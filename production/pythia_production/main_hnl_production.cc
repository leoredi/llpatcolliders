


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


const int HNL_ID = 9900012;


const std::vector<int> CHARGED_MESONS_2BODY = {
    321,   
    411,   
    431,   
    521,   
    541    
};


const std::vector<int> NEUTRAL_MESONS_3BODY = {
    421,   
    511,   
    531    
};


const std::vector<int> BARYONS_3BODY = {
    4122,  
    5122   
};


const double M_ELECTRON = 0.000511;
const double M_MUON = 0.10566;
const double M_TAU = 1.777;


const std::map<int, double> MESON_MASSES = {
    {130, 0.498},   
    {321, 0.494},   
    {411, 1.870},   
    {421, 1.865},   
    {431, 1.968},   
    {511, 5.280},   
    {521, 5.279},   
    {531, 5.367},   
    {541, 6.275}    
};


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


bool isKinematicallyAllowed2Body(double mParent, double mLepton, double mHNL) {
    return (mHNL < mParent - mLepton);
}


std::string getProductionRegime(double mHNL, const std::string& flavor = "") {
    
    if (flavor == "tau") {
        if (mHNL < 2.0) return "charm";   
        return "beauty";                   
    }
    
    if (mHNL < 0.5) return "kaon";        
    if (mHNL < 2.0) return "charm";       
    return "beauty";                       
}


std::string massToLabel(double mass) {
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(2) << mass;
    std::string label = ss.str();
    for (char& c : label) {
        if (c == '.') c = 'p';
    }
    return label;
}


int findPhysicalParent(const Event& event, int iParticle, int forbiddenId = 0) {
    if (iParticle < 0 || iParticle >= event.size()) return 0;
    
    
    int iTop = event[iParticle].iTopCopy();
    if (iTop < 0 || iTop >= event.size()) iTop = iParticle;
    
    
    int iMother = event[iTop].mother1();
    if (iMother <= 0 || iMother >= event.size()) return 0;
    
    
    if (forbiddenId != 0 && std::abs(event[iMother].id()) == std::abs(forbiddenId)) {
        return 0;
    }
    
    return event[iMother].id();
}


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
    
    
    
    
    
    

    
    
    
    
    
    double mKL = MESON_MASSES.at(130);
    double mPiCharged = 0.140;  
    if (mHNL + mLepton + mPiCharged < mKL) {
        pythia.readString("130:mayDecay = on");
        pythia.readString("130:onMode = off");
        
        pythia.readString("130:addChannel = 1 0.5 0 -211 " + lepBar + " " + hnl);
        
        pythia.readString("130:addChannel = 1 0.5 0 211 " + lep + " " + hnl);
        if (verbose) std::cout << "  K_L -> π ℓ N : ENABLED (3-body, mayDecay forced on)" << std::endl;
        nChannelsConfigured++;
    } else if (verbose) {
        std::cout << "  K_L -> π ℓ N : DISABLED (kinematically forbidden)" << std::endl;
    }

    
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

    
    double mDplus = MESON_MASSES.at(411);
    double mK0 = 0.498;  
    if (mHNL + mLepton + mK0 < mDplus) {
        
        
        pythia.readString("411:addChannel = 1 0.5 0 -311 " + lepBar + " " + hnl);
        pythia.readString("-411:addChannel = 1 0.5 0 311 " + lep + " " + hnl);
        if (verbose) std::cout << "  D± -> K0 ℓ N : ENABLED (3-body)" << std::endl;
        nChannelsConfigured++;
    } else if (verbose) {
        std::cout << "  D± -> K0 ℓ N : DISABLED (kinematically forbidden)" << std::endl;
    }

    
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

    
    double mBplus = MESON_MASSES.at(521);
    double mD0mass = MESON_MASSES.at(421);
    if (mHNL + mLepton + mD0mass < mBplus) {
        
        
        pythia.readString("521:addChannel = 1 0.5 0 -421 " + lepBar + " " + hnl);
        pythia.readString("-521:addChannel = 1 0.5 0 421 " + lep + " " + hnl);
        if (verbose) std::cout << "  B± -> D0 ℓ N : ENABLED (3-body)" << std::endl;
        nChannelsConfigured++;
    } else if (verbose) {
        std::cout << "  B± -> D0 ℓ N : DISABLED (kinematically forbidden)" << std::endl;
    }

    
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
    
    
    double mLb = 5.620;  
    double mLc = 2.286;  
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

    
    double mLambda0 = 1.115;  
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


void configureMesonDecaysToTauNu(Pythia& pythia, double mHNL, bool verbose = true) {

    if (verbose) {
        std::cout << "\n=== Configuring meson → τν decays (for fromTau mode) ===" << std::endl;
        std::cout << "HNL mass: " << mHNL << " GeV" << std::endl;
    }

    int nChannelsConfigured = 0;

    
    
    
    
    double mDs = MESON_MASSES.at(431);
    if (M_TAU < mDs) {
        pythia.readString("431:onMode = off");
        pythia.readString("431:addChannel = 1 1.0 0 -15 16");   
        pythia.readString("-431:onMode = off");
        pythia.readString("-431:addChannel = 1 1.0 0 15 -16");  
        if (verbose) std::cout << "  Ds± → τ ν : ENABLED" << std::endl;
        nChannelsConfigured++;
    }

    
    
    
    
    
    double mBplus = MESON_MASSES.at(521);
    double mD0 = MESON_MASSES.at(421);
    double mDstar0 = 2.007;  

    pythia.readString("521:onMode = off");
    pythia.readString("-521:onMode = off");

    bool bplus_d_ok = (M_TAU + mD0 < mBplus);
    bool bplus_dstar_ok = (M_TAU + mDstar0 < mBplus);

    if (bplus_d_ok && bplus_dstar_ok) {
        
        pythia.readString("521:addChannel = 1 0.4 0 -421 -15 16");   
        pythia.readString("521:addChannel = 1 0.6 0 -423 -15 16");   
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

    
    
    
    double mB0 = MESON_MASSES.at(511);
    double mDminus = MESON_MASSES.at(411);
    double mDstarMinus = 2.010;  

    pythia.readString("511:onMode = off");
    pythia.readString("-511:onMode = off");

    bool b0_d_ok = (M_TAU + mDminus < mB0);
    bool b0_dstar_ok = (M_TAU + mDstarMinus < mB0);

    if (b0_d_ok && b0_dstar_ok) {
        pythia.readString("511:addChannel = 1 0.4 0 -411 -15 16");   
        pythia.readString("511:addChannel = 1 0.6 0 -413 -15 16");   
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

    
    
    
    double mBs = MESON_MASSES.at(531);
    double mDsPlus = MESON_MASSES.at(431);
    double mDsstar = 2.112;  

    pythia.readString("531:onMode = off");
    pythia.readString("-531:onMode = off");

    bool bs_ds_ok = (M_TAU + mDsPlus < mBs);
    bool bs_dsstar_ok = (M_TAU + mDsstar < mBs);

    if (bs_ds_ok && bs_dsstar_ok) {
        pythia.readString("531:addChannel = 1 0.4 0 -431 -15 16");   
        pythia.readString("531:addChannel = 1 0.6 0 -433 -15 16");   
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

    if (verbose) {
        std::cout << "Total τν channels configured: " << nChannelsConfigured << std::endl;
        std::cout << "==========================================\n" << std::endl;
    }
}


void configureTauDecays(Pythia& pythia, double mHNL, bool verbose = true) {
    std::string hnl = std::to_string(HNL_ID);
    double mTau = M_TAU;

    if (verbose) {
        std::cout << "\n=== Configuring tau → N X decays (MODE B: fromTau) ===" << std::endl;
        std::cout << "HNL mass: " << mHNL << " GeV" << std::endl;
    }

    
    pythia.readString("15:onMode = off");
    pythia.readString("-15:onMode = off");

    
    double mPi = 0.140;   
    double mRho = 0.775;  
    double m3Pi = 3 * mPi;
    double mMu = M_MUON;  
    double mE = M_ELECTRON;  

    
    const bool allow_pi = (mHNL + mPi < mTau);      
    const bool allow_rho = (mHNL + mRho < mTau);    
    const bool allow_tripi = (mHNL + m3Pi < mTau);  

    
    const bool allow_mu = (mHNL + mMu < mTau);      
    const bool allow_e = (mHNL + mE < mTau);        

    
    
    
    double pi_weight = 0.0;
    double rho_weight = 0.0;
    double tripi_weight = 0.0;
    double mu_weight = 0.0;
    double e_weight = 0.0;

    if (allow_pi) {
        
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
            
            pi_weight = 0.70;
            mu_weight = 0.15;
            e_weight = 0.15;
        }
    } else if (allow_mu) {
        
        mu_weight = 0.50;
        e_weight = 0.50;
    } else if (allow_e) {
        
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

    
    
    

    if (verbose) {
        std::cout << "==========================================\n" << std::endl;
    }
}


int main(int argc, char* argv[]) {

    
    
    

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

    
    if (productionMode != "direct" && productionMode != "fromTau") {
        std::cerr << "Error: Invalid production mode '" << productionMode << "'" << std::endl;
        std::cerr << "Must be 'direct' or 'fromTau'" << std::endl;
        return 1;
    }

    
    if (productionMode == "fromTau" && flavor != "tau") {
        std::cerr << "Error: 'fromTau' mode only valid for tau coupling" << std::endl;
        std::cerr << "For electron/muon, use 'direct' mode only" << std::endl;
        return 1;
    }

    
    
    
    
    
    
    
    double mTauMinusE = M_TAU - M_ELECTRON;  
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

    
    std::string regime = getProductionRegime(mHNL, flavorLabel);

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
    
    
    
    
    
    Pythia pythia;

    
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
    
    
    
    
    
    
    
    
    
    
    std::ostringstream hnlDef;
    hnlDef << HNL_ID << ":new = N Nbar 2 0 0";
    pythia.readString(hnlDef.str());
    
    hnlDef.str("");
    hnlDef << HNL_ID << ":m0 = " << mHNL;
    pythia.readString(hnlDef.str());
    
    hnlDef.str("");
    hnlDef << HNL_ID << ":tau0 = 1.0e12";  
    pythia.readString(hnlDef.str());
    
    hnlDef.str("");
    hnlDef << HNL_ID << ":mayDecay = off";
    pythia.readString(hnlDef.str());
    
    
    
    

    if (flavorLabel == "tau" && productionMode == "fromTau") {
        
        
        
        
        configureMesonDecaysToTauNu(pythia, mHNL);
        configureTauDecays(pythia, mHNL);

    } else {
        
        
        
        
        

        configureMesonDecays(pythia, leptonID, mHNL, mLepton);
    }
    
    
    pythia.settings.mode("Main:numberOfEvents", nEvents);
    
    
    pythia.readString("Init:showChangedSettings = on");
    pythia.readString("Init:showChangedParticleData = on");
    pythia.readString("Next:numberCount = 10000");
    pythia.readString("Next:numberShowEvent = 0");
    
    
    if (!pythia.init()) {
        std::cerr << "Pythia initialization failed!" << std::endl;
        return 1;
    }
    
    
    
    

    std::ostringstream outFileName;
    outFileName << "HNL_" << massToLabel(mHNL) << "GeV_" << flavorLabel << "_" << regime;

    
    if (flavorLabel == "tau") {
        outFileName << "_" << productionMode;
    }

    outFileName << ".csv";
    
    std::ofstream outFile(outFileName.str());
    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open output file: " << outFileName.str() << std::endl;
        return 1;
    }
    
    
    outFile << "event,weight,hnl_id,parent_pdg,tau_parent_id,pt,eta,phi,p,E,mass,"
            << "prod_x_mm,prod_y_mm,prod_z_mm,beta_gamma" << std::endl;
    
    
    
    
    
    int nHNLfound = 0;
    int nEventsProcessed = 0;
    
    for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
        if (!pythia.next()) continue;
        nEventsProcessed++;
        
        double weight = pythia.info.weight();
        
        
        for (int i = 0; i < pythia.event.size(); ++i) {
            const Particle& p = pythia.event[i];
            
            if (std::abs(p.id()) != HNL_ID) continue;
            
            
            int parentPdg = findPhysicalParent(pythia.event, i, HNL_ID);
            int tauParentId = 0;
            if (std::abs(parentPdg) == 15) {
                int iTop = pythia.event[i].iTopCopy();
                if (iTop < 0 || iTop >= pythia.event.size()) iTop = i;
                int iTau = pythia.event[iTop].mother1();
                if (iTau > 0 && iTau < pythia.event.size()) {
                    int tauParentPdg = findPhysicalParent(pythia.event, iTau, 0);
                    if (tauParentPdg != 0) {
                        tauParentId = std::abs(tauParentPdg);
                    }
                }
            }
            
            
            double prodX = p.xProd();
            double prodY = p.yProd();
            double prodZ = p.zProd();

            
            double mass = p.m();
            if (mass <= 0.0 || !std::isfinite(mass)) {
                std::cerr << "WARNING: Invalid HNL mass " << mass
                          << " in event " << iEvent
                          << ", using input mass " << mHNL << std::endl;
                mass = mHNL;
            }
            
            if (mass <= 0.0) {
                std::cerr << "FATAL: Both p.m() and mHNL are non-positive!" << std::endl;
                return 1;
            }
            
            
            double betaGamma = p.pAbs() / mass;
            
            
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
