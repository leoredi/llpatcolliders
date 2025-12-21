// ==========================================================================
// main_alp_production.cc
//
// ALP production simulation for transverse LLP detector studies.
//
// This is a kinematic generator: it produces ALP 4-vectors from inclusive SM
// production and forces exotic decays to ALPs. Physical branching ratios are
// applied later in analysis_pbc (Stage 2).
//
// Supported production modes (forced decays):
//   - B_to_Ka        : hard bb̄, force B -> K a
//   - Z_to_gamma_a   : inclusive Z production, force Z -> γ a
//   - h_to_aa        : inclusive Higgs, force h -> a a
//   - h_to_Za        : inclusive Higgs, force h -> Z a
//
// Usage:
//   ./main_alp_production <mass_GeV> <benchmark> <mode> [nEvents]
//
// Output:
//   ALP_<mass>GeV_<benchmark>_<mode>.csv
//
// Notes:
// - We define an ALP particle (PDG 9000005) as stable in Pythia.
// - The forced-decay BRs are set to 100% for kinematics-only generation.
//
// ==========================================================================

#include "Pythia8/Pythia.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>

using namespace Pythia8;

static constexpr int ALP_ID = 9000005;

// Meson masses (GeV) for kinematic checks
static const std::map<int, double> MESON_MASSES = {
    {321, 0.493677},  // K+
    {311, 0.497611},  // K0
    {521, 5.27934},   // B+
    {511, 5.27965},   // B0
    {531, 5.36688},   // Bs
};

std::string massToLabel(double mGeV) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(2) << mGeV;
    std::string s = oss.str();
    for (char& c : s) {
        if (c == '.') c = 'p';
    }
    return s;
}

int findPhysicalParent(const Event& event, int iAlp) {
    if (iAlp < 0 || iAlp >= event.size()) return 0;

    int iTop = event[iAlp].iTopCopy();
    if (iTop < 0 || iTop >= event.size()) iTop = iAlp;

    int iMother = event[iTop].mother1();
    if (iMother <= 0 || iMother >= event.size()) return 0;

    if (std::abs(event[iMother].id()) == ALP_ID) return 0;
    return event[iMother].id();
}

bool isAllowed2Body(double mParent, double m1, double m2) {
    return (mParent > (m1 + m2));
}

void configureBtoKaDecays(Pythia& pythia, double mAlp, bool verbose = true) {
    const std::string alp = std::to_string(ALP_ID);

    const double mBplus = MESON_MASSES.at(521);
    const double mB0 = MESON_MASSES.at(511);
    const double mBs = MESON_MASSES.at(531);
    const double mKplus = MESON_MASSES.at(321);
    const double mK0 = MESON_MASSES.at(311);

    int nChannels = 0;

    if (verbose) {
        std::cout << "\n=== Configuring ALP production channels ===\n";
        std::cout << "ALP mass: " << mAlp << " GeV (PDG " << ALP_ID << ")\n";
        std::cout << "Mode: B -> K a (forced)\n";
    }

    // B+ -> K+ a
    if (isAllowed2Body(mBplus, mKplus, mAlp)) {
        pythia.readString("521:onMode = off");
        pythia.readString("521:addChannel = 1 1.0 0 321 " + alp);
        pythia.readString("-521:onMode = off");
        pythia.readString("-521:addChannel = 1 1.0 0 -321 " + alp);
        if (verbose) std::cout << "  B± -> K± a : ENABLED\n";
        nChannels++;
    } else if (verbose) {
        std::cout << "  B± -> K± a : DISABLED (kinematically forbidden)\n";
    }

    // B0 -> K0 a
    if (isAllowed2Body(mB0, mK0, mAlp)) {
        pythia.readString("511:onMode = off");
        pythia.readString("511:addChannel = 1 1.0 0 311 " + alp);
        pythia.readString("-511:onMode = off");
        pythia.readString("-511:addChannel = 1 1.0 0 -311 " + alp);
        if (verbose) std::cout << "  B0 -> K0 a : ENABLED\n";
        nChannels++;
    } else if (verbose) {
        std::cout << "  B0 -> K0 a : DISABLED (kinematically forbidden)\n";
    }

    // Bs -> K0 a (proxy; used for kinematics-only fast scans)
    if (isAllowed2Body(mBs, mK0, mAlp)) {
        pythia.readString("531:onMode = off");
        pythia.readString("531:addChannel = 1 1.0 0 311 " + alp);
        pythia.readString("-531:onMode = off");
        pythia.readString("-531:addChannel = 1 1.0 0 -311 " + alp);
        if (verbose) std::cout << "  Bs -> K0 a : ENABLED\n";
        nChannels++;
    } else if (verbose) {
        std::cout << "  Bs -> K0 a : DISABLED (kinematically forbidden)\n";
    }

    if (verbose) {
        std::cout << "Total forced channels configured: " << nChannels << "\n";
        std::cout << "==========================================\n\n";
    }
}

void configureZtoGammaADecays(Pythia& pythia, double mAlp, bool verbose = true) {
    const std::string alp = std::to_string(ALP_ID);
    if (verbose) {
        std::cout << "\n=== Configuring Z -> gamma a decays ===\n";
        std::cout << "ALP mass: " << mAlp << " GeV\n";
    }

    // Keep a physical Z width so resonance production remains well-defined.
    pythia.readString("23:mWidth = 2.4952");

    // NOTE:
    // Pythia's built-in Drell-Yan implementation selects allowed Z decay
    // final states at matrix-element level (typically fermion pairs). Exotic
    // decays like Z -> gamma + scalar are not supported as a final-state
    // selector. We therefore generate Z production with a standard decay
    // (mu+mu-) and perform the exotic 2-body decay to an ALP in our own
    // kinematics step (see event loop).
    (void)alp;
    pythia.readString("23:onMode = off");
    pythia.readString("23:onIfAny = 13 -13");
}

void configureHiggsToAa(Pythia& pythia, double mAlp, bool verbose = true) {
    const std::string alp = std::to_string(ALP_ID);
    if (verbose) {
        std::cout << "\n=== Configuring h -> a a decays ===\n";
        std::cout << "ALP mass: " << mAlp << " GeV\n";
    }
    // NOTE:
    // Pythia's SM Higgs production processes similarly select allowed decay
    // final states. For exotic decays we generate Higgs production with a
    // standard decay (bb) and perform the exotic 2-body decay to ALPs in our
    // own kinematics step (see event loop).
    (void)mAlp;
    (void)alp;
    pythia.readString("25:onMode = off");
    pythia.readString("25:onIfAny = 5 -5");
}

void configureHiggsToZa(Pythia& pythia, double mAlp, bool verbose = true) {
    const std::string alp = std::to_string(ALP_ID);
    if (verbose) {
        std::cout << "\n=== Configuring h -> Z a decays ===\n";
        std::cout << "ALP mass: " << mAlp << " GeV\n";
    }
    (void)mAlp;
    (void)alp;
    pythia.readString("25:onMode = off");
    pythia.readString("25:onIfAny = 5 -5");
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <mass_GeV> <benchmark> <mode> [nEvents]\n";
        std::cerr << "Modes: B_to_Ka, Z_to_gamma_a, h_to_aa, h_to_Za\n";
        std::cerr << "Example: " << argv[0] << " 1.00 BC10 B_to_Ka 20000\n";
        return 1;
    }

    const double mAlp = std::stod(argv[1]);
    const std::string benchmark = argv[2];
    const std::string mode = argv[3];
    const int nEvents = (argc >= 5) ? std::stoi(argv[4]) : 20000;

    if (!(mAlp > 0.0 && std::isfinite(mAlp))) {
        std::cerr << "ERROR: invalid mass " << mAlp << "\n";
        return 1;
    }
    if (nEvents <= 0) {
        std::cerr << "ERROR: invalid nEvents " << nEvents << "\n";
        return 1;
    }

    std::cout << "============================================\n";
    std::cout << "ALP Production Simulation (B -> K a)\n";
    std::cout << "============================================\n";
    std::cout << "ALP mass:    " << mAlp << " GeV\n";
    std::cout << "Benchmark:   " << benchmark << "\n";
    std::cout << "Mode:        " << mode << "\n";
    std::cout << "Events:      " << nEvents << "\n";
    std::cout << "============================================\n\n";

    Pythia pythia;

    // Card file
    std::string cardFile;
    if (mode == "B_to_Ka") {
        cardFile = "cards/alp_Bmeson.cmnd";
    } else if (mode == "Z_to_gamma_a") {
        cardFile = "cards/alp_Z.cmnd";
    } else if (mode == "h_to_aa" || mode == "h_to_Za") {
        cardFile = "cards/alp_Higgs.cmnd";
    } else {
        std::cerr << "ERROR: Unknown mode '" << mode << "'\n";
        return 1;
    }
    bool cardOK = pythia.readFile(cardFile);
    if (!cardOK) {
        std::cerr << "Warning: Could not read " << cardFile << ", trying ../cards/...\n";
        cardFile = "../cards/" + cardFile.substr(std::string("cards/").size());
        cardOK = pythia.readFile(cardFile);
    }
    if (!cardOK) {
        std::cerr << "ERROR: Could not read ALP card file for mode " << mode << "\n";
        return 1;
    }
    std::cout << "Using card file: " << cardFile << "\n";

    // Define ALP particle (spinType=1 for spin-0, chargeType=0, colType=0)
    {
        std::ostringstream def;
        def << ALP_ID << ":new = a a 1 0 0";
        pythia.readString(def.str());
    }
    {
        std::ostringstream def;
        def << ALP_ID << ":m0 = " << mAlp;
        pythia.readString(def.str());
    }
    {
        std::ostringstream def;
        def << ALP_ID << ":tau0 = 1.0e12";  // stable in Pythia (mm)
        pythia.readString(def.str());
    }
    {
        std::ostringstream def;
        def << ALP_ID << ":mayDecay = off";
        pythia.readString(def.str());
    }

    if (mode == "B_to_Ka") {
        configureBtoKaDecays(pythia, mAlp);
    } else if (mode == "Z_to_gamma_a") {
        configureZtoGammaADecays(pythia, mAlp);
    } else if (mode == "h_to_aa") {
        configureHiggsToAa(pythia, mAlp);
    } else if (mode == "h_to_Za") {
        configureHiggsToZa(pythia, mAlp);
    }

    pythia.settings.mode("Main:numberOfEvents", nEvents);
    pythia.readString("Init:showChangedSettings = on");
    pythia.readString("Init:showChangedParticleData = on");
    pythia.readString("Next:numberCount = 10000");
    pythia.readString("Next:numberShowEvent = 0");

    if (!pythia.init()) {
        std::cerr << "Pythia initialization failed!\n";
        return 1;
    }

    const std::string outName = "ALP_" + massToLabel(mAlp) + "GeV_" + benchmark + "_" + mode + ".csv";

    std::ofstream out(outName);
    if (!out.is_open()) {
        std::cerr << "Error: Could not open output file: " << outName << "\n";
        return 1;
    }

    out << "event,weight,alp_id,parent_pdg,pt,eta,phi,p,E,mass,prod_x_mm,prod_y_mm,prod_z_mm,beta_gamma\n";

    int nFound = 0;
    int nProcessed = 0;

    const double MZ_NOMINAL = 91.1876;
    auto kallen = [](double a, double b, double c) {
        return a * a + b * b + c * c - 2.0 * a * b - 2.0 * a * c - 2.0 * b * c;
    };
    auto make_isotropic_vec4 = [&](double pAbs, double mass) {
        const double u = 2.0 * pythia.rndm.flat() - 1.0;
        const double phi = 2.0 * M_PI * pythia.rndm.flat();
        const double sinT = std::sqrt(std::max(0.0, 1.0 - u * u));
        const double px = pAbs * sinT * std::cos(phi);
        const double py = pAbs * sinT * std::sin(phi);
        const double pz = pAbs * u;
        const double E = std::sqrt(pAbs * pAbs + mass * mass);
        return Vec4(px, py, pz, E);
    };
    auto find_resonance_index = [&](int pdg) -> int {
        int idx62 = -1;
        for (int i = 0; i < pythia.event.size(); ++i) {
            const Particle& p = pythia.event[i];
            if (p.id() != pdg) continue;
            if (std::abs(p.status()) == 62) idx62 = i;
        }
        if (idx62 >= 0) return idx62;

        int idxAny = -1;
        for (int i = 0; i < pythia.event.size(); ++i) {
            const Particle& p = pythia.event[i];
            if (p.id() != pdg) continue;
            idxAny = i;
        }
        return idxAny;
    };

    for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
        if (!pythia.next()) continue;
        nProcessed++;

        const double weight = pythia.info.weight();

        if (mode == "B_to_Ka") {
            for (int i = 0; i < pythia.event.size(); ++i) {
                const Particle& p = pythia.event[i];
                if (std::abs(p.id()) != ALP_ID) continue;

                const int parentPdg = findPhysicalParent(pythia.event, i);

                const double prodX = p.xProd();
                const double prodY = p.yProd();
                const double prodZ = p.zProd();

                double mass = p.m();
                if (!(mass > 0.0 && std::isfinite(mass))) mass = mAlp;

                const double betaGamma = p.pAbs() / mass;

                out << iEvent << ","
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
                    << betaGamma << "\n";

                nFound++;
            }
            continue;
        }

        if (mode == "Z_to_gamma_a") {
            const int iZ = find_resonance_index(23);
            if (iZ < 0) continue;
            const Particle& Z = pythia.event[iZ];
            const Vec4 pZ = Z.p();
            const double m0 = Z.m();
            if (!(m0 > mAlp && std::isfinite(m0))) continue;

            const double pAbs = (m0 * m0 - mAlp * mAlp) / (2.0 * m0);
            Vec4 pA = make_isotropic_vec4(pAbs, mAlp);
            pA.bst(pZ);

            const double betaGamma = pA.pAbs() / mAlp;

            out << iEvent << ","
                << weight << ","
                << ALP_ID << ","
                << 23 << ","
                << pA.pT() << ","
                << pA.eta() << ","
                << pA.phi() << ","
                << pA.pAbs() << ","
                << pA.e() << ","
                << mAlp << ","
                << Z.xProd() << ","
                << Z.yProd() << ","
                << Z.zProd() << ","
                << betaGamma << "\n";
            nFound++;
            continue;
        }

        if (mode == "h_to_aa" || mode == "h_to_Za") {
            const int iH = find_resonance_index(25);
            if (iH < 0) continue;
            const Particle& H = pythia.event[iH];
            const Vec4 pH = H.p();
            const double m0 = H.m();
            if (!(m0 > 0.0 && std::isfinite(m0))) continue;

            double pAbs = 0.0;
            if (mode == "h_to_aa") {
                if (m0 <= 2.0 * mAlp) continue;
                const double inside = 0.25 * m0 * m0 - mAlp * mAlp;
                if (inside <= 0.0) continue;
                pAbs = std::sqrt(inside);
            } else {
                if (m0 <= (MZ_NOMINAL + mAlp)) continue;
                const double lam = kallen(m0 * m0, MZ_NOMINAL * MZ_NOMINAL, mAlp * mAlp);
                if (lam <= 0.0) continue;
                pAbs = std::sqrt(lam) / (2.0 * m0);
            }

            Vec4 pA = make_isotropic_vec4(pAbs, mAlp);
            pA.bst(pH);

            const double betaGamma = pA.pAbs() / mAlp;

            out << iEvent << ","
                << weight << ","
                << ALP_ID << ","
                << 25 << ","
                << pA.pT() << ","
                << pA.eta() << ","
                << pA.phi() << ","
                << pA.pAbs() << ","
                << pA.e() << ","
                << mAlp << ","
                << H.xProd() << ","
                << H.yProd() << ","
                << H.zProd() << ","
                << betaGamma << "\n";
            nFound++;
            continue;
        }
    }

    out.close();
    pythia.stat();

    std::cout << "\n============================================\n";
    std::cout << "Summary\n";
    std::cout << "============================================\n";
    std::cout << "Events generated:  " << nProcessed << "\n";
    std::cout << "ALPs found:        " << nFound << "\n";
    std::cout << "Output file:       " << outName << "\n";
    std::cout << "============================================\n";

    return 0;
}
