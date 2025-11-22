// main_hnl_single.cc
// Runs ONE mass point passed as an argument.
// Usage: ./main_hnl_single <mass_in_GeV> [lepton_flavor]
// lepton_flavor can be: electron, muon (default), tau

#include "Pythia8/Pythia.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <sstream>

using namespace Pythia8;

// Helper function to get lepton PDG ID and neutrino PDG ID from flavor name
std::pair<int, int> getLeptonIDs(const std::string& flavor) {
  std::string lowerFlavor = flavor;
  std::transform(lowerFlavor.begin(), lowerFlavor.end(), lowerFlavor.begin(), ::tolower);

  if (lowerFlavor == "electron" || lowerFlavor == "e") {
    return {11, 12}; // electron, electron neutrino
  } else if (lowerFlavor == "muon" || lowerFlavor == "mu") {
    return {13, 14}; // muon, muon neutrino
  } else if (lowerFlavor == "tau") {
    return {15, 16}; // tau, tau neutrino
  } else {
    std::cerr << "Unknown lepton flavor: " << flavor << std::endl;
    std::cerr << "Valid options: electron, muon, tau" << std::endl;
    return {-1, -1};
  }
}

// Helper function to get lepton name
std::string getLeptonName(int leptonID) {
  if (leptonID == 11) return "electron";
  if (leptonID == 13) return "muon";
  if (leptonID == 15) return "tau";
  return "unknown";
}

// Helper function to replace all occurrences in a string
std::string replaceAll(std::string str, const std::string& from, const std::string& to) {
  size_t start_pos = 0;
  while((start_pos = str.find(from, start_pos)) != std::string::npos) {
    str.replace(start_pos, from.length(), to);
    start_pos += to.length();
  }
  return str;
}

// Helper function to generate config file from template
bool generateConfigFromTemplate(const std::string& templateFile,
                                const std::string& outputFile,
                                int leptonID,
                                int neutrinoID) {
  std::ifstream input(templateFile);
  if (!input.is_open()) {
    std::cerr << "Error: Could not open template file: " << templateFile << std::endl;
    return false;
  }

  std::ofstream output(outputFile);
  if (!output.is_open()) {
    std::cerr << "Error: Could not create output file: " << outputFile << std::endl;
    return false;
  }

  std::string line;
  while (std::getline(input, line)) {
    line = replaceAll(line, "LEPTON_ID", std::to_string(leptonID));
    line = replaceAll(line, "NEUTRINO_ID", std::to_string(neutrinoID));
    output << line << "\n";
  }

  input.close();
  output.close();
  return true;
}

// Class to handle automatic file removal using RAII
class ScopedFileRemover {
public:
  // Constructor: Remembers the filename
  ScopedFileRemover(const std::string& filename) : filename_(filename) {}

  // Destructor: Automatically removes the file when the object goes out of scope
  ~ScopedFileRemover() {
    if (!filename_.empty()) {
      std::remove(filename_.c_str());
    }
  }

  // Prevent copying to avoid confusion about who owns the file
  ScopedFileRemover(const ScopedFileRemover&) = delete;
  ScopedFileRemover& operator=(const ScopedFileRemover&) = delete;

private:
  std::string filename_;
};

int main(int argc, char* argv[]) {

  // 1. PARSE INPUT
  if (argc < 2) {
    std::cerr << "Usage: ./main_hnl_single <mass_GeV> [lepton_flavor]" << std::endl;
    std::cerr << "lepton_flavor options: electron, muon (default), tau" << std::endl;
    return 1;
  }

  double mN = std::stod(argv[1]);
  std::string leptonFlavor = "muon"; // default

  if (argc >= 3) {
    leptonFlavor = argv[2];
  }

  // Get lepton PDG IDs
  auto [leptonID, neutrinoID] = getLeptonIDs(leptonFlavor);
  if (leptonID == -1) {
    return 1;
  }

  std::string leptonName = getLeptonName(leptonID);
  std::cout << "Running simulation with " << leptonName << " (PDG ID: " << leptonID << ")" << std::endl;

  // 2. CONFIGURATION
  int nEvents = 200000;
  int llp_pdgid = 9900015;

  // Template files (in production/ directory)
  // Two production regimes following ANUBIS/MATHUSLA/PBC methodology:
  // - Meson: K + D + B combined for m_N < 5 GeV
  //   (includes kaons, charm, beauty all produced simultaneously)
  // - EW: W/Z/top for m_N ≥ 5 GeV
  std::string templateMeson = "hnl_Meson_Inclusive_Template.cmnd";
  std::string templateEW    = "hnl_HighMass_Inclusive_Template.cmnd";

  // Generated config files (temporary) - include mass to avoid race conditions
  // Store in ../tmp/ directory to keep working directory clean
  std::system("mkdir -p ../tmp");  // Ensure tmp directory exists

  std::ostringstream massStr;
  massStr << std::fixed << std::setprecision(1) << mN;
  std::string cardMeson = "../tmp/hnl_Meson_" + leptonName + "_" + massStr.str() + "_temp.cmnd";
  std::string cardEW    = "../tmp/hnl_EW_" + leptonName + "_" + massStr.str() + "_temp.cmnd";

  // --- RAII CLEANUP OBJECTS ---
  // These will automatically delete the files when main() returns or exits
  ScopedFileRemover cleanupMeson(cardMeson);
  ScopedFileRemover cleanupEW(cardEW);

  // 3. GENERATE CONFIG FILES FROM TEMPLATES
  bool mesonGenerated = generateConfigFromTemplate(templateMeson, cardMeson, leptonID, neutrinoID);
  bool ewGenerated    = generateConfigFromTemplate(templateEW, cardEW, leptonID, neutrinoID);

  if (!mesonGenerated || !ewGenerated) {
    std::cerr << "Error: Failed to generate configuration files from templates" << std::endl;
    return 1; // cleanup destructors run here automatically
  }

  // 4. SETUP PYTHIA
  Pythia pythia;

  /*
   TWO-REGIME PRODUCTION MODEL (ANUBIS/MATHUSLA/PBC-compatible):

   The HNL production hierarchy at √s = 14 TeV includes:
     Below ~0.5 GeV  → KAONS matter
     ~0.5-2 GeV      → CHARM dominates (σ_cc̄ ≈ 2×10^10 pb)
     ~2-5 GeV        → BEAUTY becomes important (σ_bb̄ ≈ 6×10^8 pb)
     > 5 GeV         → ELECTROWEAK (W/Z/top) takes over

   Instead of splitting by mass, we now produce K+D+B SIMULTANEOUSLY
   for all low masses, letting the analysis layer separate contributions
   by parent_id and apply correct BR(M→ℓN) from HNLCalc/PBC tables.
  */

  std::string productionMode;

  if (mN < 5.0) {
    // MESON REGIME: K + D + B combined
    pythia.readFile(cardMeson);
    productionMode = "Meson";
    std::cout << "Production mode: MESON (K/D/B) for m_N < 5 GeV" << std::endl;
  } else {
    // ELECTROWEAK REGIME: W/Z/top
    pythia.readFile(cardEW);
    productionMode = "EW";
    std::cout << "Production mode: ELECTROWEAK (W/Z/top) for m_N >= 5 GeV" << std::endl;
  }

  // --- Logic: Set Mass ---
  pythia.readString(std::to_string(llp_pdgid) + ":m0 = " + std::to_string(mN));
  pythia.settings.mode("Main:numberOfEvents", nEvents);
  pythia.readString("Next:numberShowEvent = 0");              // Silence event logs
  pythia.readString("Init:showChangedSettings = off");        // Reduce clutter
  pythia.readString("Init:showChangedParticleData = off");

  // 5. INIT
  if (!pythia.init()) return 1; // Destructors run here automatically

  // 6. OUTPUT FILE
  // Include production mode in filename for provenance tracking
  std::ostringstream fn;
  fn << "HNL_mass_" << std::fixed << std::setprecision(1) << mN
     << "_" << leptonName << "_" << productionMode << ".csv";

  // SAFETY: Ensure this folder exists!
  std::string csvFilename = "../output/csv/simulation/" + fn.str();

  std::ofstream myfile;
  myfile.open(csvFilename);
  // HNL is stable in Pythia (mayDecay=off), so we only save PRODUCTION vertex
  // Decay coordinates and lengths are invalid/misleading and handled in Python
  myfile << "event,weight,id,parent_id,pt,eta,phi,momentum,energy,mass,prod_x_m,prod_y_m,prod_z_m\n";

  // 7. EVENT LOOP
  int nLLPFound = 0;
  for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
    if (!pythia.next()) continue;

    // IMPORTANT: Use RELATIVE event weight, not absolute cross-section!
    // pythia.info.weight() returns a relative MC weight (phase-space reweighting)
    // Analysis normalizes to external σ from production_xsecs.get_parent_sigma_pb()
    // DO NOT use pythia.info.sigmaGen() here - that would double-count cross-section!
    double weight = pythia.info.weight();

    for (int iPrt = 0; iPrt < pythia.event.size(); ++iPrt) {
      Particle& prt = pythia.event[iPrt];

      // Select only HNLs
      if (std::abs(prt.id()) != llp_pdgid) continue;

      /*
       The HNL is stable in Pythia (mayDecay=off), so the correct physical HNL
       is simply the final-copy particle with id=9900015.
       This replaces the older isLastCopy/isFinal+status logic.
      */
      if (!prt.isFinal()) continue;

      // Parent Check
      int parent_id = 0;
      if (prt.mother1() > 0) parent_id = pythia.event[prt.mother1()].id();

      // Write only VALID data: kinematics and PRODUCTION vertex
      // Decay vertex is meaningless since HNL is stable in Pythia
      myfile << iEvent << "," << weight << "," << prt.id() << "," << parent_id << ","
             << prt.pT() << "," << prt.eta() << "," << prt.phi() << "," << prt.pAbs() << ","
             << prt.e() << "," << prt.m() << ","
             << prt.xProd()/1000.0 << "," << prt.yProd()/1000.0 << "," << prt.zProd()/1000.0 << "\n";

      nLLPFound++;
    }
  }
  myfile.close();

  // Print minimal output to stdout so parallel runner captures it
  std::cout << "Mass " << mN << " GeV (" << leptonName << "): Done. (" << nLLPFound << " HNLs)" << std::endl;

  return 0; // Destructors run here automatically
}
