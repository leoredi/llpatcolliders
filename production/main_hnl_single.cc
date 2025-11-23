// main_hnl_single.cc
// Runs ONE mass point passed as an argument.
// Usage: ./main_hnl_single <mass_GeV> [lepton_flavor]
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

// ----------------------------------------------------------------------
// Helper: map flavour string → PDG IDs
// ----------------------------------------------------------------------
std::pair<int, int> getLeptonIDs(const std::string& flavor) {
  std::string lowerFlavor = flavor;
  std::transform(lowerFlavor.begin(), lowerFlavor.end(),
                 lowerFlavor.begin(), ::tolower);

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

// Helper: map lepton PDG → name
std::string getLeptonName(int leptonID) {
  if (leptonID == 11) return "electron";
  if (leptonID == 13) return "muon";
  if (leptonID == 15) return "tau";
  return "unknown";
}

// ----------------------------------------------------------------------
// String helpers for template replacement and filenames
// ----------------------------------------------------------------------
std::string replaceAll(std::string str, const std::string& from,
                       const std::string& to) {
  size_t start_pos = 0;
  while ((start_pos = str.find(from, start_pos)) != std::string::npos) {
    str.replace(start_pos, from.length(), to);
    start_pos += to.length();
  }
  return str;
}

// Generate config file from template, replacing LEPTON_ID / NEUTRINO_ID
bool generateConfigFromTemplate(const std::string& templateFile,
                                const std::string& outputFile,
                                int leptonID,
                                int neutrinoID) {
  std::ifstream input(templateFile);
  if (!input.is_open()) {
    std::cerr << "Error: Could not open template file: "
              << templateFile << std::endl;
    return false;
  }

  std::ofstream output(outputFile);
  if (!output.is_open()) {
    std::cerr << "Error: Could not create output file: "
              << outputFile << std::endl;
    return false;
  }

  std::string line;
  while (std::getline(input, line)) {
    line = replaceAll(line, "LEPTON_ID",   std::to_string(leptonID));
    line = replaceAll(line, "NEUTRINO_ID", std::to_string(neutrinoID));
    output << line << "\n";
  }

  input.close();
  output.close();
  return true;
}

// RAII helper: delete a file automatically when leaving scope
class ScopedFileRemover {
public:
  explicit ScopedFileRemover(const std::string& filename)
    : filename_(filename) {}

  ~ScopedFileRemover() {
    if (!filename_.empty()) {
      std::remove(filename_.c_str());
    }
  }

  ScopedFileRemover(const ScopedFileRemover&) = delete;
  ScopedFileRemover& operator=(const ScopedFileRemover&) = delete;

private:
  std::string filename_;
};

// ----------------------------------------------------------------------
// Helper: turn a mass string (e.g. "0.25") into a unique label "0p25"
// so we don't get collisions between 0.25, 0.3, 0.35 etc.
// ----------------------------------------------------------------------
std::string makeMassLabel(std::string s) {
  for (char& c : s) {
    if      (c == '.') c = 'p';
    else if (c == '+') c = 'p';
    else if (c == '-') c = 'm';
  }
  return s;
}

// ----------------------------------------------------------------------
// Helper: find physical parent ID for an HNL, skipping HNL→HNL copies
// We:
//   1) Go to the "top" copy of this particle with iTopCopy()
//   2) Take its mother1() and ensure it is NOT another HNL
// ----------------------------------------------------------------------
int findPhysicalParentId(const Pythia8::Event& evt,
                         int iLlp,
                         int llp_pdgid) {
  if (iLlp < 0 || iLlp >= evt.size()) return 0;

  const Particle& p = evt[iLlp];

  // Earliest copy of this HNL in the event record
  int iTop = p.iTopCopy();
  if (iTop < 0 || iTop >= evt.size()) return 0;

  const Particle& pTop = evt[iTop];

  int iMother = pTop.mother1();
  if (iMother <= 0 || iMother >= evt.size()) return 0;

  int id = evt[iMother].id();
  if (std::abs(id) == llp_pdgid) {
    // This would mean "HNL with HNL mother" even for the top copy:
    // extremely unlikely in our setup – treat as unknown.
    return 0;
  }

  return id;
}

// ----------------------------------------------------------------------
// main()
// ----------------------------------------------------------------------
int main(int argc, char* argv[]) {

  // 1. PARSE INPUT
  if (argc < 2) {
    std::cerr << "Usage: ./main_hnl_single <mass_GeV> [lepton_flavor]"
              << std::endl;
    std::cerr << "lepton_flavor options: electron, muon (default), tau"
              << std::endl;
    return 1;
  }

  // Numeric mass for physics
  const double mN = std::stod(argv[1]);
  // Original string for file labels (no rounding)
  const std::string massArg   = argv[1];
  const std::string massLabel = makeMassLabel(massArg);

  std::string leptonFlavor = "muon"; // default
  if (argc >= 3) {
    leptonFlavor = argv[2];
  }

  auto [leptonID, neutrinoID] = getLeptonIDs(leptonFlavor);
  if (leptonID == -1) {
    return 1;
  }

  std::string leptonName = getLeptonName(leptonID);
  std::cout << "Running simulation with m_N = " << mN
            << " GeV, lepton = " << leptonName
            << " (PDG " << leptonID << ")" << std::endl;

  // 2. CONFIGURATION
  const int nEvents   = 200000;
  const int llp_pdgid = 9900015;

  // Template files (in production/ directory)
  // Two production regimes following ANUBIS/MATHUSLA/PBC methodology:
  // - Meson: K + D + B combined for m_N < 5 GeV
  // - EW:    W/Z/top for m_N ≥ 5 GeV
  std::string templateMeson = "hnl_Meson_Inclusive_Template.cmnd";
  std::string templateEW    = "hnl_HighMass_Inclusive_Template.cmnd";

  // Generated config files (temporary) – include *full* massLabel to avoid
  // race conditions & rounding collisions.
  std::system("mkdir -p ../tmp");

  std::string cardMeson =
    "../tmp/hnl_Meson_" + leptonName + "_" + massLabel + "_temp.cmnd";
  std::string cardEW =
    "../tmp/hnl_EW_"    + leptonName + "_" + massLabel + "_temp.cmnd";

  // RAII cleanup of temporary files
  ScopedFileRemover cleanupMeson(cardMeson);
  ScopedFileRemover cleanupEW(cardEW);

  // 3. GENERATE CONFIG FILES FROM TEMPLATES
  bool mesonGenerated =
    generateConfigFromTemplate(templateMeson, cardMeson,
                               leptonID, neutrinoID);
  bool ewGenerated =
    generateConfigFromTemplate(templateEW, cardEW,
                               leptonID, neutrinoID);

  if (!mesonGenerated || !ewGenerated) {
    std::cerr << "Error: Failed to generate configuration files from templates"
              << std::endl;
    return 1;
  }

  // 4. SETUP PYTHIA
  Pythia pythia;
  std::string productionMode;

  /*
    TWO-REGIME PRODUCTION MODEL:

    m_N < 5 GeV:  Meson regime — inclusive K, D, B, Λ_c, Λ_b
    m_N ≥ 5 GeV:  EW regime — W, Z, t(t̄) with W/Z → ℓ N

    Cross sections and BRs are imposed later in analysis using HNLCalc;
    here we only generate kinematic *shapes* and record event-level weights.
  */

  if (mN < 5.0) {
    pythia.readFile(cardMeson);
    productionMode = "Meson";
    std::cout << "Production mode: MESON (K/D/B) for m_N < 5 GeV" << std::endl;
  } else {
    pythia.readFile(cardEW);
    productionMode = "EW";
    std::cout << "Production mode: ELECTROWEAK (W/Z/top) for m_N ≥ 5 GeV"
              << std::endl;
  }

  // Enforce HNL mass and stability from C++, so card edits cannot break it.
  pythia.readString(std::to_string(llp_pdgid) + ":m0 = " + std::to_string(mN));
  pythia.readString("9900015:mayDecay = off");
  pythia.readString("9900015:tau0 = 1e6"); // mm, effectively stable

  // Number of events and minimal logging
  pythia.settings.mode("Main:numberOfEvents", nEvents);
  pythia.readString("Next:numberShowEvent = 0");
  pythia.readString("Init:showChangedSettings = off");
  pythia.readString("Init:showChangedParticleData = off");

  // 5. INIT
  if (!pythia.init()) {
    std::cerr << "Pythia initialization failed." << std::endl;
    return 1;
  }

  // 6. OUTPUT FILE
  std::ostringstream fn;
  fn << "HNL_mass_" << massLabel << "_"
     << leptonName << "_" << productionMode << ".csv";

  std::string csvFilename = "../output/csv/simulation/" + fn.str();

  std::ofstream myfile(csvFilename);
  if (!myfile.is_open()) {
    std::cerr << "Error: Could not open output file " << csvFilename
              << std::endl;
    return 1;
  }

  // HNL is stable in Pythia; we save PRODUCTION vertex only.
  myfile << "event,weight,id,parent_id,pt,eta,phi,momentum,energy,mass,"
         << "prod_x_m,prod_y_m,prod_z_m\n";

  // 7. EVENT LOOP
  int nLLPFound = 0;

  for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
    if (!pythia.next()) continue;

    // Relative MC weight (phase-space reweighting)
    // Absolute σ is imposed in analysis.
    const double weight = pythia.info.weight();

    for (int iPrt = 0; iPrt < pythia.event.size(); ++iPrt) {
      Particle& prt = pythia.event[iPrt];

      // Select only HNLs
      if (std::abs(prt.id()) != llp_pdgid) continue;
      if (!prt.isFinal()) continue;

      // Robust parent identification: follow copy chain back to top copy,
      // then take its non-HNL mother as the physical parent.
      int parent_id = findPhysicalParentId(pythia.event, iPrt, llp_pdgid);
      if (parent_id == 0) {
        // Skip weird/ambiguous cases (should be extremely rare).
        continue;
      }

      // For the production vertex, also use the top copy.
      int iTop = prt.iTopCopy();
      if (iTop < 0 || iTop >= pythia.event.size()) {
        iTop = iPrt;
      }
      const Particle& pTop = pythia.event[iTop];

      myfile << iEvent << "," << weight << "," << prt.id() << ","
             << parent_id << ","
             << prt.pT() << "," << prt.eta() << "," << prt.phi() << ","
             << prt.pAbs() << "," << prt.e() << "," << prt.m() << ","
             << pTop.xProd()/1000.0 << ","
             << pTop.yProd()/1000.0 << ","
             << pTop.zProd()/1000.0 << "\n";

      ++nLLPFound;
    }
  }

  myfile.close();

  std::cout << "Mass " << mN << " GeV (" << leptonName << ", "
            << productionMode << "): Done. (" << nLLPFound
            << " HNLs written to " << csvFilename << ")" << std::endl;

  return 0;
}
