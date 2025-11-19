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
#include <set>
#include <cmath>
#include <algorithm>
#include <cstdio>

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
  int nEvents = 100000;
  int llp_pdgid = 9900015;

  // Template files
  std::string templateLowMass  = "hnl_LowMass_Inclusive_Template.cmnd";
  std::string templateHighMass = "hnl_HighMass_Inclusive_Template.cmnd";

  // Generated config files (temporary) - include mass to avoid race conditions
  std::ostringstream massStr;
  massStr << std::fixed << std::setprecision(1) << mN;
  std::string cardLowMass  = "hnl_LowMass_" + leptonName + "_" + massStr.str() + "_temp.cmnd";
  std::string cardHighMass = "hnl_HighMass_" + leptonName + "_" + massStr.str() + "_temp.cmnd";

  // --- RAII CLEANUP OBJECTS ---
  // These will automatically delete the files when main() returns or exits
  ScopedFileRemover cleanupLow(cardLowMass);
  ScopedFileRemover cleanupHigh(cardHighMass);

  // 3. GENERATE CONFIG FILES FROM TEMPLATES
  bool lowMassGenerated = generateConfigFromTemplate(templateLowMass, cardLowMass, leptonID, neutrinoID);
  bool highMassGenerated = generateConfigFromTemplate(templateHighMass, cardHighMass, leptonID, neutrinoID);

  if (!lowMassGenerated || !highMassGenerated) {
    std::cerr << "Error: Failed to generate configuration files from templates" << std::endl;
    return 1; // cleanupLow and cleanupHigh destructors run here automatically
  }

  // 4. SETUP PYTHIA
  Pythia pythia;

  // --- Logic: Determine 3-body to 2-body transition threshold (lepton-dependent) ---
  double threshold_3body_to_2body;
  if (leptonID == 15) {
    // Tau: m_B - m_D - m_tau ~ 5.28 - 1.87 - 1.77 ~ 1.65 GeV
    threshold_3body_to_2body = 1.65;
  } else {
    // Electron/Muon: m_B - m_D - m_mu ~ 3.3 GeV
    threshold_3body_to_2body = 3.3;
  }

  // --- Logic: Choose Regime ---
  if (mN < 5.0) {
    pythia.readFile(cardLowMass);
    if (mN > threshold_3body_to_2body) {
      // Switch 3-body -> 2-body
      pythia.readString("521:onMode = off");
      pythia.readString("521:addChannel = 1 1.0 100 -" + std::to_string(leptonID) + " 9900015");
      pythia.readString("-521:onMode = off");
      pythia.readString("-521:addChannel = 1 1.0 100 " + std::to_string(leptonID) + " 9900015");
      pythia.readString("531:onMode = off");
      pythia.readString("531:addChannel = 1 1.0 100 -" + std::to_string(leptonID) + " 9900015");
      pythia.readString("-531:onMode = off");
      pythia.readString("-531:addChannel = 1 1.0 100 " + std::to_string(leptonID) + " 9900015");
      pythia.readString("511:onMode = off");
      pythia.readString("-511:onMode = off");
      pythia.readString("5122:onMode = off");
      pythia.readString("-5122:onMode = off");
    }
  } else {
    pythia.readFile(cardHighMass);
  }

  // --- Logic: Set Mass ---
  pythia.readString(std::to_string(llp_pdgid) + ":m0 = " + std::to_string(mN));
  pythia.settings.mode("Main:numberOfEvents", nEvents);
  pythia.readString("Next:numberShowEvent = 0"); // Silence event logs
  pythia.readString("Init:showChangedSettings = off"); // Reduce clutter
  pythia.readString("Init:showChangedParticleData = off");

  // 5. INIT
  if (!pythia.init()) return 1; // Destructors run here automatically

  // 6. OUTPUT FILE
  std::ostringstream fn;
  fn << "HNL_mass_" << std::fixed << std::setprecision(1) << mN << "_" << leptonName << ".csv";

  // SAFETY: Ensure this folder exists!
  std::string csvFilename = "csv/" + fn.str();

  std::ofstream myfile;
  myfile.open(csvFilename);
  myfile << "event,weight,id,parent_id,pt,eta,phi,momentum,mass,decay_x_m,decay_y_m,decay_z_m,L_xyz_m,L_xy_m\n";

  // 7. EVENT LOOP
  int nLLPFound = 0;
  for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
    if (!pythia.next()) continue;
    double weight = pythia.info.weight();
    std::set<int> writtenIDs;

    for (int iPrt = 0; iPrt < pythia.event.size(); ++iPrt) {
      Particle& prt = pythia.event[iPrt];
      if (abs(prt.id()) != llp_pdgid) continue;
      if (writtenIDs.count(prt.id()) > 0) continue;

      // Calcs
      double dx = (prt.xDec() - prt.xProd()) / 1000.0;
      double dy = (prt.yDec() - prt.yProd()) / 1000.0;
      double dz = (prt.zDec() - prt.zProd()) / 1000.0;
      double L_xyz = sqrt(dx*dx + dy*dy + dz*dz);
      double L_xy  = sqrt(dx*dx + dy*dy);

      // Parent Check
      int parent_id = 0;
      if (prt.mother1() > 0) parent_id = pythia.event[prt.mother1()].id();

      myfile << iEvent << "," << weight << "," << prt.id() << "," << parent_id << ","
             << prt.pT() << "," << prt.eta() << "," << prt.phi() << "," << prt.pAbs() << ","
             << prt.m() << "," << prt.xDec()/1000.0 << "," << prt.yDec()/1000.0 << "," << prt.zDec()/1000.0 << ","
             << L_xyz << "," << L_xy << "\n";

      writtenIDs.insert(prt.id());
      nLLPFound++;
    }
  }
  myfile.close();

  // Print minimal output to stdout so parallel runner captures it
  std::cout << "Mass " << mN << " GeV (" << leptonName << "): Done. (" << nLLPFound << " events)" << std::endl;
  return 0; // Destructors run here automatically
}
