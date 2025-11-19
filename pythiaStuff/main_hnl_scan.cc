// main_hnl_scan.cc
// Complete HNL Scanner
// Features: 
// 1. Automatic Mass Regime Switching
// 2. Geometric Checks (Decay Vertex)
// 3. Normalization Checks (Weight + Parent ID)

#include "Pythia8/Pythia.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <set>
#include <cmath> 

using namespace Pythia8;

int main() {

  // ========================================================================
  // 1. CONFIGURATION
  // ========================================================================
  
  std::vector<double> massPoints = { 
      1.0, 2.0, 3.0,        // Low Mass (3-body)
      3.5, 4.0, 4.5,        // Intermediate (2-body)
      10.0, 20.0, 40.0, 60.0, 80.0 // High Mass (Electroweak)
  };

  int nEventsPerPoint = 100000;
  int llp_pdgid = 9900015;

  std::string cardLowMass  = "hnl_LowMass_Inclusive_Complete.cmnd";
  std::string cardHighMass = "hnl_HighMass_Inclusive_Complete.cmnd";

  // ========================================================================
  // 2. SCAN LOOP
  // ========================================================================

  for (double mN : massPoints) {
    
    std::cout << "\n----------------------------------------------------\n";
    std::cout << " PROCESSING MASS POINT: " << mN << " GeV\n";
    std::cout << "----------------------------------------------------\n";

    Pythia pythia;

    // --- STEP A: Load Physics Regime ---
    if (mN < 5.0) {
      pythia.readFile(cardLowMass);
      if (mN > 3.3) {
        // Force 2-body decays for 3.3 < mN < 5.0
        pythia.readString("521:onMode = off");
        pythia.readString("521:addChannel = 1 1.0 100 -13 9900015"); // B+
        pythia.readString("-521:onMode = off");
        pythia.readString("-521:addChannel = 1 1.0 100 13 9900015"); // B-
        pythia.readString("531:onMode = off");
        pythia.readString("531:addChannel = 1 1.0 100 -13 9900015"); // Bs
        pythia.readString("-531:onMode = off");
        pythia.readString("-531:addChannel = 1 1.0 100 13 9900015"); // Bsbar
        pythia.readString("511:onMode = off"); 
        pythia.readString("-511:onMode = off"); 
        pythia.readString("5122:onMode = off");
        pythia.readString("-5122:onMode = off");
      }
    } else {
      pythia.readFile(cardHighMass);
    }

    // --- STEP B: Set Mass & Init ---
    std::ostringstream massCmd;
    massCmd << llp_pdgid << ":m0 = " << mN;
    pythia.readString(massCmd.str());
    pythia.settings.mode("Main:numberOfEvents", nEventsPerPoint);
    pythia.readString("Next:numberShowEvent = 0");

    if (!pythia.init()) {
      std::cerr << "   [ERROR] Pythia initialization failed.\n";
      continue; 
    }

    // --- STEP C: Output Setup ---
    std::ostringstream fn;
    fn << "HNL_mass_" << std::fixed << std::setprecision(1) << mN << ".csv";
    std::string csvFilename = "../output/csv/" + fn.str(); 

    std::ofstream myfile;
    myfile.open(csvFilename);
    
    // ADDED: 'weight' and 'parent_id'
    myfile << "event,weight,id,parent_id,pt,eta,phi,momentum,mass,decay_x_m,decay_y_m,decay_z_m,L_xyz_m,L_xy_m\n";

    // --- STEP D: Event Loop ---
    int nLLPFound = 0;

    for (int iEvent = 0; iEvent < nEventsPerPoint; ++iEvent) {
      if (!pythia.next()) continue;
      
      // Get Event Weight (Critical for Normalization!)
      double weight = pythia.info.weight();

      std::set<int> writtenIDs;

      for (int iPrt = 0; iPrt < pythia.event.size(); ++iPrt) {
        Particle& prt = pythia.event[iPrt];

        if (abs(prt.id()) != llp_pdgid) continue;
        if (writtenIDs.count(prt.id()) > 0) continue;

        // --- 1. GEOMETRY CALCS ---
        double x_dec_m = prt.xDec() / 1000.0;
        double y_dec_m = prt.yDec() / 1000.0;
        double z_dec_m = prt.zDec() / 1000.0;
        double x_prod_m = prt.xProd() / 1000.0;
        double y_prod_m = prt.yProd() / 1000.0;
        double z_prod_m = prt.zProd() / 1000.0;
        double dx = x_dec_m - x_prod_m;
        double dy = y_dec_m - y_prod_m;
        double dz = z_dec_m - z_prod_m;
        double L_xyz_m = sqrt(dx*dx + dy*dy + dz*dz);
        double L_xy_m = sqrt(dx*dx + dy*dy);

        // --- 2. PARENT ID CHECK ---
        // Find the mother (parent) of this HNL to verify production mode
        int parent_id = 0;
        int idxMother = prt.mother1();
        if (idxMother > 0 && idxMother < pythia.event.size()) {
            parent_id = pythia.event[idxMother].id();
        }

        // -------------------------
        myfile << iEvent << ","
               << weight << ","      // Event Weight
               << prt.id() << ","
               << parent_id << ","   // PDG ID of the parent
               << prt.pT() << ","
               << prt.eta() << ","
               << prt.phi() << ","
               << prt.pAbs() << ","
               << prt.m() << ","
               << x_dec_m << ","   
               << y_dec_m << ","   
               << z_dec_m << ","   
               << L_xyz_m << ","   
               << L_xy_m << "\n";  

        writtenIDs.insert(prt.id());
        nLLPFound++;
      }
    }

    myfile.close();
    std::cout << "   -> Done. Wrote " << nLLPFound << " candidates to " << csvFilename << "\n";
  }

  std::cout << "\nAll scans complete.\n";
  return 0;
}