// main144.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Christian Bierlich <christian.bierlich@fysik.lu.se>

// Keywords: analysis; hepmc; command file; command line option; root; rivet;
//           tuning

// Streamlined event generation with possibility to output ROOT files,
// output HepMC files and run RIVET analyses, all by specifying output modes
// in a cmnd file, where also the event generator settings are specified.
// The example is run with command line options, run ./main144 -h to see a
// full list. See ROOT Usage for information about ROOT output, RIVET Usage
// for information about RIVET and HepMC Interface for information about HepMC.

#include "Pythia8/Pythia.h"
#include "Pythia8/HeavyIons.h"
#include "Pythia8Plugins/InputParser.h"
#include <chrono>
#include <set>
#ifdef RIVET
#include "Pythia8Plugins/Pythia8Rivet.h"
#endif

// Use the Pythia namespace.
using namespace Pythia8;

//==========================================================================

// Implement a user supplied UserHooks derived class inside this
// wrapper, which will allow you to give settings that can be supplied
// in command files.

class UserHooksWrapper : public UserHooks {

public:

  // Add the settings you want available in the run card in this method.
  void additionalSettings(Settings* settingsIn) {
    settings = settingsIn;
    settings->addFlag("UserHooks:doMPICut", false);
    settings->addMode("UserHooks:nMPICut", 0, true, false, 0, 0);
    settings->addFlag("UserHooks:doVetoPartonLevel", false);  // NEW
  }

  // Set the LLP PDG ID to filter for
  void setLLPPdgId(int pdgid) {
    llp_pdgid = pdgid;
  }

  // Check if parton level can be vetoed.
  bool canVetoPartonLevel() final {
    return settings->flag("UserHooks:doMPICut") ||
           settings->flag("UserHooks:doVetoPartonLevel");
  }

  // Check if parton level should be vetoed.
  bool doVetoPartonLevel(const Event& process) final {
    // Existing MPI veto
    if (settings->flag("UserHooks:doMPICut") &&
        infoPtr->nMPI() < settings->mode("UserHooks:nMPICut"))
      return true;

    // NEW: LLP veto - keep only events with specified LLP PDG ID
    if (settings->flag("UserHooks:doVetoPartonLevel")) {
      bool hasLLP = false;
      for (int i = 0; i < process.size(); ++i) {
        if (abs(process[i].id()) == llp_pdgid) {
          hasLLP = true;
          break;
        }
      }
      return !hasLLP;  // Veto if NO LLP found
    }

    return false;
  }

private:

  Settings* settings{};
  int llp_pdgid = 9900015;  // Default to HNL

};

//==========================================================================

// Main example execution.

int main(int argc, char* argv[]) {

  // Parser object for command line input.
  InputParser ip("Run Pythia with cmnd file input, and get Rivet, HepMC or"
    " standard Pythia output.",
    {"./main144 [options]", "./main144 -c main144.cmnd -n 1000 -o myoutput"},
    "Additional options in cmnd file:\n"
    "\tMain:writeLog = on\n\t\tRedirect output to <-o prefix>.log.\n"
    "\tMain:writeHepMC = on \n\t\tWrite HepMC output, requires HepMC linked.\n"
    "\tMain:writeRoot = on \n\t\tWrite an LLP CSV file (LLP.csv).\n"
    "\tMain:runRivet = on \n\t\tRun Rivet analyses, requires Rivet linked.\n"
    "\tMain:rivetAnalyses = {ANALYSIS1,ANALYSIS2,...}\n "
    "\t\tComma separated list of Rivet analyses to run.\n"
    "\t\tAnalysis names can be post-fixed with analysis parameters.\n"
    "\t\tANALYSIS:parm=value:parm2=value2:...\n"
    "\tMain:rivetRunName = STRING \n\t\tAdd an optional run name to "
    "the Rivet analysis.\n"
    "\tMain:rivetIgnoreBeams = on\n\t\tIgnore beams in Rivet. \n"
    "\tMain:rivetDumpPeriod = NUMBER\n\t\tDump Rivet histograms "
    "to file evert NUMBER of events.\n"
    "\tMain:rivetDumpFile = STRING\n\t\t Specify alternative "
    "name for Rivet dump file. Default = OUT.\n");

  // Set up command line options.
  ip.require("c", "User-written command file, can use multiple times.",
    {"-cmnd"});
  ip.add("s", "-1", "Specify seed for the random number generator.",
    {"-seed"});
  ip.add("o", "main144", "Output prefix for log file, Rivet, HepMC, and ROOT.",
    {"-out"});
  ip.add("n", "-1", "Number of events. Overrides the command files.",
    {"-nevents"});
  ip.add("p", "9900015", "PDG ID of LLP to save to CSV.", {"-pdgid"});
  ip.add("l", "false", "Silence the splash screen.");
  ip.add("t", "false", "Time event generation.", {"-time"});
  ip.add("v", "false", "Print Pythia version number and exit.", {"-version"});

  // Initialize the parser and exit if necessary.
  InputParser::Status status = ip.init(argc, argv);
  if (status != InputParser::Valid) return status;

  // Print version number and exit.
  if (ip.get<bool>("v")) {
    cout << "PYTHIA version: " << PYTHIA_VERSION << endl;
    return 0;
  }

  // Get the command files.
  vector<string> cmnds = ip.getVector<string>("c");
  if (cmnds.size() == 0) {
    cout << "Please provide one or more command files with the -c option."
         << endl;
    return 1;
  }

  // Random number seed.
  string seed = ip.get<string>("s");
  // Output filename.
  string out = ip.get<string>("o");
  // Time event generation.
  bool writeTime = ip.get<bool>("t");
  // Command line number of event, overrides the one set in input .cmnd file.
  int nev = ip.get<int>("n");
  // PDG ID of LLP to save
  int llp_pdgid = ip.get<int>("p");

  // Catch the splash screen in a buffer.
  stringstream splashBuf;
  std::streambuf* sBuf = cout.rdbuf();
  cout.rdbuf(splashBuf.rdbuf());
  // The Pythia object.
  Pythia pythia;
  // Direct cout back.
  cout.rdbuf(sBuf);

  // UserHooks wrapper.
  shared_ptr<UserHooksWrapper> userHooksWrapper =
    make_shared<UserHooksWrapper>();
  userHooksWrapper->additionalSettings(&pythia.settings);
  userHooksWrapper->setLLPPdgId(llp_pdgid);
  pythia.setUserHooksPtr(userHooksWrapper);

  // Some extra parameters.
  pythia.settings.addFlag("Main:writeLog", false);
  pythia.settings.addFlag("Main:writeHepMC", false);
  pythia.settings.addFlag("Main:writeRoot", false);
  pythia.settings.addFlag("Main:runRivet", false);
  pythia.settings.addFlag("Main:rivetIgnoreBeams", false);
  pythia.settings.addMode("Main:rivetDumpPeriod", -1, true, false, -1, 0);
  pythia.settings.addWord("Main:rivetDumpFile", "");
  pythia.settings.addWord("Main:rivetRunName", "");
  pythia.settings.addWVec("Main:rivetAnalyses", {});
  pythia.settings.addWVec("Main:rivetPreload", {});

  // Read the command files.
  for (int iCmnd = 0; iCmnd < (int)cmnds.size(); ++iCmnd)
    if (!cmnds[iCmnd].empty()) pythia.readFile(cmnds[iCmnd]);

  // Set seed after reading input.
  if(seed != "-1") {
    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = "+seed);
  }

  // Read the extra parameters.
  if (nev > -1) pythia.settings.mode("Main:numberOfEvents", nev);
  int nEvent                   = pythia.mode("Main:numberOfEvents");;
  int nError                   = pythia.mode("Main:timesAllowErrors");
  bool writeLog                = pythia.flag("Main:writeLog");
  bool writeHepmc              = pythia.flag("Main:writeHepMC");
  bool writeRoot               = pythia.flag("Main:writeRoot"); // now controls CSV output
  bool runRivet                = pythia.flag("Main:runRivet");
  bool countErrors             = nError > 0;

  // Check if Rivet and HepMC are requested and available.
  bool valid = true;
#ifndef RIVET
  valid = valid && !runRivet && !writeHepmc;
  if (runRivet)
    cout << "Option Main::runRivet = on requires the Rivet library.\n";
  if (writeHepmc)
    cout << "Option Main::writeHepMC = on requires the HepMC library.\n";
#endif
  if (!valid) return 1;

  // Rivet and HepMC initialization.
#ifdef RIVET
  // Initialize HepMC.
  Pythia8ToHepMC hepmc;
  if (writeHepmc) hepmc.setNewFile(out + ".hepmc");

  // Initialize Rivet.
  Pythia8Rivet rivet(pythia, out + ".yoda");
  rivet.ignoreBeams(pythia.flag("Main:rivetIgnoreBeams"));
  rivet.dump(pythia.settings.mode("Main:rivetDumpPeriod"),
    pythia.settings.word("Main:rivetDumpFile"));

  // Load the analyses.
  vector<string> rivetAnalyses = pythia.settings.wvec("Main:rivetAnalyses");
  for (int iAna = 0; iAna < (int)rivetAnalyses.size(); ++iAna)
    rivet.addAnalysis(rivetAnalyses[iAna]);

  // Pre-load the YODA histograms.
  vector<string> rivetPreload = pythia.settings.wvec("Main:rivetPreload");
  for (int iYoda = 0; iYoda < (int)rivetPreload.size(); ++iYoda)
    rivet.addPreload(rivetPreload[iYoda]);

  // Add the run name.
  rivet.addRunName(pythia.settings.word("Main:rivetRunName"));
#endif

  // Logfile initialization.
  ofstream logBuf;
  streambuf *oldCout;
  if (writeLog) {
    oldCout = cout.rdbuf(logBuf.rdbuf());
    logBuf.open(out + ".log");
  }

  // Remove splash screen, if requested.
  ostream cnull(NULL);
  if (ip.get<bool>("l")) cnull << splashBuf.str();
  else cout << splashBuf.str();

  // If Pythia fails to initialize, exit with error.
  if (!pythia.init()) return 1;

  // Make a sanity check of initialized Rivet analyses.
#ifdef RIVET
  if (!runRivet && rivetAnalyses.size() > 0 )
    cout << "Rivet analyses are set with Main:rivetAnalyses, "
         << "but Main:runRivet = off.\n";
#endif

  // Loop over events.
  std::cout << "here" << std::endl;
  auto startAllEvents = std::chrono::high_resolution_clock::now();

  // Generate output CSV filename from input cmnd file
  string csvFilename = "LLP.csv";  // default
  // Find first non-empty command file
  for (size_t i = 0; i < cmnds.size(); ++i) {
    if (!cmnds[i].empty()) {
      string cmndFile = cmnds[i];
      // Extract base name (remove directory path)
      size_t lastSlash = cmndFile.find_last_of("/\\");
      if (lastSlash != string::npos) cmndFile = cmndFile.substr(lastSlash + 1);
      // Remove .cmnd extension if present
      size_t dotPos = cmndFile.find_last_of(".");
      if (dotPos != string::npos) cmndFile = cmndFile.substr(0, dotPos);
      csvFilename = cmndFile + "LLP.csv";
      break;
    }
  }

  // Prepend output directory path
  csvFilename = "../output/csv/" + csvFilename;

  ofstream myfile;
  myfile.open(csvFilename);
  cout << "Writing LLP data to: " << csvFilename << endl;
  myfile << "event,\tid,\tpt,\teta,\tphi,\tmomentum,\tmass\n";
  for ( int iEvent = 0; iEvent < nEvent; ++iEvent ) {
    auto startThisEvent = std::chrono::high_resolution_clock::now();

    // Exit if too many failures.
    if (!pythia.next()) {
      if (countErrors && --nError < 0) {
        pythia.stat();
        cout << " \n *-------  PYTHIA STOPPED!  -----------------------*\n"
             << " | Event generation failed due to too many errors. |\n"
             << " *-------------------------------------------------*\n";
        return 1;
      }
      continue;
    }

    // Calculate the event time.
    auto stopThisEvent = std::chrono::high_resolution_clock::now();
    auto eventTime = std::chrono::duration_cast<std::chrono::milliseconds>
      (stopThisEvent - startThisEvent);
    double tt = eventTime.count();

    // Run the Rivet analyses.
#ifdef RIVET
    if (runRivet) {
      if (writeTime) rivet.addAttribute("EventTime", tt);
      rivet();
    }
    if (writeHepmc) hepmc.writeNextEvent(pythia);
#endif

    // Write LLP candidates to CSV file (no ROOT needed).
    if (writeRoot) {
      // Track which charges (PDG IDs with sign) have been written
      std::set<int> writtenCharges;
      const int maxLLPsPerEvent = 2;  // Typically 1 HNL per event, rarely 2 from W+W-

      for (int iPrt = 0; iPrt < pythia.event.size(); ++iPrt) {
        Particle& prt = pythia.event[iPrt];

        // Keep only the desired LLP PDG ID.
        if (abs(prt.id()) != llp_pdgid) continue;

        // Skip if we already wrote a particle with this exact PDG ID (including sign)
        // This eliminates 99.98% of duplicates which have the same charge
        if (writtenCharges.count(prt.id()) > 0) continue;

        // Skip if we've already written too many LLPs for this event
        if (writtenCharges.size() >= maxLLPsPerEvent) break;

        // Write this LLP
        myfile << iEvent << ",\t"
               << prt.id() << ",\t"
               << prt.pT() << ",\t"
               << prt.eta() << ",\t"
               << prt.phi() << ",\t"
               << prt.pAbs() << ",\t"
               << prt.m() << "\n";

        // Mark this charge as written to prevent duplicates
        writtenCharges.insert(prt.id());
      }
    }
  }

  // Finalize.
  myfile.close();
  pythia.stat();

  // Print timing.
  auto stopAllEvents = std::chrono::high_resolution_clock::now();
  auto durationAll = std::chrono::duration_cast<std::chrono::milliseconds>
    (stopAllEvents - startAllEvents);
  if (writeTime) {
    cout << " \n *-------  Generation time  -----------------------*\n"
         << " | Event generation, analysis and writing to files  |\n"
         << " | took: " << double(durationAll.count()) << " ms or "
         << double(durationAll.count())/double(nEvent)
         << " ms per event     |\n"
         << " *-------------------------------------------------*\n";
  }

  // Put cout back in its place.
  if (writeLog) cout.rdbuf(oldCout);
  return 0;

}
