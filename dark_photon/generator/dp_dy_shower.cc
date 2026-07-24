// Pythia8 ISR shower of MadGraph p p > zp (dark-photon Drell-Yan) events.
//
// The parton-level 2->1 hard process produces the A' with pT = 0 (it recoils
// against nothing), so it goes down the beam and misses a transverse detector.
// Showering with initial-state radiation gives the A' its physical transverse
// recoil (Berlin/Kling approach). We read the MG5 LHE, shower with ISR on
// (MPI + hadronization off -- we only need the A' kinematics after ISR), and
// stream the final A' (pdg 3000001) "pT y" to stdout; N_events/N_zp to stderr.
//
// Usage: dp_dy_shower <lhe_file> <m_A_GeV> [seed]
#include "Pythia8/Pythia.h"
#include <cstdio>
using namespace Pythia8;

int main(int argc, char* argv[]) {
  if (argc < 3) { fprintf(stderr, "usage: dp_dy_shower <lhe> <m_A> [seed]\n"); return 2; }
  const char* lhe = argv[1];
  double mA = atof(argv[2]);
  long seed  = (argc > 3) ? atol(argv[3]) : 42;
  const int ZP = 3000001;

  Pythia pythia;
  // Read the MadGraph parton-level events (beam energies come from the LHE).
  pythia.readString("Beams:frameType = 4");
  pythia.readString(std::string("Beams:LHEF = ") + lhe);
  // Teach Pythia the dark photon so it can carry it through the shower:
  // colourless neutral spin-1, made stable (its decay is done downstream in
  // the acceptance layer, not here).
  char pbuf[128];
  snprintf(pbuf, sizeof(pbuf), "3000001:new = Zp Zp 3 0 0 %.6f 0.0 0.0 0.0 0.0", mA);
  pythia.readString(pbuf);
  pythia.readString("3000001:isResonance = false");
  pythia.readString("3000001:mayDecay = off");
  // Shower configuration: ISR gives the recoil pT; drop MPI + hadronization
  // (irrelevant to the A' momentum and much faster).
  pythia.readString("PartonLevel:ISR = on");
  pythia.readString("PartonLevel:FSR = on");
  pythia.readString("PartonLevel:MPI = off");
  pythia.readString("HadronLevel:all = off");
  pythia.readString("Print:quiet = on");
  pythia.readString("Next:numberCount = 0");
  pythia.readString("Random:setSeed = on");
  char sbuf[64]; snprintf(sbuf, sizeof(sbuf), "Random:seed = %ld", seed);
  pythia.readString(sbuf);
  if (!pythia.init()) { fprintf(stderr, "pythia init failed\n"); return 1; }

  long nEv = 0, nZp = 0;
  while (true) {
    if (!pythia.next()) {
      if (pythia.info.atEndOfFile()) break;
      continue;
    }
    ++nEv;
    for (int i = 0; i < pythia.event.size(); ++i) {
      Particle& p = pythia.event[i];
      if (!p.isFinal() || p.idAbs() != ZP) continue;
      double pT = p.pT();
      double y  = p.y();
      if (pT > 200.0 || fabs(y) > 10.0) continue;
      ++nZp;
      printf("%d %.5f %.5f\n", ZP, pT, y);
    }
  }
  fprintf(stderr, "N_EVENTS %ld\n", nEv);
  fprintf(stderr, "N_ZP %ld\n", nZp);
  return 0;
}
