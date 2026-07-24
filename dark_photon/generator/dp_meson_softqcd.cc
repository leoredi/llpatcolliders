// Pythia8 SoftQCD light-meson spectrum for the GRENDEL BC1 dark-photon
// meson-production channel. pp @ 14 TeV, inelastic. Records produced
// pi0 (111), eta (221), omega (223) by disabling their decay so they appear
// as final-state; the meson -> A' + X two-body decay is done analytically
// downstream. Prints, per species, the inelastic cross section and
// <n_meson>/event to stderr, then a "<PID> pT y" list to stdout.
#include "Pythia8/Pythia.h"
#include <cstdio>
using namespace Pythia8;

int main(int argc, char* argv[]) {
  long nEvents = (argc > 1) ? atol(argv[1]) : 20000;
  long seed    = (argc > 2) ? atol(argv[2]) : 42;

  Pythia pythia;
  pythia.readString("Beams:idA = 2212");
  pythia.readString("Beams:idB = 2212");
  pythia.readString("Beams:eCM = 14000.");
  pythia.readString("SoftQCD:inelastic = on");
  // Keep the produced mesons in the record as final-state particles.
  pythia.readString("111:mayDecay = off");   // pi0
  pythia.readString("221:mayDecay = off");   // eta
  pythia.readString("223:mayDecay = off");   // omega
  pythia.readString("Print:quiet = on");
  pythia.readString("Next:numberCount = 0");
  pythia.readString("Random:setSeed = on");
  char sbuf[64]; snprintf(sbuf, sizeof(sbuf), "Random:seed = %ld", seed);
  pythia.readString(sbuf);
  if (!pythia.init()) { fprintf(stderr, "pythia init failed\n"); return 1; }

  long nPi0 = 0, nEta = 0, nOmega = 0;
  for (long iEv = 0; iEv < nEvents; ++iEv) {
    if (!pythia.next()) continue;
    for (int i = 0; i < pythia.event.size(); ++i) {
      Particle& p = pythia.event[i];
      if (!p.isFinal()) continue;
      int id = p.idAbs();
      if (id != 111 && id != 221 && id != 223) continue;
      double pT = p.pT();
      double y  = p.y();
      if (pT > 20.0 || fabs(y) > 8.0) continue;   // soft/central window
      if (id == 111) ++nPi0; else if (id == 221) ++nEta; else ++nOmega;
      printf("%d %.5f %.5f\n", id, pT, y);
    }
  }
  double sigmaInel_mb = pythia.info.sigmaGen();     // mb
  fprintf(stderr, "SIGMA_INEL_MB %.6f\n", sigmaInel_mb);
  fprintf(stderr, "N_EVENTS %ld\n", nEvents);
  fprintf(stderr, "N_PI0 %ld\n", nPi0);
  fprintf(stderr, "N_ETA %ld\n", nEta);
  fprintf(stderr, "N_OMEGA %ld\n", nOmega);
  fprintf(stderr, "NPI0_PER_EVT %.6e\n",   double(nPi0)   / double(nEvents));
  fprintf(stderr, "NETA_PER_EVT %.6e\n",   double(nEta)   / double(nEvents));
  fprintf(stderr, "NOMEGA_PER_EVT %.6e\n", double(nOmega) / double(nEvents));
  return 0;
}
