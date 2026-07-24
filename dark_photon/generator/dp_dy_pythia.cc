// Native Pythia 8 dark-photon Drell--Yan generator.
//
// Uses the pure f fbar -> Z' piece with vector couplings chosen so that the
// Pythia Z' vertex equals epsilon * e * Q_f, i.e. a kinetically mixed dark
// photon. NNPDF4.0 is supplied through the pinned HNL LHAPDF installation.
// The A' is kept stable and its showered pT and rapidity are streamed as
// "32 pT y"; summary metadata are written to stderr.

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/LHAPDF6.h"

#include <cmath>
#include <cstdio>
#include <memory>
#include <string>

using namespace Pythia8;

static void setParm(Pythia& pythia, const char* key, double value) {
  char buf[128];
  std::snprintf(buf, sizeof(buf), "%s = %.12g", key, value);
  pythia.readString(buf);
}

int main(int argc, char* argv[]) {
  if (argc < 5) {
    std::fprintf(stderr,
      "usage: dp_dy_pythia <m_A_GeV> <n_events> <epsilon> <seed>\n");
    return 2;
  }
  const double mA = std::atof(argv[1]);
  const int nEvents = std::atoi(argv[2]);
  const double epsilon = std::atof(argv[3]);
  const int seed = std::atoi(argv[4]);
  if (!(mA > 0.) || nEvents <= 0 || !(epsilon > 0.)
      || seed <= 0 || seed > 900000000) {
    std::fprintf(stderr, "invalid argument\n");
    return 2;
  }

  Pythia pythia;
  pythia.readString("Beams:idA = 2212");
  pythia.readString("Beams:idB = 2212");
  pythia.readString("Beams:eCM = 14000.");
  pythia.readString("NewGaugeBoson:ffbar2gmZZprime = on");
  pythia.readString("Zprime:gmZmode = 3");
  pythia.readString("Zprime:universality = on");

  // Pythia normalizes the Z' vector/axial vertex by
  // e/(4 sin(theta_W) cos(theta_W)). Therefore v_f below gives
  // e/(4 sW cW) * v_f = epsilon * e * Q_f.
  const double sin2w = pythia.settings.parm("StandardModel:sin2thetaW");
  const double unit = 4. * epsilon * std::sqrt(sin2w * (1. - sin2w));
  setParm(pythia, "Zprime:vd", -unit / 3.);
  setParm(pythia, "Zprime:vu",  2. * unit / 3.);
  setParm(pythia, "Zprime:ve", -unit);
  setParm(pythia, "Zprime:vnue", 0.);
  setParm(pythia, "Zprime:ad", 0.);
  setParm(pythia, "Zprime:au", 0.);
  setParm(pythia, "Zprime:ae", 0.);
  setParm(pythia, "Zprime:anue", 0.);

  setParm(pythia, "32:m0", mA);
  setParm(pythia, "32:mMin", std::max(0.001, 0.5 * mA));
  setParm(pythia, "32:mMax", 1.5 * mA);
  pythia.readString("32:mayDecay = off");
  setParm(pythia, "PhaseSpace:mHatMin", std::max(0.001, 0.5 * mA));
  setParm(pythia, "PhaseSpace:mHatMax", 1.5 * mA);

  pythia.readString("PartonLevel:ISR = on");
  pythia.readString("PartonLevel:FSR = on");
  pythia.readString("PartonLevel:MPI = off");
  pythia.readString("HadronLevel:all = off");
  pythia.readString("Random:setSeed = on");
  setParm(pythia, "Random:seed", seed);
  pythia.readString("Print:quiet = on");
  pythia.readString("Next:numberCount = 0");

  PDFPtr pdfA = std::make_shared<LHAPDF6>(
    &pythia, &pythia.settings, &pythia.logger);
  PDFPtr pdfB = std::make_shared<LHAPDF6>(
    &pythia, &pythia.settings, &pythia.logger);
  if (!pdfA->init(2212, "NNPDF40_nlo_as_01180", 0, &pythia.logger)
      || !pdfB->init(2212, "NNPDF40_nlo_as_01180", 0, &pythia.logger)) {
    std::fprintf(stderr, "LHAPDF initialization failed\n");
    return 1;
  }
  pythia.setPDFPtr(pdfA, pdfB);
  if (!pythia.init()) {
    std::fprintf(stderr, "Pythia initialization failed\n");
    return 1;
  }

  int nAccepted = 0;
  int nZp = 0;
  for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
    if (!pythia.next()) continue;
    ++nAccepted;
    for (int i = 0; i < pythia.event.size(); ++i) {
      const Particle& particle = pythia.event[i];
      if (!particle.isFinal() || particle.id() != 32) continue;
      ++nZp;
      std::printf("32 %.10g %.10g\n", particle.pT(), particle.y());
      break;
    }
  }

  std::fprintf(stderr, "N_TRIED %d\n", nEvents);
  std::fprintf(stderr, "N_ACCEPTED %d\n", nAccepted);
  std::fprintf(stderr, "N_ZP %d\n", nZp);
  std::fprintf(stderr, "MASS_GEV %.12g\n", mA);
  std::fprintf(stderr, "EPSILON %.12g\n", epsilon);
  std::fprintf(stderr, "WIDTH_GEV %.12g\n", pythia.particleData.mWidth(32));
  std::fprintf(stderr, "SIGMA_MB %.12g\n", pythia.info.sigmaGen());
  std::fprintf(stderr, "SIGMA_ERR_MB %.12g\n", pythia.info.sigmaErr());
  return (nZp > 0) ? 0 : 1;
}
