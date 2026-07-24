# BC1 production generators

## Drell--Yan

`dy_pythia.py` is the production path. It compiles `dp_dy_pythia.cc` against
the pinned native Pythia 8.315 and HNL LHAPDF installations, generates the
pure `f fbar -> Z'` contribution with vector couplings
`g_f = epsilon e Q_f`, showers ISR, and writes the `(pT,y)` spectra consumed
by `dark_photon.production`. Docker is not required.

The driver enforces `m_A >= 1.65 GeV`, the NNPDF40 NLO grid's Q validity
floor. It must not be bypassed to fill the 0.65--1.65 GeV gap: that region
requires a resolved/nonperturbative production treatment, not frozen PDFs.

The generated production grid contains 20,000 events at each of 272 masses:
50 MeV spacing from 1.65--3 GeV, 100 MeV plus the original quarter-GeV anchor
points from 3.1--10 GeV, and 250 MeV from 10.25--50 GeV. Every NPZ records its
generator, PDF, mass-stable seed, physical width, cross section, and
integration error. The 50 GeV endpoint is deliberate: extending the
photon-like `e Q_f` approximation through the Z pole would require the full
electroweak HAHM mixing treatment instead.

`dy_madgraph.py` is an independent validation path using the HAHM UFO and the
same HNL MadGraph runner. After converting it to deterministic template cards
and the pinned LHAPDF environment:

- 10 GeV succeeds: `sigma(epsilon=0.01) = 739.3 +/- 1.1 pb`.
- 5 GeV still fails in MadEvent's 2-to-1 phase-space mapping with NaNs and an
  invalid Bjorken-x value.
- Pythia at 10 GeV gives `692.2 +/- 3.0 pb`, 6.4% below the MadGraph control.
  This is retained as a generator/electromagnetic-scheme validation
  difference; the Pythia result is not rescaled to force agreement.

MadGraph validation artifacts are written under
`data/spectra/dy_madgraph_validation/` so they cannot overwrite production
Pythia spectra.

## Light mesons

`make_meson_spectrum.py` and `dp_meson_softqcd.cc` provide the low-mass
pi0/eta/omega parent spectra used for exclusive meson decays.

The canonical BC1 sensitivity campaign scans 0.020--0.200 GeV at 1 MeV
spacing.  The extended 1.65--50 GeV Drell--Yan grid is retained as a validated
production control: its completed scan is many orders of magnitude below the
three-event sensitivity threshold and is therefore not drawn in the zoomed
physics plot.
