# ALPINIST fermionic-ALP decay-width tables (external data product)

`digitized_width_<channel>.txt` are the digitized decay-width tables of the
GKOZ paper — Dalla Valle Garcia, Kahlhoefer, Ovchynnikov, Zaporozhchenko,
*"Phenomenology of axion-like particles with universal fermion couplings —
revisited"*, arXiv:2310.03524 — as shipped by ALPINIST
(<https://github.com/jjerhot/ALPINIST>, directory
`widths/integrated_fermion_2310.03524/`, BSD 3-Clause license).

- Pinned source commit: `COMMIT_SHA.txt` (fetched 2026-07-08).
- Columns: `m_a [GeV]`, `Gamma / (1/f)^2 [GeV^3]` in the BNT
  (arXiv:1708.00443) coupling convention `g_aff = c_f m_f / f` used throughout
  `model.py` (GKOZ's own axis maps as `1/f_here = 2/f_GKOZ`).
- Mass range 0.01–3.01 GeV; physical eta/eta' mixing poles and the 2 m_c
  perturbative onset are visible in `TotalHad`.

Normalisation cross-checks performed when adopting the tables (see the
`model.py` docstring):

1. `BR(a -> mumu)` at m_a = 1 GeV comes out at 9% against the analytic
   leptonic width — GKOZ state "< 10% for m_a >~ 1 GeV".
2. The perturbative quark-level sum (pole m_c) matches the table's own
   2 m_c onset step at ~2.58 GeV to ~5%.
3. At the 3.01 GeV table ceiling the raw quark-level sum is ~1.5x the table
   (charm mass-scheme sensitivity); `model.py` therefore normalises the
   perturbative continuation to the table at the seam.

Channels kept here: `TotalHad` (total hadronic), `2Gamma` (a -> gamma gamma),
and the all-neutral / no-prompt-track channels used for the visible-fraction
split (`3Pi0`, `2K0Pi0`, `2Pi0Eta`, `2Pi0EtaPrim`), plus `3Pi` and `2PiEta`
for reference.

The b -> s a production coefficient `model.CBS_EFF` and the lepton RG factor
`model.CLL_RG` are evaluated with `../../tools/compute_cbs_alpinist.py`, a
BSD-attributed port of ALPINIST's `above_EW` RG implementation of the same
paper, at Lambda_UV = 1 TeV with universal c_f = 1.
