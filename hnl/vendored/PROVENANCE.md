# Vendored sources

## FONLL meson tables

The default production backend uses two local FONLL+LHAPDF
double-differential heavy-meson production cross sections, `dsigma/dpT/dy`,
generated for the HL-LHC pp 14 TeV setup.

- FONLL source: v1.3.3 (`c7086e49141cf6705cf7a4bc5f7d0b3a38673203`)
- LHAPDF: 6.5.6
- Process: pp at sqrt(s) = 14 TeV
- PDF set: `NNPDF40_nlo_as_01180` (LHAPDF ID 331700)
- Output: meson-level, fragmentation fraction 1
- Bottom convention: public FONLL B-hadron default, Kartvelishvili
  `alpha = 24.2`
- Charm convention: public FONLL `D0`, BCFY `r = 0.1`, with calibrated
  `D* -> D0` feeddown
- Grid: 100 pT nodes over 0..50 GeV, 100 y nodes over -3..3
- Scale: central (mu_R = mu_F = mu_0)
- Fragmentation frame: `ifrframe = 1` (y=0-frame fragmentation, public-FONLL
  default for `dsigma/dpT/dy`). Under this convention `xsecfrag` itself
  does not reference `xmh`; the only `xmh` dependence in `fragmfonll.f`
  on the `ifrframe = 1` branch is the kinematic guard at
  `szmin` (`fragmfonll.f:648`), `(pT^2 + xmh^2) * cosh(y)^2 > sh/4`.
  In the phase space of these tables (`pT <= 50 GeV`, `|y| <= 3`,
  `sh/4 = 4.9e7 GeV^2`) the guard's worst case
  `(50^2 + 5.28^2) * cosh(3)^2 ~ 2.5e5 GeV^2`
  is always satisfied, so the table is insensitive to `xmh` here —
  verified empirically by regenerating with `xmh = M_B0 / M_D0` and
  recovering byte-identical numerics. The sampler at
  `meson_sampler.py:147` reconstructs the on-shell meson four-vector
  with the physical PDG mass, which is the kinematically consistent
  thing to do for this convention. If the grid is ever extended toward
  the kinematic edge of `sh/4`, the guard would become active and this
  insensitivity claim would need to be re-verified.
- Local FONLL capacity patch: `fonllgrid.f` and `fragmfonll.f` rapidity
  array limits raised to 120 nodes so the dense grid and terminating sentinel
  are accepted; no physics formulas changed
- Files:
  - `fonll_pp14tev_nnpdf40_nlo_as_01180_fonll_meson_dsdpTdy_pt0-50_y-3to3_central_charm.dat`
  - `fonll_pp14tev_nnpdf40_nlo_as_01180_fonll_meson_dsdpTdy_pt0-50_y-3to3_central_bottom.dat`

Species splitting (D0/D+/Ds for charm; B+/B0/Bs for bottom) is applied
downstream by `hnl/production/constants.py::FRAG_C` and `FRAG_B`. The
FONLL tables supply pT-y shape and normalization before physical species
fractions. CTEQ6.6 web-generated tables remain vendored as `cteq66_legacy`
comparison inputs in `production/fonll/fonll_parser.py`.

## HNLCalc

Pure-Python computation of HNL production and decay branching ratios.

- Upstream: https://github.com/kaiserds-jr/HNLCalc
- Vendored files: `HNLCalc.py`, `alph_str.csv`, `README.md`
- Stripped from upstream snapshot: `.git/`, notebook examples, embedded
  PNG plots, build artifacts (`__pycache__/`).

Used in `hnl/production/decay_engine/` to evaluate:
- `get_2body_br` and `get_3body_dbr_*` for meson -> HNL channels
- `get_2body_br_tau` and `get_3body_dbr_tau` for tau -> HNL channels
- `integrate_3body_br` for numerical phase-space integration

All HNLCalc calls are made with unit coupling (`U^2 = 1`); the explicit
`U^2` scan is the consumer's responsibility (and not part of this PR).

### Local modifications

- Removed a duplicated `Ds+ -> K0` form-factor block in `HNLCalc.py`.
  Upstream contained two independent `if` statements matching the same
  condition (`pid0 in ["431","-431"] and pid1 in ["311","-311"]`), each
  assigning `f00` a different value: the first set `f00 = 0.747` and a
  later block set `f00 = 0.72`. Because both are plain `if`s, the second
  always overrode the first, making the `0.747` block dead code. We
  deleted the dead `0.747` block and kept the effective value
  (`f00 = 0.72`), so the numerical behavior is unchanged.

  `Ds+ -> K0` is a c -> d transition (with an anti-strange spectator),
  distinct from the `D -> K` c -> s transition immediately above it in
  `HNLCalc.py`. The retained normalization is supported by Melikhov and
  Stech (hep-ph/0001113, Table XVII), which gives
  `F_+^{Ds -> K}(0) = F_0^{Ds -> K}(0) = 0.72` in a relativistic
  dispersion constituent-quark model. For comparison, HPQCD gives
  `f_+^{D -> K}(0) = 0.747(19)` for `D -> K` (1008.4562), while an
  independent lattice calculation gives `f_+^{Ds -> K}(0) = 0.68(4)(3)`
  (0903.1664), compatible with `0.72` within its quoted uncertainties.
  The deleted `0.747` block appears to have copied the `D -> K` value
  into the `Ds -> K` channel.
