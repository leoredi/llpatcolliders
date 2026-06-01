# Provenance of the grids

## Files

    grids/fonll_pp14tev_nnpdf40_nlo_as_01180_fonll_meson_dsdpTdy_pt0-50_y-3to3_central_bottom.dat
    grids/fonll_pp14tev_nnpdf40_nlo_as_01180_fonll_meson_dsdpTdy_pt0-50_y-3to3_central_charm.dat

## How they were generated

- Code: FONLL v1.3.3 (April 2014), upstream commit
  `c7086e49141cf6705cf7a4bc5f7d0b3a38673203` in the maintained mirror,
  with the local modifications in `patches/`.
- PDF backend: LHAPDF 6.5.6.
- Driver: `scripts/generate_meson_grids.py --pdf nlo --quark {bottom,charm}`
  on commit `<tag at release time>` of this repository.

## Bottom grid

Public FONLL B-hadron default fragmentation, Kartvelishvili `(1-z) z^alpha`
with central `alpha = 24.2`, fragmentation fraction 1. Single direct
production stream from the patched `fragmfonll` executable on a 100x100
grid (pT, y) with internal quark-side fragmentation cutoffs left at FONLL
defaults.

Closure against public FONLL v1.3.2 CTEQ6.6 web form, queried 2026-05-29:

    B hadron, pT = 5 GeV, y = 0
        local CTEQ6.6 run :  7.535277e6 pb/GeV
        public CTEQ6.6    :  7.532900e6 pb/GeV
        relative          : +0.032%

The local CTEQ6.6 run uses identical code paths to the published NNPDF4.0
grid; only the LHAPDF set differs. The 3e-4 relative agreement at the
reference point bounds the systematic uncertainty introduced by the
patched code path itself.

## Charm grid

Charm uses the public FONLL D0 convention. The grid is the sum of two
contributions computed separately by the patched code:

1. Direct BCFY pseudoscalar (mode 5), r = 0.1, fragmentation fraction 1.
2. Direct BCFY vector (mode 4), r = 0.1, fragmentation fraction 1,
   convolved with a collinear `D* -> D0` two-body decay using PDG masses
   and BR(D*0 -> D0 pi0) + BR(D*0 -> D0 gamma) + ...

The two streams are combined with a single global weight,
`charm_feeddown_vector_to_direct_weight = 1.21525336412`, fitted to nine
public FONLL CTEQ6.6 D0 reference points at y = 0 and
pT = 1, 5, 8, 15, 22, 29, 36, 43, 50 GeV. After fitting, the maximum
absolute relative residual over those nine points is 1.94e-3.

Closure on the direct stream alone, for reference:

    D* vector,  pT = 5 GeV, y = 0, r = 0.1
        local CTEQ6.6 run :  4.0129e7 pb/GeV
        public CTEQ6.6    :  4.0129e7 pb/GeV
        relative          :  match at the displayed precision

    D0 direct pseudoscalar, pT = 5 GeV, y = 0, r = 0.1 (no feeddown)
        local CTEQ6.6 run :  3.679271e7 pb/GeV
        public CTEQ6.6    :  3.522100e7 pb/GeV
        relative          : +4.46%

The +4.46% offset on the direct stream is what the calibrated feeddown
absorbs. Without the feeddown convolution the grid is not consistent with
the public FONLL D0 convention and should not be used.

## Convention notes

- The output grid is `pT y dsigma/dpT/dy` in pb/GeV. The public web form
  also offers `dsigma/dpT^2/dy` in pb/GeV^2; that differs by a factor of
  `1/(2 pT)`. Consumers must check which convention they expect.
- FONLL returns `(q + qbar)/2`. Studies that count produced mesons must
  multiply by 2.
- Fragmentation frame: `ifrframe = 1` (y=0 frame), the public FONLL
  default for `dsigma/dpT/dy`. On this branch the meson mass `xmh`
  enters only the `szmin` kinematic guard
  (`patches/fragmfonll.f` around line 648). The guard reads
  `(pT^2 + xmh^2) cosh(y)^2 > sh/4`. Across the released grid
  (`pT <= 50 GeV`, `|y| <= 3`, `sh/4 = 4.9e7 GeV^2`) the worst-case
  left-hand side is approximately `2.5e5 GeV^2`, two orders of magnitude
  below the threshold, and the grid is therefore empirically insensitive
  to `xmh` in this phase space. Regenerating with `xmh = M_B0 / M_D0`
  yields byte-identical numerics. Consumers that put the meson on shell
  at sampling time using the PDG mass are kinematically consistent.
  Extending the grid toward the kinematic edge would re-activate the
  guard and require re-verification.

## Charm species interpretation

The released charm grid is a *D0-convention* shape, calibrated to the
public FONLL D0 output. Total charm-meson cross sections downstream
should still use measured species fractions f(D0), f(D+), f(Ds), f(Lambda_c);
the present grid only provides the pT-y shape under the D0 convention.

## What the patches do

See `patches/README.md` for the textual description. The two changes are:

1. Raise the rapidity-grid array dimension `nymax` (and `nymx`) from the
   FONLL default of 100 to 120 in `fonllgrid.f` and `fragmfonll.f`. This
   is needed so a 100-node rapidity grid plus the terminating sentinel
   is accepted; no physics formulas are altered.
2. Add fragmentation options 4 (BCFY vector) and 5/8 (BCFY pseudoscalar)
   in `fragmfonll.f::fragfun`. The BCFY forms are the standard ones from
   Braaten, Cheung, Fleming, Yuan, Phys. Rev. D51 (1995) 4819,
   hep-ph/9408231.
