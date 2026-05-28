# Vendored sources

## FONLL meson tables

Two double-differential heavy-meson production cross sections,
`dsigma/dpT/dy`, generated from the FONLL web interface at LPTHE.

- URL: http://www.lpthe.jussieu.fr/~cacciari/fonll/fonllform.html
- FONLL version: v1.3.2
- Process: pp at sqrt(s) = 14 TeV
- PDF set: CTEQ6.6
- Output: meson-level (`meson = D0`), fragmentation fraction 1
- Grid: 100 pT bins in [0, 50] GeV x 100 y bins in [-3, 3]
- Scale: central (mu_R = mu_F = mu_0)
- Files:
  - `fonll_pp14tev_cteq66_fonll_meson_dsdpTdy_pt0-50_y-3to3_central_charm.dat`
  - `fonll_pp14tev_cteq66_fonll_meson_dsdpTdy_pt0-50_y-3to3_central_bottom.dat`

Species splitting (D0/D+/Ds for charm; B+/B0/Bs for bottom) is applied
downstream by `hnl/production/constants.py::FRAG_C` and `FRAG_B`. The
single FONLL fragmentation choice (`meson = D0`) supplies the pT and y
*shape* only; species fractions are an external input.

PDF choice rationale: matches the Curtin-group MATHUSLA reference files
(`MATHUSLA_LLPfiles_RHN_U{e,mu,tau}`) so that the resulting GARGOYLE
sensitivity uses the same production baseline as the published MATHUSLA
curves. An NNPDF4.0 update is intentionally deferred (one-file swap).

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
