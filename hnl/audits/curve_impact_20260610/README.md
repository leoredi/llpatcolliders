# Missing-input curve-impact audit

This audit estimates how the known missing inputs in `REMAINING_WORK.md` could
move the current GRENDEL HNL exclusion boundaries. It covers all 116 configured
masses for `Ue`, `Umu`, and `Utau` (348 flavor/mass points) and 26 effects.

The baseline is the saved `full_20260606_all` run. That run predates the new
`Bbaryon` channel and the current 11-source induced-tau implementation, so those
changes are included as a finite-statistics current-code rerun proxy rather
than silently folded into the baseline.

Calculations used:

```text
/Volumes/sandbox/conda/envs/hnl/bin/python -P
```

## Read this first

- Percentages are relative shifts in an exclusion boundary `U^2`.
- A negative lower-boundary shift improves reach; a positive shift weakens it.
- A positive upper-boundary shift extends the large-mixing edge.
- `NaN` means that the boundary is open, absent, or lost in that scenario.
- Most rows are broad scenarios or reweighting proxies, not confidence
  intervals and not measured biases.
- Detector-efficiency and background rows show sensitivity to assumptions
  that do not yet exist, not an uncertainty estimate.
- The combined physics-model RSS is an engineering envelope. Correlations and
  unknowable detector/background tails are not represented.
- The requested 0.3 MeV spacing is finer than the configured mass grid, whose
  minimum spacing is 15 MeV. Detailed results are therefore provided at every
  actual mass point, with a separate 0.3 GeV summary.

## Data products

- `per_mass_boundary_impacts.csv`: every effect at every available
  flavor/mass/boundary point.
- `summary_0p3GeV_bins.csv`: medians and extrema in 0.3 GeV mass bins.
- `baseline_channel_fractions.csv`: production and accepted-yield channel
  fractions used by the reweighting study.
- `csv_precision_6_vs_8.json`: focused serialization-precision test.
- `diagnostics.json`: baseline run and table dimensions.

The detailed CSV columns include the baseline boundary, whether it is open,
the effect classification, low/high scenario shifts, boundary-existence flags,
and the scenario definition.

## Main findings

- **Charged-kaon transport below about 0.4 GeV:** retaining only 10--70% of
  the current prompt-IP charged-kaon yield weakens the lower boundary by about
  19--212%. This is the largest identified production-model effect at low mass.
- **Neutral kaons:** the additive proxy can improve the lower boundary by up
  to about 19%.
- **Direct form-factor/rate model:** the legacy-rate envelope moves resolved
  boundaries by up to about 53%.
- **HNL hadron/parton matching near 1 GeV:** the proxy reaches about 22%.
- **Missing charm baryons:** the additive proxy reaches about 17%.
- **Exclusive-channel closure and weak cascades:** the proxy reaches about
  15%.
- **Missing multi-hadron tau modes:** the proxy reaches about 6%.
- **Production Monte Carlo statistics:** median boundary effect about 2.6%,
  maximum about 15.2%.
- **Higgs-mediated production:** below 0.14% in the upper-bound proxy and
  negligible for this setup.
- **CSV precision:** degrading eight to six significant digits changed no
  geometry-hit decisions in the tested samples; the largest resolved boundary
  shift was 0.00157%. Eight significant digits are adequate at the tested
  points.

The current-code rerun proxy improves the `Utau` lower boundary by roughly
16--30% below 1 GeV. For `Ue` and `Umu`, the added bottom-baryon contribution
reaches roughly 11% near 2 GeV.

## Combined model envelope

The typical lower-boundary engineering envelope is:

| Region | Approximate relative shift |
|---|---:|
| `Ue`/`Umu`, kaon-dominated region | `-34%` to `+234%` |
| 0.5--3 GeV | roughly `-26%` to `-34%`, and `+40%` to `+54%` |
| around 4 GeV | roughly `-38%` to `-44%`, and `+28%` to `+61%` |

This table is not a combined experimental uncertainty. It excludes
unquantified detector and background tails and combines selected physics
proxies in log-boundary space.

## Open upper boundaries

The baseline has 30 open `Ue`, 30 open `Umu`, and 34 open `Utau` upper
boundaries. The saved `U^2 = 0.1` values at those points are scan endpoints,
not physical exclusion edges. The detailed table leaves their baseline
boundary blank and sets `baseline_open = True`.

## Scenario construction

Directly evaluated or bounded:

- current-code baryon and induced-tau rerun proxy using up to 5,000
  ray-cast events per new channel;
- threshold response for 5- and 10-event alternatives to the current
  three-event threshold;
- weighted production-MC variance;
- Higgs-production upper-bound proxy;
- six-versus-eight-significant-digit serialization test.

Literature-informed ranges:

- FONLL charm and bottom normalization/shape placeholders;
- fragmentation-fraction extrapolation;
- flat electroweak K-factor and mixed W/Z treatment;
- induced-tau branching inputs;
- visible branching fraction.

Unknowable until the missing input exists:

- detector/reconstruction efficiency;
- visible-decay topology acceptance;
- charged-kaon transport;
- kaon and Bc spectra/normalizations;
- neutral-kaon production;
- a real background threshold.

Proxy reweights:

- missing charm baryons;
- exclusive-channel closure and weak cascades;
- spin/residual-angle effects;
- missing multi-hadron tau modes;
- HNL lifetime/width and hadron/parton matching.

For exact ranges and notes, use `per_mass_boundary_impacts.csv`; it is the
authoritative output of this audit.

## Reference basis for broad proxy ranges

- FONLL heavy-flavor framework: https://arxiv.org/abs/1205.6344
- LHCb bottom-hadron fragmentation fractions and kinematic dependence:
  https://arxiv.org/abs/1902.06794
- ALICE charm fragmentation fractions:
  https://arxiv.org/abs/2308.04877
- GeV-scale HNL production and decay treatment:
  https://arxiv.org/abs/1805.08567
