# GRENDEL BC4 uncertainty bundle

Version-controlled theory uncertainty inputs for the canonical BC4 sensitivity
curve in `../bc4_island.csv`.

Files:

- `bc4_uncertainty_curves.csv` — all 110 independently simulated compact curves.
- `bc4_uncertainty_band.csv` — central contour, separate FONLL ribbon, alternate
  decay-model contour, combined ribbon, open-edge flags, and every component in
  dex for both contour edges.
- `bc4_uncertainty_variations.json` — pinned FONLL grid checksums, campaign
  configuration, combination prescription, output checksum, and limitations.

The FONLL contribution uses all 109 coherent bottom grids from the NNPDF4.0 NLO
campaign: central, the six non-central points of the standard seven-point scale
set, 100 Monte Carlo PDF replicas, and `m_b = 4.5, 5.0 GeV`. Every variation
uses a fresh 200,000-event parent pool, producing 600,000 scalar events per
mass after the three B-species contributions. Scalar kinematics, geometry,
decay/reconstruction MC, and the sensitivity scan are all independent between
variations; no importance reweighting or detector-outcome reuse is used.

Combination is in `log10(sin^2 theta)`: asymmetric scale envelope, sample
standard deviation over PDF replicas, and maximum absolute bottom-mass shift.
The decay-model contribution comes from a separate fresh central-FONLL
production and full detector simulation using analytic LO-ChPT widths below
2 GeV plus spectator widths above 2 GeV, compared with the matched/dispersive
Winkler central run. It is an alternate-model envelope, not a Gaussian error.
Independent FONLL and model components are added in quadrature per direction
and rebased onto the canonical central curve.

Reproduce from the repository root in the `llpatcolliders_FONLL` environment:

```bash
python -m scalar.uncertainty_band run --workers 2 \
  --grid-dir /Volumes/sandbox/projects/aaaPHYSICSaaa/shared/NNPDF40/fonll-local/output \
  --scratch-dir /Volumes/GRENDEL/extra_space/bc4_uncertainty
python -m scalar.uncertainty_band status
python -m scalar.uncertainty_band collect
```

The heavy resumable run tree stays in external scratch and is not published.
Each vector, result, and aggregate curve is written atomically and pinned by
checksums. This band does not include FONLL alpha-s companions, detector
response systematics, or background uncertainty.
