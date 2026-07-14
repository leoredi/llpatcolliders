# GRENDEL BC10 variation-envelope bundle

Version-controlled theory/model variation inputs for the canonical BC10
sensitivity curve in `../bc10_sensitivity.csv`.

Files:

- `bc10_uncertainty_variations.csv` - all 114 exact full-run contours: central,
  six scale choices, 100 NNPDF replicas, two bottom-mass choices, three gluon
  decay surrogates, and two `C_bs` scheme choices.
- `bc10_single_source_variation_envelope.csv` - the 99-row plot input with the
  central contour, one-source-at-a-time display envelope, open/missing flags,
  named source extrema, PDF percentiles, and PDF log-coupling standard
  deviation for audit.
- `UNCERTAINTY_MANIFEST.json` - exact variation registry, input and output
  checksums, code/template provenance, combination prescription, and headline
  shifts.

The display object is explicitly a
`single_source_variation_envelope`, not a confidence interval. Scale,
bottom-mass, `C_bs`, and `a -> gg` surrogate alternatives use their named
pointwise extrema. The 100 NNPDF replicas use their 16th and 84th percentiles
in `log10(1/f)`; their raw extrema do not define the display envelope. The
overall edge is the outermost of these source intervals with one source varied
at a time. Sources are not added in quadrature.

Every FONLL and `C_bs` variation uses a fresh 600,000-parent production pool
and full geometry, reconstruction, and sensitivity scan. The `u`, `d`, and
`s` gluon-surrogate templates use exact central production-vector reuse because
their production physics is identical, followed by independent downstream
simulation. No importance reweighting or detector-outcome reuse is used.

Central insensitive rows and the excluded eta/eta-prime pole windows remain
explicit NaN gaps. Heavy resumable run trees stay in external scratch. After a
non-central run validates and records its raw vector/geometry hashes, those
regenerable trees are compacted; the sensitivity CSV, markers, provenance, and
logs remain. Central raw output is retained for exact gluon-variant reuse.

Reproduce in the `llpatcolliders_FONLL` environment using the commands in
`alp_fermion/UNCERTAINTY.md`, then collect with:

```bash
python -m alp_fermion.combine_uncertainty_band \
  --scratch-root /scratch/bc10_uncertainty \
  --grid-dir /path/to/fonll-nnpdf40/output
```
