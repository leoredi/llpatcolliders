# Published GRENDEL HNL curve (canonical handoff)

This is the **stable, canonical** GRENDEL HNL sensitivity curve that downstream
repositories consume. It is a *copy* of one chosen analysis run, promoted to a
fixed path so consumers never hardcode a per-run `tmp/runs/<tag>/...` directory
(which changes every run and gets garbage-collected).

## Files

- `grendel_hnl_sensitivity.csv` — the curve. Same schema as the analysis output
  (`hnl/analysis/run_sensitivity.py`); the columns the consumer requires are
  `mass_GeV, flavor, u2_min, u2_max` (plus `has_sensitivity`, `u2_max_open`).
- `MANIFEST.json` — provenance: source run, git sha, cut, luminosity, threshold,
  per-flavor headline reach, and the high-mass-closure physics note.

## To re-publish (after a better run)

Copy the chosen run's `hnl_sensitivity.csv` here and update `MANIFEST.json`:

```bash
cp tmp/runs/<tag>/analysis/hnl_sensitivity.csv \
   data/published/grendel_hnl_sensitivity.csv
# then refresh MANIFEST.json (run tag, sha, date, csv_sha256, headline_reach)
```

## Current contents

- Source run: `central_newgrids_20260623/analysis_exact_100_betafix` (P>100 MeV,
  timing-chi2 `beta=p/E` fix -- implementation `79c7337`, canonical publication
  `fad3a8b`; reuses the same production 4-vectors as the pre-fix run). The
  conservative P>600 MeV variant is `exact_600`.
- Reach: `|U_e|^2` 5.9e-9 @0.40 GeV, `|U_mu|^2` 8.4e-9 @0.37 GeV,
  `|U_tau|^2` 3.0e-7 @2.8 GeV; window 0.2-3.6 GeV (closure is the cτ∝1/m^5
  lifetime law, not a production cutoff — see `MANIFEST.json`).
- **Timing beta-fix note:** the previous publish computed the 4-hit timing chi2
  with `beta=1` for every daughter; it now uses each track's true `beta=p/E`, which imposes a
  species-dependent effective momentum floor above the 100 MeV cut. No mass points
  are lost and the closure is unchanged; `u2_min` degrades modestly (median
  +3.7/+4.5/+9% for Ue/Umu/Utau, worst +33% at low mass).
- The tracked diagnostics in `bundle/` were fully re-derived after this fix:
  111 exact-hit FONLL variations plus refreshed decay-model, Bc-normalization,
  and channel-composition scans.  Their campaign and post-processing revisions
  are `6482484` and `2ab1f98`; see `bundle/MANIFEST.json`.  These are separate
  theory/model diagnostics, not confidence bands on the central paper contour.

## Consumer

`shared/curves_PBC` reads this file through `digitize/paths.py`; override the
default with `HNL_GRENDEL_CSV` (`CURVES_PBC_GRENDEL_CSV` is a legacy alias).
It overlays the PBC BC7 contours and the HNLimits community compilation to
produce the final comparison plots. **This repo only produces the GRENDEL
curve; the comparison/final figures are made there.**
