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

- Source run: `central_newgrids_20260623/analysis_exact_100` (P>100 MeV, the
  deepest/latest exact run). The conservative P>600 MeV variant is `exact_600`.
- Reach: `|U_e|^2` 5.0e-9 @0.40 GeV, `|U_mu|^2` 6.6e-9 @0.34 GeV,
  `|U_tau|^2` 2.8e-7 @2.8 GeV; window 0.2-3.6 GeV (closure is the cτ∝1/m^5
  lifetime law, not a production cutoff — see `MANIFEST.json`).

## Consumer

`../curves_PBC` reads this file (default in `digitize/paths.py`, overridable via
`CURVES_PBC_GRENDEL_CSV`) and overlays the PBC BC7 contours + the HNLimits
community compilation to produce the final comparison plots. **This repo only
produces the GRENDEL curve; the comparison/final figures are made there.**
