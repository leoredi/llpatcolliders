# Published GRENDEL BC4 curve (canonical handoff)

This is the **stable, canonical** GRENDEL BC4 (Higgs-portal dark scalar)
sensitivity curve, following the same publish pattern as
`hnl/data/published/`: a *copy* of one chosen analysis run promoted to a
fixed, committed path so consumers (the paper, comparison plots) never read a
git-ignored `tmp/` file.

## Files

- `bc4_island.csv` — the curve. Same schema as the analysis output
  (`scalar/run_sensitivity.py`); the columns a consumer requires are
  `mass_GeV, u2_min, u2_max` (plus `has_sensitivity`), where `u2` is
  `sin^2 theta`.
- `MANIFEST.json` — provenance: source run, git sha, csv sha256, inputs,
  luminosity, threshold, headline reach, and the closure note.

## To re-publish (after a better run)

```bash
cp scalar/tmp/bc4_island.csv scalar/data/published/bc4_island.csv
# then refresh MANIFEST.json (sha, date, csv_sha256, headline_reach)
```

## Current contents

- Source: the 2026-07-15 high-statistics campaign with a 1.2M-scalar high-pT
  importance production pool per mass, 100 decay/reconstruction samples per
  entering scalar, and the timing-chi2 daughter `beta = p/E` correction.
- Reach: the sensitive grid spans `m_S = 0.14-3.70 GeV`; the upper boundary is
  open above the scan ceiling for the electron-only `0.14-0.20 GeV` extension.
  The deepest lower edge is `sin^2 theta = 7.3988e-12` at 0.975 GeV, near the
  f0(980) enhancement. Log-yield interpolation places the final closure at
  approximately 3.798 GeV.

## Numerical convergence

The importance proposal gives at least 2249 effective events on every finite
sensitive lower edge and at least 2049 at every sensitive peak. Independent
six-million-event controls at 3.75, 3.80, and 3.85 GeV give peak yields of
3.925, 2.858, and 2.069, respectively; their interpolated closure is 3.792 GeV.
Exact values, seeds, hashes, and the production configuration are pinned in
`MANIFEST.json`.
