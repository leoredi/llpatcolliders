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

- Source: the 2026-07-09 full-grid rerun with the timing-chi2 beta = p/E fix
  propagated through the BC4 acceptance (supersedes the 2026-07-02 pre-fix
  run).
- Reach: island spans `m_S` 0.22–3.70 GeV; deepest `sin^2 theta = 7.4e-12`
  at 0.975 GeV (the f0(980) width peak).
