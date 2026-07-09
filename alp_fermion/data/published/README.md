# Published GRENDEL BC10 curve (canonical handoff)

This is the **stable, canonical** GRENDEL BC10 (fermiophilic ALP) sensitivity
curve, following the same publish pattern as `hnl/data/published/`: a *copy*
of one chosen analysis run promoted to a fixed, committed path so consumers
(the paper, comparison plots) never read a git-ignored `tmp/` file.

## Files

- `bc10_sensitivity.csv` — the curve. Same schema as the analysis output
  (`alp_fermion/sensitivity.py`); the columns a consumer requires are
  `mass_GeV, invf_min, invf_max` (plus `has_sensitivity`). Couplings are
  `1/f` in GeV^-1, BNT convention (`1/f = 2/f_GKOZ`).
- `MANIFEST.json` — provenance: source run, git sha, csv sha256, inputs,
  luminosity, threshold, headline reach, and the closure note.

## To re-publish (after a better run)

```bash
cp alp_fermion/tmp/analysis/bc10_sensitivity.csv \
   alp_fermion/data/published/bc10_sensitivity.csv
# then refresh MANIFEST.json (sha, date, csv_sha256, headline_reach)
```

## Current contents

- Source: the 2026-07-08 sensitivity rerun on the PR #15-synced `hnl/`
  subtree (timing-chi2 beta = p/E fix included).
- Reach: island spans `m_a` 0.22–2.5 GeV; deepest `1/f = 1.1e-8 GeV^-1`
  at 1.4 GeV.
