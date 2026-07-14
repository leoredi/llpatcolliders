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
- `bundle/` — exact full-run theory/model variation curves, the compact
  one-source-at-a-time display envelope, and their provenance manifest.

## To re-publish (after a better run)

```bash
cp alp_fermion/tmp/analysis/bc10_sensitivity.csv \
   alp_fermion/data/published/bc10_sensitivity.csv
# then refresh MANIFEST.json (sha, date, csv_sha256, headline_reach)
```

## Current contents

- Source: the 2026-07-15 high-statistics campaign with a 1.2M-event high-pT
  importance production pool, 20,000 full-branching Pythia templates per
  supported mass, 60 decay/reconstruction samples per entering ALP, and the
  decoded arXiv:2501.04525 SensCalc v1.3.3 widths, exclusive branching ratios,
  and exact three-body matrix elements.
- Reach: the deepest lower edge is `1/f = 1.2073e-8 GeV^-1` at 1.23 GeV; the
  final high-mass component closes at approximately 3.293 GeV.
- Topology: unsupported pole rows and finite insensitive rows are retained, so
  plotters preserve the five disconnected runs `0.22-0.52`, `0.56-0.94`,
  `0.98-1.25`, `1.40-1.41`, and `1.47-3.28` GeV.

## Numerical convergence

The high-pT proposal corrects the sparse survival tail that affected the old
upper edge. Independent six-million-event nominal-pool controls at 0.48, 0.50,
0.58, and 0.76 GeV agree with the published upper edge within 2.0%. Controls at
1.25, 1.40, 3.20, and 3.30 GeV independently reproduce the resonance pockets
and the final closure state. In the sensitive heavy-pseudoscalar region the
minimum peak event effective sample size is 1459; the narrow components are
therefore decay-model structure, not a small-simulation artifact. Exact values
and the production configuration are pinned in `MANIFEST.json`.
