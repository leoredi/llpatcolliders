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
  timing-χ² β=p/E fix — PR #15 `f3d69bc`; reuses the same production 4-vectors as
  the pre-fix run). The conservative P>600 MeV variant is `exact_600`.
- Reach: `|U_e|^2` 5.9e-9 @0.40 GeV, `|U_mu|^2` 8.4e-9 @0.37 GeV,
  `|U_tau|^2` 3.0e-7 @2.8 GeV; window 0.2-3.6 GeV (closure is the cτ∝1/m^5
  lifetime law, not a production cutoff — see `MANIFEST.json`).
- **Timing β-fix note:** the previous publish computed the 4-hit timing χ² with
  β=1 for every daughter; it now uses each track's true β=p/E, which imposes a
  species-dependent effective momentum floor above the 100 MeV cut. No mass points
  are lost and the closure is unchanged; `u2_min` degrades modestly (median
  +3.7/+4.5/+9% for Ue/Umu/Utau, worst +33% at low mass). The FONLL/decay-model/Bc
  ribbons in `bundle/` are **not re-derived** — they are dex (log₁₀U²) half-widths,
  which are β-fix-invariant to second order, and `plot_money` re-anchors them onto
  this central automatically.

## Consumer

`../curves_PBC` reads this file (default in `digitize/paths.py`, overridable via
`CURVES_PBC_GRENDEL_CSV`) and overlays the PBC BC7 contours + the HNLimits
community compilation to produce the final comparison plots. **This repo only
produces the GRENDEL curve; the comparison/final figures are made there.**
