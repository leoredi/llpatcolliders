# GRENDEL HNL money-plot bundle

Version-controlled inputs to reproduce `money_plot.pdf` (the GRENDEL HNL
sensitivity projection, single-flavor, P>100 MeV cut) on a fresh clone -- no
`tmp/` run tree required. `analysis/plot_money.py` reads these by default.

Files:
- `../grendel_hnl_sensitivity.csv` -- central exclusion curve (red contour + fill)
- `hnl_band_fonll.csv` -- FONLL heavy-flavor production band (orange, both edges)
- `decay_model_band.csv` -- HNL total-width/lifetime duality band (blue, upper edge)
- `bc_nuisance_band.csv` -- Bc direct-normalization band (teal, lower edge)
- `channel_breakdown_u2min.csv` -- per-channel u2_min composition (provenance; not drawn)

Reproduce (hnl conda env):

    cd hnl
    python -m analysis.plot_money --out-dir tmp/money_repro

Writes `money_plot.pdf`/`.png` + `run_metadata.json` from the tracked CSVs above.
The central curve closes on real refined grid points (3.62-3.70 GeV) where the
peak signal yield crosses N=3 -- no synthetic pinch/extrapolation.
