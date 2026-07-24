# Existing BC1 exclusions

`existing_exclusion_senscalc.npz` is a lossless compressed export of the
`Baseline` entry in SensCalc v1.3.3
`contours/DP/Constraints-DP.json`, revision
`0bca050633aae16e148d47f21840fa07ff4b8724` (BSD-3-Clause).

The source coordinates are dark-photon mass in GeV and kinetic mixing
`epsilon`. `dark_photon.plot_exclusion` squares the ordinate when overlaying
the envelope on the GRENDEL `epsilon^2` result. Regenerate with:

```bash
python -m dark_photon.generator.export_existing_constraints
```

The AxionLimits dark-photon collection was reviewed but not used: its main
tables address ultralight dark-photon dark matter in eV and do not share the
visible-decay BC1 assumptions of this plot.
