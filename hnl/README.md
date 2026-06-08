# `hnl/` - Heavy Neutral Lepton production and sensitivity

This package runs the full GRENDEL HNL chain for proton-proton collisions at
14 TeV:

1. produce HNL four-vectors from mesons, kaons, induced taus, prompt taus,
   and electroweak W/Z processes;
2. combine the production channels for each flavor and mass;
3. ray-cast the HNL trajectories through the GRENDEL detector geometry;
4. scan the active-flavor mixing and produce sensitivity tables and plots.

The supported flavor hypotheses are `Ue`, `Umu`, and `Utau`. The committed
mass grid contains 116 points from 0.2 to 10 GeV.

## Layout

```text
hnl/
|-- README.md
|-- environment.yml
|-- config_mass_grid.py
|-- run_all.py                 production orchestrator
|-- run_analysis.py            sensitivity entry point
|-- run_full_all.sh            complete three-flavor workflow
|-- analysis/                  decay probability, acceptance, contours, plots
|-- geometry/
|   `-- grendel_geometry.py    GRENDEL fiducial mesh and ray casting
|-- data/
|   `-- ctau/                  committed lifetime and visible-BR tables
|-- production/
|   |-- combine_channels.py
|   |-- decay_engine/          meson, kaon, and tau decays
|   |-- fonll/                 heavy-meson spectrum parser and sampler
|   `-- madgraph/
|       |-- run_tau_production.py
|       |-- run_wz_production.py
|       `-- run_wz_sharded.py
|-- tests/
|-- vendored/
|   |-- HNLCalc/
|   |-- SM_HeavyN_CKM_AllMasses_LO/
|   `-- FONLL tables
`-- tmp/                       generated runs, caches, logs; git-ignored
```

## Physics Chain

Production includes:

- charged and semileptonic decays of `K`, `D`, `Ds`, `B`, `Bs`, and `Bc`;
- induced taus from `Ds -> tau nu` and `B+ -> tau nu`;
- prompt taus from `W -> tau nu` and `gamma*/Z -> tau tau`;
- `W/Z -> ell N` through MadGraph.

The meson spectra come from the committed 14 TeV FONLL tables. HNLCalc
provides production branching ratios, differential three-body rates, HNL
lifetimes, visible branching fractions, and tau-decay rates.

For a production row generated at unit mixing, the downstream yield is

```text
N_signal = L_int * sum_i [
    weight_i * U_alpha^2 * P_decay_i(U_alpha^2) * acceptance_i
]
```

The analysis implements this scan at `3000 fb^-1` and uses
`N_signal >= 3` as its default threshold.

### Decay Kinematics

Two-body decays use exact rest-frame kinematics followed by a Lorentz boost.
The polarized two-body tau modes use the configured longitudinal analyzing
power.

Meson three-body decays are sampled with
`decay_3body_weighted_dq2dE`, using the HNLCalc differential
`dBR/(dq2 dE)` expression. Leptonic tau three-body decays use
`decay_3body_weighted_dE` with the HNLCalc `dBR/dE` expression. These are not
flat phase-space samplers. The remaining approximation is the treatment of
angular and spin correlations that are not contained in those one- or
two-dimensional differential rates.

### Stored Weights

All production CSVs are headerless:

```text
weight, E, px, py, pz
```

`weight` is in pb and factors out the active mixing, so it corresponds to
`U_alpha^2 = 1`. Four-momenta are in GeV in the CMS lab frame.

Direct meson weights use

```text
2 * sigma_FONLL * fragmentation_fraction * BR(parent -> N + X) / N_sampled
```

where the factor of two converts the FONLL quark-plus-antiquark convention to
the total rate. `Bc` uses its separately configured cross section.

Induced-tau weights additionally include `BR(parent -> tau nu)`. Prompt-tau
and W/Z rows use their MadGraph event weights and the configured electroweak
K-factor. The kaon channel uses the approximate inclusive kaon flux documented
in `production/constants.py`.

## Setup

Create and activate the committed environment:

```bash
conda env create -f hnl/environment.yml
conda activate hnl
```

The direct runtime dependencies include `numpy`, `pandas`, `scipy`, `sympy`,
`mpmath`, `particle`, `matplotlib`, `tqdm`, `trimesh`, `rtree`, `numba`, and
`cycler`. `pytest` is included for tests. The full runner checks these imports
before starting expensive production.

MadGraph5_aMC@NLO v3.6.6 is required for prompt-tau and W/Z production. Its
executable is resolved in this order:

1. `$HNL_MG5_EXE`;
2. `hnl/vendored/MG5_aMC_v3_6_6/bin/mg5_aMC`;
3. the shared project-level install;
4. the legacy sibling-checkout install.

LHAPDF is supplied by the conda environment. The PDF data directory can be
overridden with `$HNL_LHAPDF_DATA`.

## Full Run

From the repository root:

```bash
conda activate hnl
cd hnl
./run_full_all.sh
```

The script performs the complete three-flavor chain. It first runs non-W/Z
production, then runs W/Z production in independent mass shards, combines all
channels, and runs the analysis.

By default it uses:

- 6 non-W/Z production workers;
- 12 cores for prompt-tau MadGraph generation;
- 8 simultaneous W/Z shards with one MG5 core each;
- 3 analysis workers.

Analysis workers load large combined CSVs and geometry arrays. On an 18 GB
development machine, 12 analysis workers caused severe memory pressure and a
userspace watchdog restart. Three workers completed reliably. This limit is
about memory, not CPU utilization.

The runner creates a unique tag such as `full_20260608_231500`, preventing a
new invocation from silently overwriting an older run. To resume a run, reuse
its tag:

```bash
HNL_RUN_TAG=full_20260608_231500 ./run_full_all.sh
```

W/Z production passes `--skip-existing`, so non-empty completed mass points
are retained. Other useful overrides are:

```bash
HNL_PYTHON=/path/to/python \
HNL_PRODUCTION_WORKERS=6 \
HNL_PROMPT_TAU_CORES=12 \
HNL_WZ_JOBS=8 \
HNL_WZ_CORES_PER_JOB=1 \
HNL_ANALYSIS_WORKERS=3 \
./run_full_all.sh
```

The target W/Z CPU use is `HNL_WZ_JOBS * HNL_WZ_CORES_PER_JOB`. Avoid setting
both values high enough to oversubscribe the machine.

Generated files are written below:

```text
tmp/runs/<tag>/llp_4vectors/
tmp/runs/<tag>/analysis/
tmp/cache/
tmp/wz_shards/
```

Top-level stage logs are written as `tmp/<tag>_*.log`.

## Manual Operation

### Production without W/Z

```bash
python -u run_all.py \
  --flavor Ue Umu Utau \
  --no-wz \
  --skip-combine \
  --workers 6 \
  --prompt-tau-nb-core 12
```

`run_all.py` parallelizes channel/flavor jobs. Its built-in W/Z path processes
the masses serially inside each flavor job; the full runner deliberately uses
the sharded W/Z driver instead.

### Sharded W/Z Production

```bash
python -u production/madgraph/run_wz_sharded.py \
  --flavor Umu \
  --jobs 8 \
  --nb-core 1 \
  --skip-existing
```

Run this once per flavor. Every shard has an isolated MadGraph work directory,
which avoids process-tree collisions.

### Combine

```bash
python -m production.combine_channels --flavor Ue Umu Utau
```

### Analysis

```bash
python -u run_analysis.py --flavor Ue Umu Utau --workers 3
```

The analysis writes:

```text
tmp/runs/<tag>/analysis/hnl_sensitivity.csv
tmp/runs/<tag>/analysis/hnl_exclusion.pdf
tmp/runs/<tag>/analysis/hnl_exclusion.png
tmp/runs/<tag>/analysis/scan_status.json
tmp/runs/<tag>/analysis/run_metadata.json
```

If an exclusion boundary is not crossed inside the configured `U^2` scan, the
CSV stores a `NaN` boundary and sets `u2_min_open` or `u2_max_open`. The plot
fills to the corresponding axis edge, omits a false closing line, and marks
the open direction.

Use `--plot-only` to regenerate the plot from an existing sensitivity CSV and
`--force-geometry` to rebuild cached ray intersections.

## Paths and Files

Path defaults are centralized in `production/paths.py`. Supported overrides
include:

```text
HNL_RUN_TAG
HNL_RUN_DIR
HNL_TMP_DIR
HNL_MG5_WORK_DIR
HNL_TAU_POOL_CSV
HNL_LLP_VECTORS_DIR
HNL_ANALYSIS_DIR
```

Mass filenames use three decimal places with `p` replacing the decimal point:
`1.025 GeV -> mN_1p025.csv`. Use
`config_mass_grid.parse_mass_from_filename` rather than parsing labels by
hand.

Zero-byte channel files are intentional threshold sentinels. They mean that
the channel is closed at that mass. Non-empty files near a threshold can have
fewer rows than the requested pool because only a subset of parent species is
kinematically allowed.

## Tests

From the repository root:

```bash
conda activate hnl
python -P -m pytest hnl/tests/ -q
```

The tests cover FONLL parsing and sampling, mass labels, channel combination,
two- and three-body kinematics, tau and kaon production, W/Z CSV conversion,
open contour handling, and a meson-production smoke path. `-P` prevents any
stale ignored MG5 parser cache named `py.py` from shadowing pytest's
compatibility module.

## Physics Inputs and Limitations

- The default heavy-meson backend is the committed NNPDF4.0 NLO FONLL grid.
  `HNL_FONLL_SET=cteq66_legacy` selects the older alternate backend.
- The FONLL tables stop at `pT = 50 GeV` and contain one central scale choice.
  PDF and scale uncertainties are not propagated.
- Charm and bottom species share one heavy-flavor shape per table; species
  fractions are applied in the event weights.
- Bottom fragmentation fractions are static and omit explicit heavy-baryon
  production channels.
- The kaon flux is a parametrized soft-QCD approximation and dominates the
  normalization uncertainty at the lightest masses.
- Meson and tau three-body energy distributions are HNLCalc-weighted, but
  complete multidimensional matrix-element spin correlations are not modeled.
- The W/Z K-factor is a constant approximation.
- The analysis uses a zero-background three-event threshold. It does not
  include detector backgrounds or systematic uncertainties.

The FONLL and HNLCalc source details, versions, and local modifications are
recorded in `vendored/PROVENANCE.md`.
