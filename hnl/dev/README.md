# `hnl/dev/` - internal development utilities

This directory contains exploratory or diagnostic tools that are useful while
developing the package but are not required to reproduce the production and
analysis chain.

## `debug_plots.py`

This script plots yield versus mass, channel fractions, and HNL spectra from a
completed run:

```bash
conda activate hnl
cd hnl
python dev/debug_plots.py --channel-fraction Ue
python dev/debug_plots.py --channel-fraction all
```

It reads `tmp/runs/<tag>/llp_4vectors/` and writes below
`tmp/runs/<tag>/analysis/debug_plots/`.

Production and analysis code must not import from `hnl/dev/`.

## Environment

Use the committed package environment:

```bash
conda env create -f hnl/environment.yml
conda activate hnl
python hnl/dev/debug_plots.py --help
```

For a one-shot invocation without activation:

```bash
conda run -n hnl python hnl/dev/debug_plots.py --help
```

## PR Preparation

The publishable package includes production, analysis, geometry, static HNL
inputs, tests, and top-level entry points. A branch assembled from another
development branch must include all of them:

```bash
git fetch origin
git checkout -b pr/<topic> origin/main

git checkout hnl-production-mesons-ew -- \
    hnl/README.md \
    hnl/.gitignore \
    hnl/environment.yml \
    hnl/config_mass_grid.py \
    hnl/run_all.py \
    hnl/run_analysis.py \
    hnl/run_full_all.sh \
    hnl/analysis/ \
    hnl/geometry/ \
    hnl/data/ctau/ \
    hnl/data/README.md \
    hnl/production/ \
    hnl/tests/ \
    hnl/vendored/

git status --short
conda activate hnl
python -P -m pytest hnl/tests/ -q
```

Do not include `hnl/tmp/`, generated plots, caches, logs, or local MadGraph
work trees.

## Placement Rule

Files needed by another user to reproduce the result belong in one of:

- `hnl/production/`
- `hnl/analysis/`
- `hnl/geometry/`
- `hnl/data/`
- `hnl/tests/`
- `hnl/vendored/`
- a documented top-level entry point

Exploratory utilities that are not part of the reproducible chain belong in
`hnl/dev/`.
