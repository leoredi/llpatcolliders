# `hnl/dev/` — internal scratch, plots, validation scripts

Material that travels with our fork (`myfork = leoredi/llpatcolliders`) and
helps us keep working, but is **never** part of an upstream PR
(`origin = exoticdarksectors/llpatcolliders`).

> If you are reading this in an upstream PR, something went wrong in PR prep —
> flag it. `hnl/dev/` should not be on a PR branch.

## What's here

### `debug_plots.py`
Yield-vs-mass per channel, channel-fraction stacks, and HNL spectra at a fixed
mass. Reads `hnl/output/llp_4vectors/`, writes PNGs to
`hnl/output/debug_plots/`. Run it after `run_all.py` to see whether the
production curves moved.

```bash
python dev/debug_plots.py --channel-fraction Ue
python dev/debug_plots.py --channel-fraction all
```

### `validation_plot.py`
Overlays our `combined` curve against the **MATHUSLA RHN reference 4-vector
files** for each flavor. The factor-of-2 disagreement check on shared channels
(B, D, tau) is exactly this overlay — a clean apples-to-apples on B+D+tau
means we are not double-counting `q + qbar`.

```bash
python dev/validation_plot.py
python dev/validation_plot.py --ref /custom/path/to/MATHUSLA_root
```

It auto-discovers the reference data at
`llpatcolliders_FONLL/vendored/MATHUSLA_LLPfiles_RHN_U{e,mu,tau}/` one or two
levels above `hnl/`; pass `--ref` only if your checkout layout differs. The
expected per-flavor sub-tree is
`All_RHN_*/RHN_*_LLPweight4vector{Bmeson,Dmeson,Tau,WZ}list_mN_*.csv`.

## The one rule

Production code never imports from `hnl/dev/`. Imports the other direction
(dev → production) are fine and expected.

Anything diagnostic, exploratory, or one-off goes here. Don't delete scripts
just because you stopped using them — they cost nothing on disk and may
save you a re-derivation later.

## Environment

Everything in `hnl/` (production and dev alike) runs in the
`llpatcolliders_FONLL` conda env, which has `numpy`, `scipy`, `sympy`,
`mpmath`, `particle`, `numba`, `matplotlib`, and `pytest`:

```bash
conda activate llpatcolliders_FONLL
# or, for one-shot calls:
/Volumes/sandbox/conda/envs/llpatcolliders_FONLL/bin/python dev/<script>.py
```

The base conda env does **not** have the FONLL dependencies and silently
skips the e2e smoke test — always use the env above.

## Workflow — fork vs upstream PR

The fork branch (`hnl-production-mesons-ew`) holds the full state: production
code, tests, `dev/`, scratch — push to it freely. The PR branch is built
fresh from `origin/main` at PR time and only contains the publishable subset.
The two never merge.

### PR-prep recipe

```bash
git fetch origin
git checkout -b pr/<topic> origin/main

# Copy across only the publishable subset (note: hnl/dev/ is absent on purpose)
git checkout hnl-production-mesons-ew -- \
    hnl/README.md \
    hnl/.gitignore \
    hnl/config_mass_grid.py \
    hnl/run_all.py \
    hnl/production/ \
    hnl/tests/ \
    hnl/vendored/

# Safety: hnl/dev/ must NOT be staged
git status hnl/dev   # should print: pathspec did not match any file

# Tests must pass against the same env that produced our current results
/Volumes/sandbox/conda/envs/llpatcolliders_FONLL/bin/python -m pytest hnl/tests/ -q

git commit -m "hnl: <topic>"
git push myfork pr/<topic>

gh pr create \
    --repo exoticdarksectors/llpatcolliders \
    --base main --head leoredi:pr/<topic> \
    --title "hnl: <topic>" --body "<see commit body / link to internal notes>"
```

The `git checkout <branch> -- <paths>` step is what enforces the "minimum
publishable" rule: nothing outside the named paths can leak onto the PR
branch, no matter how messy the dev branch gets.

## Rule of thumb for adding new things

> "Does this file need to exist for an outside reader to reproduce the
> published result?"
>
> - **Yes** → it goes under `hnl/production/`, `hnl/tests/`, or
>   `hnl/vendored/`. It will be in the next PR.
> - **No** → it goes here. It stays on the fork.
