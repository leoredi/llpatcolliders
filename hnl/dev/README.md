# `hnl/dev/` — internal scratch, plots, validation scripts

This directory holds **developer-only** material that travels with our personal
fork (`myfork = leoredi/llpatcolliders`) but is **never** part of an upstream
PR (`origin = exoticdarksectors/llpatcolliders`).

If you are reading this in an upstream PR, something went wrong in PR prep —
flag it. `hnl/dev/` should not exist on a PR branch.

## What lives here

| File | What it does |
|---|---|
| `debug_plots.py` | Internal diagnostics — yield-vs-mass per channel, channel-fraction stacks, spectra at a fixed `m_N`. Writes PNGs to `hnl/output/debug_plots/`. |
| `validation_plot.py` | Overlays our `combined` curve against MATHUSLA RHN reference 4-vector files. The factor-of-2-disagreement check on shared channels (B, D, tau) lives here. |

Both scripts read from `hnl/output/llp_4vectors/` and write to
`hnl/output/debug_plots/`; both are gitignored by the top-level `output/` rule.

## Rules

- **No production code imports from `hnl/dev/`.** Anything `hnl/production/` or
  `hnl/tests/` depends on must live outside this directory. Imports the other
  way (dev → production) are fine and expected.
- **Add anything diagnostic, exploratory, or one-off here.** Plot scripts,
  validation overlays, scratch notebooks, ad-hoc CSV diff tools. If it would
  embarrass us in front of a reviewer, it goes in `hnl/dev/`.
- **Never delete a script just because you stopped using it.** They are cheap;
  the cost is the seconds to scroll past them. Past-us may have learned
  something we want to recover.

## Workflow: fork vs upstream PR

```
                                 origin/main
                                     │
                  myfork/hnl-production-mesons-ew  ←  daily work, full state
                                                      (production + tests + dev/ + scratch)
                                     │
                       at PR-prep time, fresh branch from origin/main:
                                     │
                          myfork/pr/<topic>         ←  curated subset, no dev/
                                     │
                                     ▼
                         PR → exoticdarksectors:main
```

The dev branch is **not** the source for the PR. The PR branch is built fresh,
file-by-file (or directory-by-directory), so we never accidentally ship
scratch material.

## PR-prep recipe

When the work on the dev branch is ready to be proposed upstream, do this:

```bash
# from the hnl-production-mesons-ew dev branch, in repo root:
git fetch origin
git checkout -b pr/<topic> origin/main

# copy the publishable subset across, excluding hnl/dev/ and hnl/output/
git checkout hnl-production-mesons-ew -- \
    hnl/README.md \
    hnl/.gitignore \
    hnl/config_mass_grid.py \
    hnl/run_all.py \
    hnl/production/ \
    hnl/tests/ \
    hnl/vendored/

# sanity: hnl/dev/ must NOT be staged
git status hnl/dev   # should say "did not match any file"

# run the suite from the env that has the dependencies
/Volumes/sandbox/conda/envs/llpatcolliders_FONLL/bin/python -m pytest hnl/tests/ -q

git commit -m "hnl: <topic>"
git push myfork pr/<topic>
# open PR: leoredi:pr/<topic>  →  exoticdarksectors:main
```

The `git checkout <branch> -- <paths>` step is what enforces the "minimum
publishable" rule: nothing outside the named paths makes it onto the PR
branch.

## What goes in the PR vs what stays here

| Path | PR? | Why |
|---|---|---|
| `hnl/{README.md, .gitignore, config_mass_grid.py, run_all.py}` | yes | entry points and contract |
| `hnl/production/` | yes | the actual physics pipeline |
| `hnl/tests/` | yes | regression coverage |
| `hnl/vendored/` | yes | FONLL tables, HNLCalc, HeavyN UFO model + PROVENANCE |
| `hnl/dev/` | **no** | this directory |
| `hnl/output/` | no | gitignored everywhere |

If the answer to "does this code need to exist for someone to reproduce the
published result?" is **no**, the file belongs in `hnl/dev/`.

## Pushing freely to the fork

The dev branch is yours. Commit small, push often:

```bash
git push myfork hnl-production-mesons-ew
```

There is no review gate on the fork; the gate is on the PR branch.
