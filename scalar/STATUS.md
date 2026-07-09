# BC4 (light dark scalar) — status: GRENDEL curve FINAL

_Last updated 2026-07-09. Worktree `aaaPHYSICSaaa/bc4-scalar` (top level, sibling of `bc10-alp`), branch `bc4-scalar`, package `scalar/`._

**2026-07-09 β-fix republish.** The PR #15 timing-χ² fix (daughter β = p/E
instead of β = 1) is now propagated through the BC4 acceptance (commit
7f792dc) and the full 82-mass grid was rerun (`tmp/rerun_betafix2.log`,
finished 00:58). The canonical curve now lives at
`scalar/data/published/bc4_island.csv` (+ MANIFEST): island spans
**0.22–3.70 GeV**, deepest **sin²θ = 7.4e-12 @ 0.975 GeV**. Versus the
2026-07-02 pre-fix run: the lower edge is ~4% shallower (consistent with the
HNL β-fix median) and the 3.80 GeV point drops below N = 3, closing the span
at 3.70 GeV. The numbers below in "Results" refer to the superseded pre-fix
run and are kept for the exclusive-vs-inclusive comparison.

**2026-07-08 close-out.** The deliverable is the GRENDEL exclusion curve; the
former TODOs are resolved as follows (details in the TODO section):
- **gg→S: out of scope for BC4.** The competitor projections define the
  benchmark as B-meson production only — Evans (arXiv:1708.08503, the MATHUSLA
  scalar treatment: inclusive BR(B→X_s S) ≈ 6.2 (1−m_S²/m_B²)² sin²θ, kaons a
  "small correction", no gluon fusion) and the CODEX-b physics case
  (arXiv:1911.00481: b→s penguin dominant; direct Higgs production only enters
  the non-minimal quartic scenario = BC5). Our inclusive b→X_s S with Winkler
  dispersive widths meets or exceeds that standard; the ~3.8 GeV meson ceiling
  is the benchmark's own ceiling.
- **Competitor overlays: dropped by decision** (GRENDEL curve only).
- **Committed** as 6dc4a74 on `bc4-scalar`, pushed to leoredi/llpatcolliders.
The curve in `scalar/tmp/bc4_exclusion.{png,pdf}` / `bc4_island.csv` (full run
2026-07-02) is the final BC4 result.

## What BC4 is
Coupling-controlled, model-complete GRENDEL sensitivity for the PBC **BC4** benchmark
(Higgs-portal dark scalar S mixing with the SM Higgs, mixing angle sin²θ). A single
coupling sets production, lifetime cτ, and all decay BRs simultaneously → a **closed
island** in (m_S, sin²θ), requiring N_signal ≥ 3 at 3000 fb⁻¹, background-free.
Reuses the `higgs/` shared reco (`grendel_geometry` + `reco_common`) and the `hnl/`
FONLL b-hadron sampler — never copied, imported.

## Current physics (both refinements applied 2026-07-02)
1. **Decay widths — Winkler dispersive.** Hadronic widths from digitized Winkler
   (arXiv:1809.01876) Fig. 4 (`scalar/data/winkler_widths.csv`, log-log interp in
   `model._winkler_width`). `M_SPECTATOR = 2.0` gates the below-2 GeV (ππ/KK/4π)
   vs above-2 GeV (gg/ss/cc) channels to avoid a double-count seam. This moved the
   hadronic pinch from a spurious ~1.7 GeV LO-ChPT artifact onto the physical
   **f₀(980)/2m_K peak at ~0.98 GeV**.
2. **Production — inclusive `b → X_s S`.** Switched from exclusive `B → K S` to the
   inclusive spectator rate (`model.br_B_to_Xs_S`, Winkler eq. A7, ~5.3·sin²θ),
   summed over **B⁺, B⁰, B_s** (the HNL inclusive bottom set; b-baryons omitted).
   Rationale: GRENDEL reconstructs only the S vertex — the prompt X_s system is
   invisible, so we sum over it inclusively, exactly as for HNLs. S-spectrum
   kinematics unchanged (recoil = m_K, spectrum <2% sensitive to it at a 5.3 GeV
   parent); only the normalization moves.

## Results
| | exclusive B→K S | **inclusive b→X_s S (current)** |
|---|---|---|
| deepest reach (sin²θ) | 4.7e-11 @ 2.8 GeV | **7.1e-12 @ 0.975 GeV** |
| deepening at matched mass | — | ~3.5× (√10 from rate) |
| sensitive mass span | 0.22 – 3.5 GeV | **0.22 – 3.8 GeV** |
| f₀(980) region (~1 GeV) | island closed (gap) | **open, and deepest** |

The f₀(980) pinch flipped from a closed gap into the deepest point: the width peak
gives the best lower edge (sin²θ_min ∝ 1/√Γ_total in the long-lifetime regime), and
inclusive rate finally lifts peak_N there (4 → 34) above threshold to exploit it.
Above ~3.8 GeV: rate-starved (peak_N → 1 as m_S → m_B) — meson ceiling is real.

Outputs: `scalar/tmp/bc4_exclusion.{png,pdf}`, `scalar/tmp/bc4_island.csv`.
Superseded (exclusive) outputs preserved in `scalar/tmp/pre_inclusive/`.

## TODO — all resolved 2026-07-08
- [x] **gg→S direct production** — RESOLVED: out of scope. The BC4 standard set
      by the competitor studies (Evans arXiv:1708.08503 / MATHUSLA; CODEX-b
      arXiv:1911.00481) is B-meson production only; none include gg→S. BC4 is
      held as the ≤3.8 GeV inclusive-meson benchmark.
- [x] **Competitor / existing-bound overlays** — DROPPED by decision: the
      GRENDEL curve is the deliverable. (`plot_exclusion` still picks up
      `scalar/data/competitors/*.csv` if that is ever revisited.)
- [x] **Commit / PR decision** — committed as 6dc4a74 on `bc4-scalar`, pushed
      to leoredi/llpatcolliders. Upstream PR deferred until the paper decides
      which benchmarks it carries.

## Related (portal-wide)
- **BC10 fermiophilic ALP** — worktree `../bc10-alp`, package `alp_fermion/`. Built,
  closed island, reach 1/f ≈ 2e-8 GeV⁻¹; rate-starved above ~2.5 GeV. NOT yet given
  the inclusive-production treatment applied to BC4 here — candidate for the same fix.
- Done & leading already on `main`: BC5 (higgs portal h→SS) and BC6/7/8 (HNL).
