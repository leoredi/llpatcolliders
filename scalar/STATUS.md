# BC4 (light dark scalar) — status & todo

_Last updated 2026-07-03. Worktree `../bc4-scalar`, branch `bc4-scalar`, package `scalar/`._

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

## TODO
- [ ] **gg→S direct production** (only way past the ~3.8 GeV meson ceiling).
      Non-mesonic, θ²-suppressed, Higgs-like. Needs a real σ(gg→S) source — do NOT
      fabricate (see `EXTERNAL_INPUTS_NEEDED.md`). Decision pending: build it, or
      hold BC4 as the ≤3.8 GeV inclusive-meson benchmark.
- [ ] **Competitor / existing-bound overlays** for the plot: CHARM, LHCb (B→Kμμ),
      MATHUSLA, CODEX-b, ANUBIS, SHiP. Digitized CSVs → `scalar/data/competitors/`
      (`plot_exclusion` already skips them gracefully when absent).
- [ ] **Commit / PR decision** — branch `bc4-scalar` is uncommitted; outputs in `tmp/`.
      Port onto `main` as a purely-additive package mirroring `hnl/` when ready.

## Related (portal-wide)
- **BC10 fermiophilic ALP** — worktree `../bc10-alp`, package `alp_fermion/`. Built,
  closed island, reach 1/f ≈ 2e-8 GeV⁻¹; rate-starved above ~2.5 GeV. NOT yet given
  the inclusive-production treatment applied to BC4 here — candidate for the same fix.
- Done & leading already on `main`: BC5 (higgs portal h→SS) and BC6/7/8 (HNL).
