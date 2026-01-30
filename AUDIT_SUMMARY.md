# HNL Analysis Pipeline: Consolidated Audit Summary

**Date:** 2026-01-26
**Sources:** AUDIT_REPORT.md (Claude Opus 4.5), AUDIT_VERIFICATION.md (Gemini CLI)

---

## Overview

This document consolidates the findings from multiple independent audits of the HNL sensitivity analysis pipeline. All 20 identified issues were verified by at least two reviewers.

**Verification Summary:**
- 20 issues identified
- 19 fully confirmed by all reviewers
- 1 partially confirmed (B4: float precision - lower practical risk than initially stated)

**Resolution Progress (as of 2026-01-29):**
- 8 issues resolved (B1, B3, B6, B7, B8, C7, C9 + B1 was N/A)
- 2 marked low priority (B4, B10)
- Remaining: documentation/validation tasks and minor maintainability items

---

## Physics Biases (10 Issues)

| ID | Issue | Severity | Status | Notes |
|----|-------|----------|--------|-------|
| C1 | W/Z BR sums all |U|² instead of per-flavor | None | N/A | Only single-flavor scans used (100/010/001) |
| C2 | 3-body decays use phase space (meMode=0) | Moderate | Verified | Documented approximation; O(10-20%) shape uncertainty |
| C3 | τ→NX uses representative channel weights | Minor/Moderate | Verified | Kinematic sampling weights, not physical BRs |
| C4 | Decay file mass matching uses nearest neighbor | Minor | **Mitigated** | Doubled mass grid density below 2 GeV; pending regen |
| C5 | Track separation uses midpoint approximation | Minor | Verified | Documented; weak dependence for drainage gallery |
| C6 | Kaon cross-section highly uncertain (~50 mb) | Moderate | Verified | Dominant systematic for m_N < 0.5 GeV |
| C7 | Missing K_L/K_S production channels | Minor | **Resolved** | Added K_L → πℓN; K_S omitted (τ_S/τ_L ≈ 1/570 suppression) |
| C8 | Majorana vs Dirac factor of 2 | None | Verified | Correctly implemented via `--dirac` flag |
| C9 | HNLCalc MC integration uses nsample=100 | Minor | **Resolved** | Increased to nsample=1000 (~3% noise) |
| C10 | Form factors use single-pole dominance | Minor | Verified | Standard approach; O(10-20%) systematic |

---

## Coding Issues (10 Issues)

| ID | Issue | Severity | Status | Notes |
|----|-------|----------|--------|-------|
| B1 | Race condition in summary CSV writes | None | **Resolved** | Removed unused MadGraph summary CSV code |
| B2 | Unsafe eval() in HNLCalc | Minor | Verified | Low risk - expressions are internal |
| B3 | Division-by-zero clamping (λ→1e-9) | Minor | **Resolved** | Now warns when clamping triggers |
| B4 | Float precision in mass filenames | Minor | **Low priority** | Lower practical risk; 2-decimal format mitigates |
| B5 | Hardcoded magic numbers in geometry | Minor | Verified | Maintainability concern; coordinates undocumented |
| B6 | LHE parser defaults missing parent to W+ | Moderate | **Resolved** | Now uses process ID from LHE header to determine W/Z |
| B7 | Inconsistent column naming (boost_gamma/beta_gamma) | Minor | **Resolved** | Migrated fully to `beta_gamma`; legacy handling removed |
| B8 | tqdm in library code | Minor | **Resolved** | Auto-disables in non-TTY; `--no-progress` flag added |
| B9 | Global mutable mesh cache | Minor | Verified | Not thread-safe; single-geometry limitation |
| B10 | Broad exception catch in charge lookup | Minor | **Low priority** | Fails safe to neutral; fallback list covers common cases |

---

## B6 Quantification (From Production Data)

Analysis of 279 MadGraph EW output files (27.9M events):

| Parent PDG | Count | Percentage |
|------------|-------|------------|
| W+ (24) | 11,590,178 | 41.54% |
| W- (-24) | 8,795,681 | 31.53% |
| Z (23) | 7,514,141 | 26.93% |
| Unknown/0 | 0 | 0.00% |

**Conclusion:** No parent_pdg=0 entries observed in existing outputs. The W+ fallback was rarely triggered in practice.

**Resolution (2026-01-29):** `lhe_to_csv.py` now parses the LHE header (`<MGProcCard>` and `<init>` blocks) to build a process ID → parent PDG mapping. When the parent boson is absent from the particle list (off-shell), the correct W+/W-/Z is inferred from the event's `idprup` field. No more blind defaulting to W+.

---

## Prioritized Recommendations

### High Priority

1. **Document kaon cross-section uncertainty** - Dominant systematic for m_N < 0.5 GeV
2. ~~**Consider adding K_L/K_S channels** - Captures full kaon-regime sensitivity~~ ✓ Added K_L
3. ~~**Add file locking to CSV writes** - For parallel-safe production runs~~ ✓ Removed unused summary CSV
4. ~~**Increase HNLCalc nsample** - From 100 to 1000 for better BR precision~~ ✓ Done

### Medium Priority

5. **Validate phase-space vs matrix-element kinematics** - Compare to known distributions
6. **Create geometry abstraction layer** - Enable different detector configurations
7. ~~**Consolidate column naming** - Migrate fully to `beta_gamma`~~ ✓ Done
8. **Add integration tests** - Compare pipeline output to benchmark results

### Low Priority

9. **Move magic numbers to config files** - Coordinates, transformations
10. **Add logging framework** - Replace scattered print statements
11. **Refactor HNLCalc** - Replace string-based eval() with symbolic math

---

---

## Files Reviewed

**Production:**
- `production/pythia_production/main_hnl_production.cc`
- `production/madgraph_production/scripts/run_hnl_scan.py`
- `production/madgraph_production/scripts/lhe_to_csv.py`

**Analysis:**
- `analysis_pbc/config/production_xsecs.py`
- `analysis_pbc/models/hnl_model_hnlcalc.py`
- `analysis_pbc/geometry/per_parent_efficiency.py`
- `analysis_pbc/decay/decay_detector.py`
- `analysis_pbc/decay/rhn_decay_library.py`
- `analysis_pbc/limits/run.py`
- `analysis_pbc/limits/expected_signal.py`
- `analysis_pbc/limits/combine_production_channels.py`
- `analysis_pbc/HNLCalc/HNLCalc.py`

**Configuration:**
- `config_mass_grid.py`

**Visualization:**
- `money_plot/plot_money_island.py`

---

*End of Consolidated Audit Summary*
