# HNL Analysis Pipeline: Verified Relevant Audit Summary

**Date:** 2026-02-06
**Verification:** This document has been generated after a thorough, independent verification of each original audit point against the current codebase. Only issues that are confirmed to be still present are included.

---

## Overview

This document summarizes the audit findings that remain relevant and unresolved in the HNL sensitivity analysis pipeline. Many issues from previous audits have been confirmed as **resolved** and are omitted here for clarity. The items listed below are either known, documented approximations in the physics model or outstanding code quality concerns.

---

## Relevant Physics Model Approximations & Uncertainties

These are not bugs, but documented, intentional choices in the physics modeling that introduce systematic uncertainties.

| ID | Issue | Severity | Status | Verification Notes |
|----|-------|----------|--------|-------|
| C1 | **W/Z Boson BR Formula Uses Total \|U\|²** | Documentation | Verified | The code in `hnl_model_hnlcalc.py` correctly calculates the *inclusive* W/Z branching ratio by summing all `|U_l|²`. This is correct for single-flavor benchmarks but is a point of clarification for users. |
| C2 | **Phase-Space Kinematics in 3-Body Decays** | Moderate | Verified | `main_hnl_production.cc` intentionally uses `meMode=0` (pure phase space) for 3-body semileptonic decays. This is a documented approximation known to introduce O(10-20%) shape uncertainty. |
| C3 | **τ → N X Uses Representative Channels** | Minor/Moderate | Verified | `main_hnl_production.cc` uses a weighted mixture of representative decay channels (π, ρ, 3π, etc.) for `τ→NX` kinematics. These weights are not physical Branching Ratios, which are applied later. This is a documented source of kinematic shape uncertainty. |
| C4 | **Decay File Mass Matching Uses Nearest Neighbor** | Minor | Mitigated | `rhn_decay_library.py` still selects pre-computed decay kinematics by finding the file with the nearest mass, without interpolation. This can cause kinematic discontinuities. The issue is mitigated by a denser grid and warnings, but the core logic is unchanged. |
| C5 | **Track Separation Uses Midpoint Approximation** | Minor | Verified | The function `compute_separation_pass_static` in `decay_detector.py` uses the midpoint of the HNL's path as a fixed decay position for efficiency. This is a documented, ctau-independent approximation. |
| C6 | **Kaon Cross-Section Uncertainty** | Moderate | Verified | The kaon production cross-section in `production_xsecs.py` is explicitly commented as being `"very approximate, soft QCD dominated"`. This remains a dominant systematic uncertainty for low-mass HNLs. |
| C10| **Form Factor Model Dependence** | Minor | Verified | `HNLCalc.py` uses hardcoded form factor parametrizations (e.g., single-pole dominance) for many semileptonic decays. This is a standard but model-dependent approach, introducing O(10-20%) systematic uncertainty. |

---

## Relevant Code Quality & Maintainability Issues

These are outstanding technical issues that do not cause incorrect physics results but violate best practices and could affect future maintenance or extension of the code.

| ID | Issue | Severity | Status | Verification Notes |
|----|-------|----------|--------|-------|
| B2 | **Unsafe `eval()` in HNLCalc** | Minor | Verified | The core integration functions in `HNLCalc.py` (e.g., `integrate_3body_br`) still use the built-in `eval()` on formula strings. While partially mitigated by a safe evaluator in the `hnl_model_hnlcalc.py` wrapper, the underlying `eval()` remains. |
| B5 | **Hardcoded Magic Numbers in Geometry** | Minor | Verified | The `build_drainage_gallery_mesh` function in `per_parent_efficiency.py` contains a large, hardcoded list of coordinates and transformation values defining the detector geometry, making it difficult to audit or modify. |
| B9 | **Global Mutable State for Mesh Caching** | Minor | Verified | A non-thread-safe global variable (`_MESH_CACHE`) is used in `expected_signal.py` to cache the detector mesh. This singleton pattern works for the current single-geometry use case but is not robust for parallel execution or multi-geometry scenarios. |
| B10| **Broad Exception Handling in Charge Lookup** | Minor | Verified | The `_charge_from_pdg` function in `decay_detector.py` uses a broad `except Exception:` block, which can silently mask underlying errors (e.g., from an unknown PDG ID) by defaulting to a neutral charge. |

---
*End of Verified Relevant Audit Summary*