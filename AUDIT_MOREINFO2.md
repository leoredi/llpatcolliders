## Part 4: Verification of Audit Report

**Date:** 2026-01-26
**Reviewer:** Gemini CLI

---

This section documents the verification of the issues raised in the initial audit report. Each point from the original report was checked against the codebase.

### Physics Biases Verification

-   **C1. W/Z Boson BR Formula:** **Verified.** The audit correctly identifies that the `W/Z` production branching ratio calculation sums the squared mixing elements. This is the physically correct approach for calculating the *inclusive* production rate of HNLs from `W/Z` decays. The implementation is consistent with the goal of computing a total production yield.

-   **C2. Phase-Space Kinematics in 3-Body Decays:** **Verified.** The audit correctly identifies that 3-body semileptonic decays are configured in Pythia using `meMode=0` (pure phase space). This is explicitly stated in the source code comments in `production/pythia_production/main_hnl_production.cc`. This is a known approximation.

-   **C3. τ → N X Representative Channels:** **Verified.** The `fromTau` production mode in `production/pythia_production/main_hnl_production.cc` uses a weighted mixture of representative decay channels for kinematic sampling, not as physical branching ratios, just as the audit states. The final normalization is applied later.

-   **C4. Decay File Mass Matching Uses Nearest Neighbor:** **Verified.** The function `_nearest_entry` in `analysis_pbc/decay/rhn_decay_library.py` selects a pre-computed decay kinematics file by finding the one with the closest mass. This is a "nearest-neighbor" approach, as correctly identified.

-   **C5. Track Separation Uses Midpoint Approximation:** **Verified.** The function `compute_separation_pass_static` in `analysis_pbc/decay/decay_detector.py` uses the midpoint of the HNL's path as a fixed reference decay position. This is an explicit, documented approximation for computational efficiency.

-   **C6. Kaon Cross-Section Uncertainty:** **Verified.** The file `analysis_pbc/config/production_xsecs.py` defines the kaon production cross-section with a comment explicitly stating it is a "very approximate, soft QCD" value. The audit is correct that this is a dominant systematic uncertainty.

-   **C7. Missing Neutral Kaon Channels:** **RESOLVED.** K_L (130) has been added to both Pythia production (`main_hnl_production.cc`: K_L → πℓN) and cross-sections (`production_xsecs.py`: SIGMA_KL_PB ≈ 0.5 × SIGMA_KAON_PB). K_S (310) intentionally omitted due to τ_S/τ_L ≈ 1/570 suppression.

-   **C8. Majorana vs Dirac Factor of 2:** **Verified.** The `dirac` flag, set via `--dirac`, is passed to `expected_signal_events` function in `analysis_pbc/limits/expected_signal.py`, where the total expected signal is multiplied by 2.0 if the flag is true. This is correctly implemented.

-   **C9. HNLCalc Integration Uses Monte Carlo with Fixed nsample=100:** **Verified.** The function `integrate_3body_br` in `analysis_pbc/HNLCalc/HNLCalc.py` uses a simple Monte Carlo method with a default of `nsample=100`, which can lead to statistical noise as noted.

-   **C10. Form Factor Model Dependence:** **Verified.** The `HNLCalc.py` script hardcodes form factor parametrizations for semileptonic decays, primarily using a single-pole dominance model, which introduces a model-dependent systematic uncertainty as stated.

### Coding Mistakes & Quality Issues Verification

-   **B1. Race Condition in Summary CSV Writing:** **Verified.** While `run_hnl_scan.py` is serial, if a user runs multiple instances in parallel for different flavours, a race condition exists when initializing and writing to the summary CSV file due to a lack of file locking.

-   **B2. Unsafe eval() in HNLCalc:** **Verified.** `HNLCalc.py` uses `eval()` to calculate branching ratios. The wrapper in `hnl_model_hnlcalc.py` uses a safe evaluator for the first level of string evaluation, but the final formula string is still processed by the unsafe `eval()` in `HNLCalc.py`. The audit's assessment is accurate.

-   **B3. Division by Zero Protection May Mask Issues:** **Verified.** The line `lam = np.where(lam <= 1e-9, 1e-9, lam)` in `analysis_pbc/limits/expected_signal.py` silently clamps the decay length to prevent division-by-zero errors, which can mask upstream data quality issues.

-   **B4. Float Precision in Mass String Formatting:** **Verified.** The function `format_mass_for_filename` in `config_mass_grid.py` uses f-string formatting which, due to floating-point inaccuracies, could lead to filename mismatches in edge cases.

-   **B5. Hardcoded Magic Numbers in Geometry:** **Verified.** The function `build_drainage_gallery_mesh` in `analysis_pbc/geometry/per_parent_efficiency.py` contains a large, hardcoded list of coordinates and transformation values, reducing readability and maintainability.

-   **B6. Missing Parent PDG Fallback in LHE Parser:** **Verified.** `lhe_to_csv.py` contains a fallback that assigns `W+` as the parent if the true parent is not found in the LHE record, which can lead to misattribution of off-shell Z events.

-   **B7. Inconsistent Column Naming (Legacy Support):** **Verified.** The codebase contains logic in multiple files (e.g., `combine_production_channels.py`, `per_parent_efficiency.py`) to handle both `boost_gamma` and `beta_gamma` column names for backward compatibility, adding complexity.

-   **B8. tqdm Progress Bar in Library Code:** **Verified.** The library module `analysis_pbc/geometry/per_parent_efficiency.py` directly uses `tqdm`, which is poor practice for code that may be used in parallel or non-interactive environments.

-   **B9. Global Mutable State in Mesh Caching:** **Verified.** The code in `analysis_pbc/limits/expected_signal.py` uses a non-thread-safe, mutable global variable (`_MESH_CACHE`) for caching the detector mesh, which is inefficient for multiprocessing.

-   **B10. Unhandled Exception in Particle Charge Lookup:** **Verified.** The function `_charge_from_pdg` in `analysis_pbc/decay/decay_detector.py` uses a broad `except Exception:` block, which can hide underlying issues by silently treating particles as neutral.

---
*End of Verification Report*
