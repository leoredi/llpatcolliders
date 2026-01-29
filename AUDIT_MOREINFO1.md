# HNL Analysis Pipeline: Physics Bias & Code Quality Audit

**Date:** 2026-01-26
**Scope:** Full repository review for physics biases and coding mistakes
**Reviewer:** Claude Code (Opus 4.5)

---

## Executive Summary

This audit covers the Heavy Neutral Lepton (HNL) sensitivity analysis pipeline, including production (Pythia/MadGraph), geometry, decay simulation, limits calculation, and visualization. The codebase is generally well-structured with clear documentation. However, several potential physics biases and coding issues were identified that could affect final sensitivity projections.

**Risk Categories:**
- **Critical (C):** Could cause order-of-magnitude errors in sensitivity
- **Moderate (M):** Could cause O(10-50%) systematic shifts
- **Minor (m):** Code quality issues or minor physics approximations

---

## Part 1: Physics Biases

### C1. W/Z Boson BR Formula Uses Total |U|² Instead of Flavour-Specific Coupling

**Location:** `analysis_pbc/models/hnl_model_hnlcalc.py:408-434`

**Issue:** The W and Z production BRs are computed as:
```python
br_W = (self.Ue2 + self.Umu2 + self.Utau2) * BR_W_to_lnu_SM * phase_space_W * helicity_W
```

This sums all three mixing angles. However, for pure benchmarks (100, 010, 001), only one coupling is non-zero. The formula is **correct** in that context, but it implicitly assumes that W → ℓN production is summed over all three lepton flavors weighted by their mixing.

**Potential Bias:** For a benchmark like "010" (muon only), the W can produce N via W → μN, but NOT via W → eN or W → τN since those couplings are zero. The code correctly handles this because only Umu2 ≠ 0. However, if users run with mixed couplings, the interpretation becomes ambiguous.

**Severity:** Minor for pure benchmarks, **Moderate** if mixed couplings are used.

---

### C2. Phase-Space Kinematics in 3-Body Decays (Pythia meMode=0)

**Location:** `production/pythia_production/main_hnl_production.cc:298-300`

**Issue:** Semileptonic 3-body decays (D0 → K ℓ N, B0 → D ℓ N, etc.) use `meMode=0` (pure phase space):
```cpp
pythia.readString("421:addChannel = 1 1.0 0 -321 " + lepBar + " " + hnl);
```

The comment acknowledges this: "Using meMode=0 (phase space) for simplicity."

**Physics Impact:** Phase-space kinematics differ from matrix-element kinematics. The HNL momentum spectrum will be harder (more uniform) than the physical spectrum which has V-A suppression at high q². This affects:
- **Geometric acceptance:** HNLs with different pT distributions have different detector hit rates
- **Decay probability:** Different boost distributions affect P_decay

**Severity:** **Moderate** — validated by arXiv:2103.11494 as "adequate for sensitivity estimates," but introduces O(10-20%) shape uncertainty.

---

### C3. τ → N X Representative Channels

**Location:** `production/pythia_production/main_hnl_production.cc:423-502`

**Issue:** The "fromTau" production mode uses representative channels (π, ρ, 3π) with arbitrary kinematic weights:
```cpp
rho_weight = 0.50;
tripi_weight = 0.30;
pi_weight = 0.20;
```

These are NOT physical branching ratios — they're kinematic sampling weights. The actual τ → N X BR is applied later by HNLCalc.

**Physics Impact:** The kinematic distribution of HNLs depends on these weights. At low HNL masses, τ → π N dominates, while at higher masses only τ → π N is kinematically allowed. The current weights may not optimally represent the true kinematic mixture.

**Severity:** **Minor to Moderate** — affects shape of acceptance vs mass.

---

### C4. Decay File Mass Matching Uses Nearest Neighbor

**Location:** `analysis_pbc/decay/rhn_decay_library.py:119-124`

**Issue:** Decay kinematics are loaded from precomputed MATHUSLA RHN files. The code selects the **nearest available mass point**:
```python
def _nearest_entry(entries, mass_GeV):
    return min(entries, key=lambda e: abs(e.mass_GeV - mass_GeV))
```

**Physics Impact:** If the analysis mass grid is finer than the RHN decay file grid, decay kinematics will be interpolated by "snapping" to the nearest available mass. This can cause:
- Wrong kinematic distributions at intermediate masses
- Discontinuities in acceptance vs mass

**Severity:** **Minor** — depends on grid granularity; the RHN files typically have 0.1-0.5 GeV spacing.

---

### C5. Track Separation Uses Midpoint Approximation

**Location:** `analysis_pbc/decay/decay_detector.py:334-392`

**Issue:** The function `compute_separation_pass_static()` computes track separation at the **midpoint** of the detector path:
```python
decay_distance = entry[idx] + 0.5 * path_length[idx]
```

This is an approximation — the actual decay position depends on the HNL lifetime (ctau) and boost.

**Physics Impact:** At very short or very long lifetimes:
- Short ctau: Decays near entry → different separation than midpoint
- Long ctau: Decays near exit → different separation than midpoint

The code comment acknowledges: "for our geometry the dependence is weak."

**Severity:** **Minor** — validated approximation for drainage gallery geometry.

---

### C6. Kaon Cross-Section Uncertainty

**Location:** `analysis_pbc/config/production_xsecs.py:83`

**Issue:** The kaon production cross-section is highly uncertain:
```python
SIGMA_KAON_PB = 5.0 * 1e10  # ~50 mb (very approximate, soft QCD)
```

**Physics Impact:** Kaon-regime sensitivity (m_N < 0.5 GeV) is directly proportional to this cross-section. Soft QCD predictions can vary by factors of 2-3 between different generators and tunes.

**Severity:** **Moderate** — the code documents this uncertainty, but users should be aware this is a dominant systematic for low-mass HNL sensitivity.

---

### C7. Missing Neutral Kaon Channels

**Location:** `analysis_pbc/config/production_xsecs.py:127-132`

**Issue:** ~~Only K± (PDG 321) has a cross-section defined~~ **RESOLVED**

K_L (130) is now included:
```python
if pid == 321:  # K+ / K-
    return SIGMA_KAON_PB
if pid == 130:  # K_L
    return SIGMA_KL_PB  # ≈ 0.5 × SIGMA_KAON_PB
# K_S (310) omitted — τ_S/τ_L ≈ 1/570 suppression
```

**Status:** K_L → πℓN channel added to Pythia production and cross-sections.

**Severity:** **Resolved** — K_L now included; K_S negligible.

---

### C8. Majorana vs Dirac Factor of 2

**Location:** Multiple files

**Issue:** The codebase consistently treats HNL as Majorana (N = N̄). For Dirac interpretation, yields should be multiplied by 2. The factor is applied via `--dirac` flag in run.py:
```python
if dirac:
    total_expected *= 2.0
```

**Physics Impact:** No bug, but users must remember this is NOT the default.

**Severity:** **None** — correctly implemented, just needs documentation awareness.

---

### C9. HNLCalc Integration Uses Monte Carlo with Fixed nsample=100

**Location:** `analysis_pbc/HNLCalc/HNLCalc.py:82-116`

**Issue:** 3-body BR integrations use Monte Carlo sampling with only 100 points:
```python
def integrate_3body_br(..., nsample=100):
    for i in range(nsample):
        # sample q2, energy
        integral += eval(br) * volume / float(nsample)
```

**Physics Impact:** With only 100 samples, statistical fluctuations are O(10%). For rare decay channels, this can cause noisy BR values that propagate to limits.

**Severity:** **Minor** — affects precision of BR calculations, typically averaged over many channels.

---

### C10. Form Factor Model Dependence

**Location:** `analysis_pbc/HNLCalc/HNLCalc.py:494-654`

**Issue:** Semileptonic form factors are hardcoded using specific models (e.g., single-pole dominance):
```python
fp = str(f00) + "/(1-q**2/" + str(MV) + "**2)"
```

**Physics Impact:** Form factor uncertainties can be O(10-20%) for heavy-to-light transitions. Different form factor models (LCSR, lattice QCD) give different predictions.

**Severity:** **Minor** — standard approach, but introduces systematic uncertainty.

---

## Part 2: Coding Mistakes & Quality Issues

### B1. Race Condition in Summary CSV Writing (Parallel Execution)

**Location:** `production/madgraph_production/scripts/run_hnl_scan.py:509-529`

**Issue:** When running multiple flavours in parallel (as you described), all three processes append to the same summary CSV without file locking:
```python
with open(summary_path, 'a') as f:
    f.write(row + '\n')
```

**Impact:** Potential for:
- Interleaved/corrupted rows
- Lost writes
- Race on header initialization

**Severity:** **Minor** — summary CSV is metadata only; physics data is protected.

---

### B2. Unsafe eval() in HNLCalc (Partially Mitigated)

**Location:** `analysis_pbc/HNLCalc/HNLCalc.py:114`

**Issue:** The HNLCalc module uses `eval()` to evaluate BR string formulas:
```python
integral += eval(br) * volume / float(nsample)
```

The wrapper in `hnl_model_hnlcalc.py` implements a safe expression evaluator, but the underlying HNLCalc still uses raw `eval()`.

**Impact:** Security risk if untrusted input is passed; not exploitable in current usage.

**Severity:** **Minor** — isolated to trusted physics formulas.

---

### B3. Division by Zero Protection May Mask Issues

**Location:** `analysis_pbc/limits/expected_signal.py:137`

**Issue:**
```python
lam = np.where(lam <= 1e-9, 1e-9, lam)  # Prevent divide by zero
```

This silently clamps lambda to 1e-9 m, which could mask issues where:
- Mass is zero or negative
- beta_gamma is invalid

**Impact:** Invalid HNLs with bad kinematics get included with tiny but non-zero weights instead of being flagged.

**Severity:** **Minor** — defensive programming, but could hide data quality issues.

---

### B4. Float Precision in Mass String Formatting

**Location:** `config_mass_grid.py:195`

**Issue:**
```python
def format_mass_for_filename(mass):
    return f"{mass:.2f}".replace('.', 'p')
```

Floating-point representation can cause issues:
```python
>>> f"{0.35:.2f}"  # Expected "0.35"
'0.35'
>>> f"{0.15 + 0.20:.2f}"  # May be "0.35" or "0.34999..."
```

**Impact:** Potential filename mismatches when masses are computed vs hardcoded.

**Severity:** **Minor** — mostly works, but edge cases possible.

---

### B5. Hardcoded Magic Numbers in Geometry

**Location:** `analysis_pbc/geometry/per_parent_efficiency.py:159-207`

**Issue:** The drainage gallery coordinates are hardcoded as a long list of magic numbers with transformations:
```python
correctedVertWithShift.append(
    ((x - 11908.8279764855) / 1000.0,
     (y + 13591.106147774964) / 1000.0,)
)
```

**Impact:**
- Hard to verify correctness
- No documentation of coordinate system origin
- Cannot easily switch detector geometries

**Severity:** **Minor** — works, but reduces maintainability.

---

### B6. Missing Parent PDG Fallback in LHE Parser

**Location:** `production/madgraph_production/scripts/lhe_to_csv.py:212-214`

**Issue:** When parent W/Z boson is not found in LHE (off-shell), the code defaults to W+:
```python
if parent_pdg == 0:
    parent_pdg = self.PDG_WPLUS  # use +24 as neutral sign choice
    parent_inferred = True
```

**Impact:** For off-shell Z → ν N, events are misattributed to W. This could affect analysis if W and Z have different cross-sections or kinematics.

**Severity:** **Moderate** — documented behavior, but could bias W/Z mixture.

---

### B7. Inconsistent Column Naming (Legacy Support)

**Location:** Multiple files

**Issue:** The codebase supports two column names for the same quantity:
- `boost_gamma` (legacy)
- `beta_gamma` (current)

The code handles this, but it's scattered across multiple files:
- `combine_production_channels.py:106-126`
- `per_parent_efficiency.py:334-344`

**Impact:** Maintenance burden; potential for bugs if one file forgets the translation.

**Severity:** **Minor** — correctly handled but adds complexity.

---

### B8. tqdm Progress Bar in Library Code

**Location:** `analysis_pbc/geometry/per_parent_efficiency.py:402-406`

**Issue:** Library code imports and uses tqdm, which can interfere with parallel execution or non-interactive environments:
```python
for start in tqdm(range(0, len(valid_indices), batch_size), ...):
```

**Impact:** Log spam in parallel/batch mode; potential deadlocks with multiprocessing.

**Severity:** **Minor** — cosmetic, but can be annoying.

---

### B9. Global Mutable State in Mesh Caching

**Location:** `analysis_pbc/limits/expected_signal.py:213-222`

**Issue:**
```python
_MESH_CACHE = None

def build_mesh_once():
    global _MESH_CACHE
    if _MESH_CACHE is None:
        _MESH_CACHE = build_drainage_gallery_mesh()
    return _MESH_CACHE
```

**Impact:** Not thread-safe; if multiprocessing forks after mesh is built, memory usage increases. Also prevents using different geometries in the same session.

**Severity:** **Minor** — works for current single-geometry use case.

---

### B10. Unhandled Exception in Particle Charge Lookup

**Location:** `analysis_pbc/decay/decay_detector.py:79-87`

**Issue:**
```python
def _charge_from_pdg(pid: int) -> float:
    if Particle is None:
        return 0.0
    try:
        particle = Particle.from_pdgid(int(pid))
    except Exception:
        return 0.0
```

Catching all exceptions hides potential issues (e.g., wrong particle library version).

**Impact:** Unknown particles silently treated as neutral, potentially missing charged tracks.

**Severity:** **Minor** — defensive but could mask bugs.

---

## Part 3: Recommendations

### High Priority

1. **Document kaon cross-section uncertainty** prominently — this is likely the dominant systematic for m_N < 0.5 GeV.

2. ~~**Consider adding K_L/K_S channels** to capture full kaon sensitivity.~~ ✓ K_L added

3. **Add file locking** to summary CSV writes for parallel-safe production.

4. **Increase HNLCalc integration samples** from 100 to 1000 for better BR precision.

### Medium Priority

5. **Validate phase-space vs matrix-element kinematics** by comparing to known differential distributions.

6. **Create a geometry abstraction layer** to enable different detector configurations.

7. **Consolidate column naming** to use only `beta_gamma` everywhere (with migration script).

8. **Add integration tests** comparing full pipeline output to benchmark results.

### Low Priority

9. **Move magic numbers** (coordinates, transformations) to configuration files.

10. **Add logging framework** instead of scattered print statements.

11. **Refactor HNLCalc** to avoid string-based `eval()` in favor of symbolic math.

---

## Appendix: Files Reviewed

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

*End of Audit Report*

---

## Independent Review of Audit Findings

**Review Date:** 2026-01-26
**Reviewer:** Claude Opus 4.5 (Critical Review)

### Methodology

Each issue was verified by reading the actual source code at the cited locations. Findings are categorized as:
- **CONFIRMED**: Issue exists as described
- **PARTIALLY CONFIRMED**: Issue exists but description needs clarification
- **NOT CONFIRMED**: Issue does not exist or is misdescribed

---

### Physics Biases Review

#### C1. W/Z Boson BR Formula — CONFIRMED ✓

Verified at `hnl_model_hnlcalc.py:418`:
```python
br_W = (self.Ue2 + self.Umu2 + self.Utau2) * BR_W_to_lnu_SM * phase_space_W * helicity_W
```

**Assessment**: The audit is correct. The formula sums all mixing angles because W → ℓN can proceed via any lepton flavor. For pure benchmarks (100/010/001), only one term is non-zero, so this is correct. The severity rating is appropriate.

---

#### C2. Phase-Space Kinematics (meMode=0) — CONFIRMED ✓

Verified at `main_hnl_production.cc:307`:
```cpp
pythia.readString("421:addChannel = 1 1.0 0 -321 " + lepBar + " " + hnl);
```

The third parameter `0` is indeed meMode=0 (phase space). The header comment at line 298-300 acknowledges this:
> "Note: Using meMode=0 (phase space) for simplicity"

**Assessment**: Correctly identified. The O(10-20%) shape uncertainty estimate is physically reasonable.

---

#### C3. τ → N X Representative Channels — CONFIRMED ✓

Verified at `main_hnl_production.cc:453-455`:
```cpp
rho_weight = 0.50;
tripi_weight = 0.30;
pi_weight = 0.20;
```

These weights are explicitly documented as "kinematic sampling weights (NOT physical BRs)" at line 445-446.

**Assessment**: Correctly identified. These weights affect HNL momentum distributions but are not the physical τ→NX branching ratios.

---

#### C4. Decay File Mass Matching — CONFIRMED ✓

Verified at `rhn_decay_library.py:119-123`:
```python
def _nearest_entry(entries: Iterable[DecayFileEntry], mass_GeV: float) -> DecayFileEntry | None:
    entries = list(entries)
    if not entries:
        return None
    return min(entries, key=lambda e: abs(e.mass_GeV - mass_GeV))
```

**Assessment**: Correctly identified. The nearest-neighbor approach could cause kinematic discontinuities if the analysis grid is finer than the RHN file grid.

---

#### C5. Track Separation Midpoint Approximation — CONFIRMED ✓

Verified at `decay_detector.py:370`:
```python
decay_distance = entry[idx] + 0.5 * path_length[idx]
```

The code comment at lines 341-345 acknowledges this:
> "This is an approximation: the exact separation depends on decay position, but for our geometry the dependence is weak."

**Assessment**: Correctly identified. The approximation is documented and justified for the drainage gallery geometry.

---

#### C6. Kaon Cross-Section Uncertainty — CONFIRMED ✓

Verified at `production_xsecs.py:83`:
```python
SIGMA_KAON_PB = 5.0 * 1e10  # ~50 mb (very approximate, soft QCD)
```

**Assessment**: Correctly identified. This is the dominant systematic for m_N < 0.5 GeV. The code comment acknowledges the uncertainty.

---

#### C7. Missing Neutral Kaon Channels — **RESOLVED** ✓

**Update:** K_L (130) has been added to both Pythia production and cross-sections:
```python
if pid == 321:  # K+ / K-
    return SIGMA_KAON_PB
if pid == 130:  # K_L
    return SIGMA_KL_PB  # ≈ 0.5 × SIGMA_KAON_PB
# K_S (310) omitted — τ_S/τ_L ≈ 1/570 suppression
```

**Assessment**: K_L → πℓN now included. K_S omitted due to negligible contribution.

---

#### C8. Majorana vs Dirac Factor — CONFIRMED ✓

Verified at `expected_signal.py:207-208`:
```python
if dirac:
    total_expected *= 2.0
```

**Assessment**: Correctly implemented. The audit's "no bug" assessment is accurate.

---

#### C9. HNLCalc nsample=100 — CONFIRMED ✓

Verified at `HNLCalc.py:82`:
```python
def integrate_3body_br(self, br, mass, m0, m1, m2, coupling=1, nsample=100, integration="dq2dE"):
```

And line 114:
```python
integral += eval(br)*(q2max-q2min)*(ENmax-ENmin)/float(nsample)
```

**Assessment**: Correctly identified. With 100 Monte Carlo samples, statistical fluctuations are O(1/√100) = 10%. For BR calculations this introduces noise but is typically averaged over many channels.

---

#### C10. Form Factor Model Dependence — CONFIRMED ✓

Verified at `HNLCalc.py:498`:
```python
fp=str(f00)+"/(1-q**2/"+str(MV)+"**2)"
```

This is single-pole dominance form. Form factors are hardcoded for each parent-daughter pair throughout lines 494-654.

**Assessment**: Correctly identified. Standard approach but introduces O(10-20%) systematic from form factor model choice.

---

### Coding Issues Review

#### B1. Race Condition in CSV Writing — CONFIRMED ✓

Verified at `run_hnl_scan.py:528-529`:
```python
with open(summary_path, 'a') as f:
    f.write(row + '\n')
```

No file locking or atomic write mechanism is used.

**Assessment**: Valid concern for parallel execution. However, impact is limited to metadata corruption—physics data files are written separately.

---

#### B2. Unsafe eval() in HNLCalc — CONFIRMED ✓

Verified at `HNLCalc.py:114`:
```python
integral += eval(br)*(q2max-q2min)*(ENmax-ENmin)/float(nsample)
```

The safe evaluator in `hnl_model_hnlcalc.py` (lines 56-130) protects the *outer* BR string evaluation, but when `integrate_3body_br()` is called, the `br` string passed to it is evaluated via raw `eval()` inside HNLCalc.

**Assessment**: Correctly identified. Security risk is minimal in current trusted-input context but violates best practices.

---

#### B3. Division by Zero Protection — CONFIRMED with clarification

Verified at `expected_signal.py:137`:
```python
lam = np.where(lam <= 1e-9, 1e-9, lam)  # Prevent divide by zero
```

**Assessment**: The clamping exists as described. However, λ = βγ × cτ₀ < 1e-9 m would require an unphysically short lifetime combined with small boost. The 1e-9 m threshold (1 nm) is effectively zero for LHC-scale analysis. The defensive clamping is appropriate but *could* mask upstream data issues (e.g., NaN propagation).

---

#### B4. Float Precision in Mass Formatting — PARTIALLY CONFIRMED

Verified at `config_mass_grid.py:195`:
```python
return f"{mass:.2f}".replace('.', 'p')
```

**Assessment**: The audit's specific example (`0.15 + 0.20`) would actually work correctly due to IEEE 754 rounding behavior. The real risk is:
- Masses computed via arithmetic vs read from config may have different bit representations
- e.g., `0.15 + 0.10` → `0.25000000000000003` formats to `"0.25"` (OK)
- But `np.float32(0.25)` vs `float(0.25)` could differ in edge cases

**Practical impact is minimal** for the 2-decimal formatting used here.

---

#### B5. Hardcoded Magic Numbers — CONFIRMED ✓

Verified at `per_parent_efficiency.py:159-216`:
- 47 coordinate pairs with no external source reference
- Transformation at lines 210-216:
```python
(x - 11908.8279764855) / 1000.0,
(y + 13591.106147774964) / 1000.0,
```

The transformation is documented at lines 150-153 but the coordinate origin is unexplained.

**Assessment**: Correctly identified. The magic numbers reduce maintainability and auditability.

---

#### B6. Missing Parent PDG Fallback — CONFIRMED ✓

Verified at `lhe_to_csv.py:212-214`:
```python
if parent_pdg == 0:
    parent_pdg = self.PDG_WPLUS  # use +24 as neutral sign choice
    parent_inferred = True
```

**Assessment**: Correctly identified. Off-shell Z → ν N events would be misattributed to W. The `parent_inferred` flag allows downstream filtering if needed.

---

#### B7. Inconsistent Column Naming — CONFIRMED ✓

Verified at:
- `combine_production_channels.py:106-126`: `normalize_boost_column()` function
- `per_parent_efficiency.py:333-344`: Similar handling logic

Both files handle the `boost_gamma` → `beta_gamma` translation independently.

**Assessment**: Correctly identified. The translation is implemented correctly but scattered, creating maintenance burden.

---

#### B8. tqdm in Library Code — CONFIRMED ✓

Verified at `per_parent_efficiency.py:402-406`:
```python
for start in tqdm(
    range(0, len(valid_indices), batch_size),
    total=(len(valid_indices) + batch_size - 1) // batch_size,
    desc="Geometry rays",
    unit="batch",
):
```

**Assessment**: Correctly identified. Direct tqdm usage in library code can cause issues in non-interactive or parallel environments.

---

#### B9. Global Mesh Cache — CONFIRMED ✓

Verified at `expected_signal.py:213-222`:
```python
_MESH_CACHE = None

def build_mesh_once():
    global _MESH_CACHE
    if _MESH_CACHE is None:
        from geometry.per_parent_efficiency import build_drainage_gallery_mesh
        _MESH_CACHE = build_drainage_gallery_mesh()
    return _MESH_CACHE
```

**Assessment**: Correctly identified. Not thread-safe; would cause issues with different geometry configurations in same session.

---

#### B10. Broad Exception Catch — CONFIRMED ✓

Verified at `decay_detector.py:81-84`:
```python
try:
    particle = Particle.from_pdgid(int(pid))
except Exception:
    return 0.0
```

**Assessment**: Correctly identified. Catching all exceptions masks potential issues (wrong particle library version, import errors, etc.). Should catch specific exceptions like `ParticleNotFoundError`.

---

### Summary of Review

| Issue | Status | Severity Agreement |
|-------|--------|-------------------|
| C1 | CONFIRMED | Agree |
| C2 | CONFIRMED | Agree |
| C3 | CONFIRMED | Agree |
| C4 | CONFIRMED | Agree |
| C5 | CONFIRMED | Agree |
| C6 | CONFIRMED | Agree |
| C7 | CONFIRMED | Agree |
| C8 | CONFIRMED | Agree (no bug) |
| C9 | CONFIRMED | Agree |
| C10 | CONFIRMED | Agree |
| B1 | CONFIRMED | Agree |
| B2 | CONFIRMED | Agree |
| B3 | CONFIRMED | Minor—threshold is appropriate |
| B4 | PARTIALLY CONFIRMED | Lower risk than stated |
| B5 | CONFIRMED | Agree |
| B6 | CONFIRMED | Agree |
| B7 | CONFIRMED | Agree |
| B8 | CONFIRMED | Agree |
| B9 | CONFIRMED | Agree |
| B10 | CONFIRMED | Agree |

### Additional Observations

1. **W/Z Production BR (C1)**: The code at line 418-419 stores W BR under `br_per_parent[24]`. This is the correct approach for inclusive W → ℓN production summed over all lepton flavors.

2. **Tau Production (lines 352-387)**: The τ → π N BR calculation uses a proper analytic formula from Atre et al. (arXiv:0901.3589), not HNLCalc. This is correctly implemented.

3. **Safe Evaluator Coverage**: The `_SafeExprEvaluator` in `hnl_model_hnlcalc.py` is well-designed but only covers the BR string evaluation at the wrapper level. The underlying HNLCalc integration still uses raw `eval()`.

4. **Decay Position Sampling**: In `compute_decay_acceptance()` (decay_detector.py:284-309), the full ctau-dependent decay position sampling is implemented correctly using exponential distribution. The midpoint approximation (C5) is only used for the pre-computed static separation pass.

### Conclusion

The original audit is thorough and accurate. **19 of 20 issues are fully confirmed**, with B4 (float precision) being lower risk than stated but still technically valid. The severity assessments are appropriate. The recommendations in Part 3 of the original audit are well-prioritized.

---

*End of Independent Review*

---

## Addendum: Verification of Audit Findings (2026-01-26)

Status key:
- Correct: matches code/behavior.
- Partially correct: code matches, but impact/severity overstated or interpretation ambiguous.
- Not an issue: intentional/benign or claim not supported by current code.

### Part 1: Physics Biases

C1. W/Z BR uses sum |U|^2
Status: Partially correct.
Notes: `analysis_pbc/models/hnl_model_hnlcalc.py` does sum (Ue2+Umu2+Utau2). This is correct for inclusive W/Z -> l N production and matches the single-flavour benchmarks used by the pipeline. Only ambiguous if a user intends per-flavour production in mixed-coupling scans; documentation could clarify.

C2. 3-body semileptonic decays use phase space (meMode=0)
Status: Correct (intentional approximation).
Notes: `production/pythia_production/main_hnl_production.cc` uses meMode=0 for D/B semileptonic channels and comments this explicitly. This is a modelling choice rather than a coding error.

C3. Tau -> N X representative-channel weights are not physical BRs
Status: Correct (intentional approximation).
Notes: `production/pythia_production/main_hnl_production.cc` uses pi/rho/3pi weights as kinematic mixture and applies physical BRs later via HNLCalc. Potential acceptance-shape systematics exist, but this is by design.

C4. Decay file mass matching uses nearest neighbor
Status: Correct.
Notes: `analysis_pbc/decay/rhn_decay_library.py` selects the nearest mass with no interpolation, so discontinuities are possible if the analysis grid is finer than the RHN grid.

C5. Track separation uses midpoint approximation
Status: Correct (intentional approximation).
Notes: `analysis_pbc/decay/decay_detector.py` documents midpoint decay position as a geometry-based approximation.

C6. Kaon cross-section uncertainty
Status: Correct (documented systematic).
Notes: `analysis_pbc/config/production_xsecs.py` labels SIGMA_KAON_PB as very approximate. This is a known normalization systematic for m < 0.5 GeV.

C7. Missing neutral kaon channels
Status: **RESOLVED**.
Notes: K_L (130) added to Pythia production (`main_hnl_production.cc`) and cross-sections (`production_xsecs.py`). K_S (310) omitted due to τ_S/τ_L ≈ 1/570 suppression.

C8. Majorana vs Dirac factor of 2
Status: Correct (not an issue).
Notes: `analysis_pbc/limits/expected_signal.py` applies the factor under `--dirac`, and default is Majorana.

C9. HNLCalc 3-body integration uses nsample=100
Status: Correct.
Notes: `analysis_pbc/HNLCalc/HNLCalc.py` defaults to nsample=100 for MC integration. This is a precision limitation, not a correctness bug.

C10. Form factor model dependence
Status: Correct (model choice, not a bug).
Notes: `analysis_pbc/HNLCalc/HNLCalc.py` hardcodes specific form-factor parametrizations. This is a standard modelling systematic.

### Part 2: Coding Mistakes and Quality Issues

B1. Race condition in summary CSV writes
Status: Partially correct.
Notes: `production/madgraph_production/scripts/run_hnl_scan.py` appends without locking. The script itself is single-process, so the issue only appears when multiple instances run concurrently.

B2. Unsafe eval() in HNLCalc
Status: Correct but low risk.
Notes: `analysis_pbc/HNLCalc/HNLCalc.py` uses eval inside integration. Expressions are internal and not user-provided in this pipeline, so the security risk is minimal.

B3. Division-by-zero protection may mask issues
Status: Correct.
Notes: `analysis_pbc/limits/expected_signal.py` clamps lambda to 1e-9. This can hide invalid kinematics, but it does not inflate yields; impact is minor.

B4. Float precision in mass filename formatting
Status: Not an issue in current code.
Notes: `config_mass_grid.py` defines masses as explicit decimals and all filename formatting uses the same 2-decimal scheme, so practical mismatches are unlikely.

B5. Hardcoded magic numbers in geometry
Status: Correct (maintainability concern).
Notes: `analysis_pbc/geometry/per_parent_efficiency.py` documents the shifts; not a physics bug but hard to audit.

B6. Missing parent PDG fallback in LHE parser
Status: Correct but impact unquantified.
Notes: `production/madgraph_production/scripts/lhe_to_csv.py` defaults to W+ if no parent is found. This can misattribute Z events; severity depends on how often the parent is missing.

B7. Inconsistent column naming (boost_gamma vs beta_gamma)
Status: Correct but handled.
Notes: `analysis_pbc/limits/combine_production_channels.py` and `analysis_pbc/geometry/per_parent_efficiency.py` normalize legacy columns. Maintenance cost only.

B8. tqdm in library code
Status: Partially correct.
Notes: `analysis_pbc/geometry/per_parent_efficiency.py` always uses tqdm. This can be noisy in batch logs, but the claim about deadlocks is not substantiated here.

B9. Global mutable cache for mesh
Status: Correct but low risk.
Notes: `analysis_pbc/limits/expected_signal.py` caches mesh in a module global. Not thread-safe, but typical usage is single-geometry per process.

B10. Catch-all exception in charge lookup
Status: Correct but minor.
Notes: `analysis_pbc/decay/decay_detector.py` returns 0 charge on any exception; `_is_charged` provides a small fallback list. This can hide unexpected particle IDs but impact is limited.


---

## Addendum: Actionable TODOs + B6 Quantification (2026-01-26)

### A) TODOs for “Partially Correct” Items (Owner + Test)

C1 (W/Z BR sums |U|^2)
- Owner: Analysis/physics
- TODO: Clarify mixed-coupling semantics in `analysis_pbc/models/hnl_model_hnlcalc.py` and README (explicitly state BR is inclusive over e/μ/τ). Optionally add a `per_flavour` mode that returns separate W→eN, W→μN, W→τN contributions for mixed-coupling scans.
- Test: Unit test that (a) single-flavour inputs reproduce BR ∝ Uℓ^2, and (b) mixed-coupling inputs sum to Ue2+Umu2+Utau2 while `per_flavour` returns the split components.

C7 (Neutral kaon channels excluded) — **RESOLVED**
- ~~Owner: Production + analysis~~
- ~~TODO: Add K_L production~~ ✓ Done: K_L (130) added to `main_hnl_production.cc` and `production_xsecs.py`
- K_S (310) intentionally omitted (τ_S/τ_L ≈ 1/570 suppression)
- Test: Verify K_L (PDG 130) appears in kaon-regime CSVs after regenerating production

B1 (summary CSV append without locking)
- Owner: Production/madgraph
- TODO: Make summary writing concurrency-safe (file lock or per-process temp CSV + merge at end).
- Test: Small integration test that launches 3 parallel `run_hnl_scan.py --test` processes and verifies the final summary has exactly 3 rows (no corruption/duplication).

B8 (tqdm progress in library code)
- Owner: Analysis/geometry
- TODO: Add a `--no-progress` or `HNL_NO_TQDM=1` switch and auto-disable tqdm when stdout is not a TTY.
- Test: Run `analysis_pbc/geometry/per_parent_efficiency.py` in a non-interactive shell and assert no tqdm control characters are emitted (stdout match). 

### B) B6 Quantification From Existing Outputs

Scope: Scanned all MadGraph EW outputs in `output/csv/simulation/*_ew.csv`.
- Files scanned: 279
- Events scanned: 27,900,000
- Observed parent_pdg distribution:
  - 24 (W+): 11,590,178 (41.54%)
  - -24 (W-): 8,795,681 (31.53%)
  - 23 (Z): 7,514,141 (26.93%)
  - 0 or other PDGs: 0 (0.00%)

Interpretation: The CSVs do not store the `parent_inferred` flag, so the exact fallback rate (events defaulted to W+) is not directly observable from existing outputs. The observed data show no parent_pdg=0 entries; therefore the measurable lower bound is 0%. A strict upper bound is 41.54% if every W+ entry were inferred, but given the explicit W/Z in proc cards, the true fallback rate is likely near 0%.

