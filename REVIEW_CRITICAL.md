# HNL LLP Pipeline Critical Review

**Date:** 2025-12-10
**Reviewer:** Claude Code (Opus 4.5)
**Scope:** Bugs, failure modes, and bad assumptions in production and analysis code

---

## Summary Matrix

|                          | **A. Production (Pythia + MadGraph)**                              | **B. Analysis (geometry / limits / plots)**                         |
|--------------------------|---------------------------------------------------------------------|----------------------------------------------------------------------|
| **1. Code / SW Engineering** | parent_pdg=0 data loss (logged); race condition in geometry cache; fragile regex for xsec extraction → [§1A](#1a-code--production) | Slow row-by-row iteration in geometry; eval() security pattern in HNLCalc wrapper; eps2 grid resolution → [§1B](#1b-code--analysis) |
| **2. Physics / Phenomenology** | τ→πN only approximation; fromTau BR=0 (missing tau in production_brs); Λc channels not configured; phase-space 3-body kinematics → [§2A](#2a-physics--production) | Majorana factor=2 never applied; no detector/reco efficiency or backgrounds; W/Z BR approximations → [§2B](#2b-physics--analysis) |

---

## Updates (2025-12-XX, session with Codex)

- Fixed tau BR zeroing: `analysis_pbc/models/hnl_model_hnlcalc.py` now includes τ→πN BR for PDG 15, so fromTau samples contribute.
- Fixed off-shell W/Z losses: `production/madgraph_production/scripts/lhe_to_csv.py` defaults missing parent_pdg to 24 and logs inferred parents to prevent BR=0 drops.
- Added Λc semileptonic production: `production/pythia_production/main_hnl_production.cc` now forces Λc→ΛℓN when kinematically allowed.
- Improved tau fromTau kinematics: added τ→ρN and τ→3πN representative modes (weights split) in `production/pythia_production/main_hnl_production.cc`.
- Enabled form-factor semileptonic ingestion: combine/analysis accept and prefer `*_ff` CSVs; documented the optional third production step in `production/madgraph_production/FORMFACTOR_SEMILEPTONIC.md` and noted it in `run_hnl_scan.py` docstring. No MG generator/cards provided yet—user must supply MG FF samples.

## Updates (2025-12-11, session with Claude Opus 4.5)

- Fixed geometry cache race condition: `analysis_pbc/limits/run_serial.py` now uses atomic writes (tempfile + os.replace) instead of check-then-write pattern.
- Added Majorana/Dirac flag: `run_serial.py --dirac` applies ×2 yield factor for Dirac HNL interpretation; defaults to Majorana.
- Removed K_S/K_L from cross-section lookup: `analysis_pbc/config/production_xsecs.py` no longer returns σ for PDG 310/130 since these aren't configured in Pythia (fixes §2A.5 inconsistency).
- Made Pythia card files mandatory: `production/pythia_production/main_hnl_production.cc` now exits with error if card files not found, instead of silently using different physics settings (fixes §1A.4).

---

## §1A: Code / Production

### 1A.1 — parent_pdg=0 causes signal loss (MadGraph LHE)

**File:** `production/madgraph_production/scripts/lhe_to_csv.py:196-207`

**Issue:** When W/Z bosons are off-shell (controlled by `bw_cut` in MadGraph), they may not appear in the LHE record. The parser sets `parent_pdg=0`. Downstream, `production_brs()` has no entry for PDG=0, so these events contribute BR=0.

**Logging:** Warnings ARE printed at `lhe_to_csv.py:312-314` and `u2_limit_calculator.py:222-225`, but the latter only prints at the first eps2 scan point (1e-12), so can be missed in logs.

**Risk:** For high HNL masses near m_W threshold, some events could have off-shell W, causing underestimated sensitivity. Severity depends on `bw_cut` setting in MadGraph.

**Fix:** Either:
1. Force parent_pdg=24 for all EW events (known to be W/Z production)
2. Add fallback in `production_brs()`: `br_per_parent.setdefault(0, br_per_parent.get(24, 0.0))`

---

### 1A.2 — Geometry cache race condition in parallel mode

**File:** `analysis_pbc/limits/run_serial.py:80-82`

**Issue:**
```python
if not geom_csv.exists():
    geom_df.to_csv(geom_csv, index=False)  # Race: multiple workers check & write
```
Two parallel workers can simultaneously check `exists()=False`, then both write. Last write wins, potentially with incomplete data if IO interrupted.

**Fix:** Use atomic write pattern:
```python
import tempfile
with tempfile.NamedTemporaryFile(delete=False, dir=geom_csv.parent) as tmp:
    geom_df.to_csv(tmp.name, index=False)
    os.rename(tmp.name, geom_csv)  # Atomic on POSIX
```

---

### 1A.3 — Cross-section regex extraction fragile

**File:** `production/madgraph_production/scripts/run_hnl_scan.py:408-439`

**Issue:** Multiple regex patterns tried in sequence. If MadGraph output format changes (new MG5 version), extraction silently returns `None`, logged as warning but run continues.

**Fix:** Add assertion or structured parser for the banner XML format.

---

### 1A.4 — Pythia fallback uses different physics silently

**File:** `production/pythia_production/main_hnl_production.cc:530-557`

**Issue:** If card files not found, code falls back to hardcoded settings with only a "Warning" but different Tune:pp, different process selection. Users may not notice physics differs.

**Fix:** Make missing card file a hard error, or print clear "USING FALLBACK SETTINGS" banner.

---

### 1A.5 — findPhysicalParent safety checks

**File:** `production/pythia_production/main_hnl_production.cc:180-197`

**Issue:** Returns 0 if mother chain is broken or particle index out of bounds (lines 181, 189). These events get `parent_pdg=0` in CSV → downstream BR=0 → lost signal.

**Severity:** LOW. In normal Pythia operation, HNLs from meson decays should always have valid mothers. This is primarily a safety check against corrupted events.

**Fix (optional):** Add fallback to infer parent from production regime (e.g., beauty regime → assume B meson PDG). Or simply log when this happens to verify it's rare.

---

## §1B: Code / Analysis

### 1B.1 — Slow row-by-row geometry iteration

**File:** `analysis_pbc/geometry/per_parent_efficiency.py:323-378`

**Issue:** Using `for idx, row in df.iterrows()` is O(N) with DataFrame overhead per row. For 100k+ HNLs, this is very slow.

**Impact:** Geometry preprocessing becomes the bottleneck; users wait hours unnecessarily.

**Fix:** Vectorize using trimesh batch ray-casting:
```python
directions = np.vstack([eta_phi_to_direction(e, p) for e, p in zip(df['eta'], df['phi'])])
locations, _, _ = mesh.ray.intersects_location(ray_origins=origins, ray_directions=directions)
```

---

### 1B.2 — eval() usage in HNLCalc wrapper

**File:** `analysis_pbc/models/hnl_model_hnlcalc.py:194-203`

**Issue:** Uses Python `eval()` on strings from HNLCalc. While restricted namespace is used, this is a security anti-pattern. If HNLCalc were compromised (e.g., malicious update), arbitrary code execution is possible.

**Risk:** Low (HNLCalc is vendored/trusted), but code review would flag this.

**Fix:** Parse the BR formula structure explicitly rather than using eval().

---

### 1B.3 — Weight sanity check (working correctly)

**File:** `analysis_pbc/limits/u2_limit_calculator.py:343-352`

**Behavior:**
- `w_max > 1e6`: prints ERROR, returns None (correctly stops processing)
- `w_max > 1000`: prints WARN, continues (appropriate for unusual but possibly valid weights)

**Status:** This is working as intended. The two-tier approach (ERROR/stop vs WARN/continue) is reasonable.

**Note:** Original review incorrectly stated this "logs ERROR but continues" - that was wrong. The ERROR case correctly returns None.

---

### 1B.4 — eps2 grid resolution may miss sensitivity edges

**File:** `analysis_pbc/limits/u2_limit_calculator.py:248`

**Issue:** Fixed 100-point log grid from 10^-12 to 10^-2. For very weak couplings, the first excluded point may not be at the true boundary.

**Fix:** Add adaptive refinement near exclusion boundary.

---

### 1B.5 — Plot styling minor issues

**File:** `money_plot/plot_money_island.py`

**Issues (MINOR):**
1. Line 55-58: Tau xlim [0.5, 50] is reasonable - data simply won't exist below the kinematic threshold
2. No error bars or systematic uncertainty bands shown
3. Legend labels could be clearer ("Too long-lived" / "Too prompt" vs "Lower bound" / "Upper bound")

**Severity:** LOW. These are cosmetic issues, not physics bugs. The xlim is fine - matplotlib will show empty space below the first data point.

**Fix:** Add uncertainty bands if systematic uncertainties are quantified. Otherwise, leave as-is.

---

## §2A: Physics / Production

### 2A.1 — fromTau events have zero signal weight (missing tau BR)

**Files:**
- `analysis_pbc/models/hnl_model_hnlcalc.py:257-304` (production_brs() missing tau)
- `analysis_pbc/config/production_xsecs.py:147` (σ_tau IS defined)

**Issue:** For tau coupling at m < 1.64 GeV, the dominant production is through the cascade:
```
B/Ds → τ ν,  then  τ → π N
```
This "fromTau" mode is correctly simulated by Pythia and correctly combined by `combine_production_channels.py`. However, fromTau events have `parent_pdg = 15` (tau), and:

- `production_xsecs.py` **HAS** tau cross-section: `σ_tau = σ_Ds × BR(Ds→τν)` ✓
- `production_brs()` **MISSING** tau BR: no entry for PDG 15, returns 0.0 ✗

Result: `N_sig = L × σ_tau × BR(τ→πN) × ε = ... × 0.0 = 0`

Even though fromTau events are in the combined file, they contribute **zero signal** because `production_brs()` doesn't have tau decay channels.

**Note:** The `"_fromTau" not in f.name` filter in run_serial.py is a safety net to prevent double-counting if combine wasn't run. It's working as designed.

**Fix:** Add tau (PDG 15) to `production_brs()` with BR(τ → N X) formula, similar to how W/Z are added:
```python
# τ → π N (and other hadronic modes)
m_tau = 1.777
m_pi = 0.140
if mass < m_tau - m_pi:
    # BR(τ → π N) ∝ |U_τ|² × phase_space × form_factors
    # Reference: arXiv:1805.08567 Eq. 2.19
    br_tau = self.Utau2 * compute_tau_to_piN_br(mass, m_tau, m_pi)
    br_per_parent[15] = br_tau
```

---

### 2A.2 — τ → πN only for fromTau mode

**File:** `production/pythia_production/main_hnl_production.cc:424-436`

**Issue:** Uses only τ → π N channel. Real tau HNL decays include τ → ρN (dominant), τ → 3πN, τ → ℓννN, etc. Using only π channel biases the HNL momentum distribution.

**Impact:** Geometric acceptance ε_geom will be systematically biased (likely pessimistic since π is softer than ρ).

**Fix:** Add τ → ρ N channel (PDG ρ⁺ = 213) with appropriate BR weighting, or use external decay tools.

---

### 2A.3 — Λc channels not configured

**File:** `production/pythia_production/main_hnl_production.cc:104-108, 217-390`

**Issue:** BARYONS_3BODY defines Λc (PDG 4122), but `configureMesonDecays()` doesn't add Λc → Λ ℓ N or Λc → p K ℓ N channels. Only Λb → Λc is configured.

**Impact:** Missing charm baryon production channel. Λc cross-section is ~6% of ccbar, so this is a few-percent effect.

**Fix:** Add Λc decay channels in configureMesonDecays().

---

### 2A.4 — Phase-space 3-body kinematics

**File:** `production/pythia_production/main_hnl_production.cc:297-299`

**Issue:** Uses meMode=0 (pure phase space) for semileptonic decays like B → D ℓ N. Real matrix elements have form-factor dependence that affects the lepton/HNL momentum distributions.

**Impact:** HNL η/pT distributions slightly different from reality. May affect geometric acceptance by 10-20%.

**Fix:** Use proper form-factor-weighted decay (requires ExternalDecay plugin or post-hoc reweighting).

---

### 2A.5 — K_S/K_L production but no decay channels

**File:**
- `analysis_pbc/config/production_xsecs.py:112-113` (KS/KL have 0.5× K± cross-section)
- `production/pythia_production/main_hnl_production.cc` (no 310/130 channels configured)

**Issue:** The cross-section lookup supports K_S (310) and K_L (130), but Pythia doesn't configure their decay to HNL. These neutral kaons won't produce HNLs → inconsistency.

**Fix:** Either remove 310/130 from xsec lookup, or add the corresponding Pythia channels.

---

## §2B: Physics / Analysis

### 2B.1 — Majorana factor of 2 never applied

**File:**
- `production/pythia_production/main_hnl_production.cc:51-55` (documents factor 2 for Dirac)
- No application anywhere in analysis

**Issue:** Documentation says "For Dirac HNL interpretation (N ≠ N̄), multiply final yields by factor 2". But:
1. For Majorana HNL: this factor shouldn't be applied (currently correct)
2. For Dirac HNL: factor 2 should be applied somewhere

The code doesn't distinguish Majorana vs Dirac, and PBC benchmarks assume Majorana. This is fine for Majorana-only, but if users want Dirac interpretation, there's no mechanism.

**Fix:** Add `--majorana/--dirac` flag to analysis that applies ×2 factor if Dirac.

---

### 2B.2 — No detector/reconstruction efficiency

**File:** `analysis_pbc/limits/u2_limit_calculator.py:72-120`

**Issue:** The signal formula is:
```
N_sig = L × σ × BR × ε_geom × P_decay
```
Missing: reconstruction efficiency ε_reco (track finding, vertex reconstruction), trigger efficiency ε_trigger, material effects.

**Impact:** Sensitivity is **over-estimated** by factor of ~2-5 depending on detector assumptions.

**Fix:** Add configurable ε_reco factor. CMS drainage gallery has no tracking, so might assume ε=1 for decay-in-flight, but this should be documented.

---

### 2B.3 — No background model

**File:** Throughout analysis

**Issue:** Analysis assumes B=0 (zero background). Real limit setting requires:
1. Background estimation from data sidebands or MC
2. Statistical treatment (e.g., CLs method, not simple N > 2.996)

**Impact:** Limits are "projected sensitivity" not "expected limits". This is fine for physics case studies but shouldn't be called "exclusion limits" without caveats.

**Fix:** At minimum, document that these are B=0 projections. Better: add systematic uncertainty band.

---

### 2B.4 — W/Z BR formulas are approximations

**File:** `analysis_pbc/models/hnl_model_hnlcalc.py:266-302`

**Issue:** The W → ℓN and Z → νN branching ratios use simplified formulas:
```python
br_W = |U|² × BR_SM × (1 - r²)² × (1 + r²)
```
This is leading-order. Missing: NLO QCD corrections (~10%), finite-width effects, interference terms.

**Impact:** ~10-15% systematic uncertainty on EW production channel.

**Fix:** Compare with MadGraph cross-section and add correction factor if needed.

---

### 2B.5 — σ_K very approximate

**File:** `analysis_pbc/config/production_xsecs.py:67`

**Issue:** `SIGMA_KAON_PB = 5.0 * 1e10` is marked "very approximate, soft QCD" in comments. Kaon production cross-section varies by ~factor 2 depending on pT cut and rapidity region.

**Impact:** K-regime limits (m < 0.5 GeV) have ~factor 2 systematic uncertainty from cross-section alone.

**Fix:** Use measured K⁺ production cross-section from ALICE/CMS with explicit kinematic acceptance.

---

## Top 5 Issues

| # | Tag | File(s) | Issue | Minimal Fix |
|---|-----|---------|-------|-------------|
| **1** | [PHYSICS] | `hnl_model_hnlcalc.py:production_brs()` | **fromTau events get BR=0** — tau (PDG 15) missing from production_brs(), so signal=0 even when included | Add tau decay BR formula to production_brs() like W/Z |
| **2** | [BOTH] | `lhe_to_csv.py:196-207`, `production_brs()` | **parent_pdg=0 causes signal loss** — off-shell W/Z events get BR=0 (warnings logged but can be missed) | Force parent_pdg=24 for all EW events, or add fallback in BR lookup |
| **3** | [CODE] | `run_serial.py:80-82` | **Race condition in geometry cache** — parallel workers can corrupt cache files | Use atomic write (tempfile + rename) |
| **4** | [PHYSICS] | `main_hnl_production.cc:424-436` | **τ → πN only** — missing ρ channel biases acceptance | Add τ → ρN (PDG 213) channel |
| **5** | [PHYSICS] | `u2_limit_calculator.py:72` | **No ε_reco or background** — limits are over-optimistic projections | Document as "B=0 projected sensitivity"; add ε_reco parameter |

---

## Notes

- **AGENTS.md vs Code:** The combine step works as documented. The bug is in `production_brs()` missing tau BR, causing fromTau events to contribute zero signal despite being correctly merged.
- **Ground truth:** This review treats the code as ground truth per instructions.
- **Severity:** Issues ranked by impact on physics results, not code quality.

---

## Revision History

**v2 (same day):** Corrections after re-review:
- §2A.1: Originally claimed `combine_production_channels.py` skips fromTau. **Corrected:** combine step works correctly; the bug is `production_brs()` missing tau BR.
- §1A.1: Changed "silent signal loss" to "signal loss (logged)" — warnings exist but are easily missed.
- §1A.5: Reduced severity to LOW — this is a safety check, not a likely failure mode.
- §1B.3: Corrected — weight sanity check works correctly (ERROR case returns None, doesn't continue).
- §1B.5: Reduced to MINOR — tau xlim is actually reasonable.
