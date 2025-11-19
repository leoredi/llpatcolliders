# HNL Coupling Limit Analysis - Technical Documentation

**Date**: 2025-11-19
**Purpose**: Complete documentation of the BR vs lifetime → |U_ℓ|² vs mass conversion for HNL searches
**Status**: ⚠️ RESULTS UNDER REVIEW - Coupling limits appear too weak (0.02-0.6 vs expected ~10⁻⁵-10⁻⁷)

## Executive Summary

This document details the complete physics analysis pipeline for converting detector-level exclusion limits (branching ratio vs lifetime) to the standard HNL parameter space (mixing parameter |U_ℓ|² vs mass).

**Current Issue**: The calculated coupling limits are orders of magnitude too weak:
- Observed: |U_μ|² ~ 0.02-0.6 for masses 15-39 GeV
- Expected: |U_μ|² ~ 10⁻⁵-10⁻⁷ (based on typical LLP detector sensitivities)

This suggests a bug in either:
1. The physics formulas (BR-lifetime-coupling relations)
2. The numerical conversion algorithm
3. The detector sensitivity calculation (upstream)

---

## 1. Physics Background

### 1.1 Heavy Neutral Leptons (HNLs)

HNLs are hypothetical particles that mix with Standard Model neutrinos through mixing parameter |U_ℓ|², where ℓ = e, μ, τ.

**Key relations:**
- Production: BR(W → ℓ N) ∝ |U_ℓ|²
- Decay: Γ_N ∝ |U_ℓ|² (lifetime τ_N ∝ 1/|U_ℓ|²)

### 1.2 Detector-Level Exclusion

The detector analysis (`decayProbPerEvent.py`) calculates:
- For each lifetime τ: What branching ratio BR(W → ℓ N) can we exclude?
- Output: BR_limit(τ) - the minimum BR we can exclude at each lifetime

**Method:**
```
For lifetime τ:
  1. Calculate decay probability P_decay(τ) for each particle reaching detector
  2. Calculate event-level probability (accounting for multiple particles)
  3. Find BR such that expected signal = 3 events (95% CL)
  4. BR_limit = 3 / (N_events × ε_geom × P_decay(τ))
```

### 1.3 Coupling Conversion

To convert to |U_ℓ|² limits, we use the physical relations:
- BR(|U|², m) - production depends on coupling and mass
- τ(|U|², m) - lifetime depends on coupling and mass

For each coupling value, we compute:
- Expected BR = BR(|U|², m)
- Expected lifetime = τ(|U|², m)

The coupling is excluded if: **BR(|U|², m) > BR_limit(τ(|U|², m))**

---

## 2. Implementation Details

### 2.1 Production Branching Ratio

**Location**: `hnl_coupling_limit.py`, function `br_w_to_hnl()`

```python
def br_w_to_hnl(coupling_sq, m_N, lepton='mu'):
    """
    Calculate BR(W → ℓ N) as a function of |U_ℓ|^2 and m_N

    Formula:
    BR(W → ℓ N) ≈ (1/9) × |U_ℓ|^2 × f(m_N)

    Where:
    - 1/9 ≈ BR(W → ℓ ν) for each lepton family
    - f(m_N) is the phase space factor
    """
    m_lepton = M_MU if lepton == 'mu' else M_TAU

    # Phase space factor
    f = phase_space_factor(m_N, M_W, m_lepton)

    # BR(W → ℓ ν) ≈ 1/9 for each generation
    br_w_lnu = 1.0 / 9.0

    # BR(W → ℓ N) = |U_ℓ|^2 × f × BR(W → ℓ ν)
    br = coupling_sq * f * br_w_lnu

    return br
```

**Phase space factor**:
```python
def phase_space_factor(m_N, m_W=M_W, m_lepton=M_MU):
    """
    Phase space factor for W → ℓ N decay

    Formula:
    f(m_N, m_W, m_ℓ) = [(1 - x_N^2)^2 - x_ℓ^2(1 + x_N^2)]

    where x_N = m_N/m_W and x_ℓ = m_ℓ/m_W
    """
    x_N = m_N / m_W
    x_l = m_lepton / m_W

    # Kinematic threshold
    if m_N + m_lepton >= m_W:
        return 0.0

    # Phase space factor (proportional to decay width)
    f = (1 - x_N**2)**2 - x_l**2 * (1 + x_N**2)

    return max(f, 0.0)
```

**Numerical check** for m_N = 15 GeV, m_μ = 0.106 GeV, m_W = 80.4 GeV:
```
x_N = 15 / 80.4 = 0.187
x_μ = 0.106 / 80.4 = 0.00132
f = (1 - 0.187²)² - 0.00132² × (1 + 0.187²)
  = (1 - 0.035)² - 1.74e-6 × 1.035
  = 0.931 - 1.8e-6
  ≈ 0.931
```

**⚠️ ISSUE #1**: Is the factor of 1/9 correct?

The formula assumes BR(W → ℓ N) ≈ |U_ℓ|² × f × (1/9). This might be incorrect:
- Standard Model: BR(W → ℓ ν) ≈ 0.108 (10.8%, not 11.1%)
- The normalization might need to include additional factors (CKM mixing, etc.)

---

### 2.2 HNL Lifetime

**Location**: `hnl_coupling_limit.py`, function `hnl_lifetime()`

```python
def hnl_lifetime(coupling_sq, m_N):
    """
    Calculate HNL lifetime as a function of |U_ℓ|^2 and m_N

    Formula:
    Γ_N ≈ C × |U_ℓ|^2 × G_F^2 × m_N^5

    Calibrated using CMS 2024 data: for m=10 GeV, |U|^2 = 5×10^-7 gives cτ = 17 mm

    This gives C ≈ 1.7×10^-3
    """
    C = 1.7e-3

    # Γ_N in GeV
    gamma_N = coupling_sq * C * G_F**2 * m_N**5

    # Lifetime τ = ℏ / Γ
    if gamma_N > 0:
        tau = HBAR_GEV_S / gamma_N  # in seconds
    else:
        tau = np.inf

    return tau
```

**Constants used**:
```python
HBAR_GEV_S = 6.582119569e-25  # ℏ in GeV·s
G_F = 1.1663787e-5           # Fermi constant in GeV^-2
```

**Calibration check** (CMS 2024: m=10 GeV, |U|²=5×10⁻⁷, cτ=17 mm):
```
Γ_N = 5×10⁻⁷ × 1.7×10⁻³ × (1.166×10⁻⁵)² × 10⁵
    = 5×10⁻⁷ × 1.7×10⁻³ × 1.36×10⁻¹⁰ × 10⁵
    = 5×10⁻⁷ × 1.7×10⁻³ × 1.36×10⁻⁵
    = 1.156×10⁻¹⁴ GeV

τ = ℏ/Γ = 6.582×10⁻²⁵ / 1.156×10⁻¹⁴
  = 5.69×10⁻¹¹ s

cτ = c × τ = 3×10⁸ × 5.69×10⁻¹¹
   = 0.01707 m = 17.07 mm ✓
```

Calibration checks out!

**⚠️ ISSUE #2**: Is the calibration constant C universal?

The value C ≈ 1.7×10⁻³ was calibrated for m=10 GeV. Does it hold for all masses 15-71 GeV?
- Decay width formulas often have mass-dependent form factors
- May need mass-dependent C(m) instead of constant

---

### 2.3 Coupling Limit Calculation

**Location**: `hnl_coupling_limit.py`, function `compute_coupling_limit()`

```python
def compute_coupling_limit(mass_gev, br_limits, lifetimes, lepton='mu'):
    """
    Convert BR vs lifetime exclusion to coupling^2 limit

    Algorithm:
    1. Create interpolator: BR_limit(τ)
    2. For each trial coupling |U|²:
       a. Calculate expected BR(|U|², m)
       b. Calculate expected τ(|U|², m)
       c. Check if BR(|U|², m) > BR_limit(τ(|U|², m))
    3. Find minimum excluded coupling
    """
    # Create interpolation function for BR_limit(τ)
    log_tau = np.log10(lifetimes)
    log_br = np.log10(br_limits)

    # Remove any infinite or nan values
    valid = np.isfinite(log_tau) & np.isfinite(log_br)

    br_limit_interp = interp1d(log_tau[valid], log_br[valid],
                               kind='linear', bounds_error=False,
                               fill_value=(log_br[valid][0], log_br[valid][-1]))

    # Scan over coupling values
    log_coupling_sq_range = np.linspace(-12, 0, 1000)  # |U|^2 from 10^-12 to 1
    coupling_sq_range = 10**log_coupling_sq_range

    excluded = []

    for coupling_sq in coupling_sq_range:
        # Calculate expected BR and lifetime for this coupling
        br_expected = br_w_to_hnl(coupling_sq, mass_gev, lepton)
        tau_expected = hnl_lifetime(coupling_sq, mass_gev)

        # Get BR limit at this lifetime
        log_tau_exp = np.log10(tau_expected)
        log_br_lim = br_limit_interp(log_tau_exp)
        br_lim = 10**log_br_lim

        # Check if excluded
        is_excluded = (br_expected > br_lim)
        excluded.append(is_excluded)

    # Find the boundary (smallest coupling that's excluded)
    excluded = np.array(excluded)
    if np.any(excluded):
        excluded_indices = np.where(excluded)[0]
        if len(excluded_indices) > 0:
            boundary_idx = excluded_indices[0]
            coupling_sq_limit = coupling_sq_range[boundary_idx]
            return coupling_sq_limit

    return np.nan
```

**⚠️ ISSUE #3**: Algorithm logic check

The algorithm finds the **first** (minimum) coupling that gets excluded. But is this correct?

For increasing |U|²:
- BR(|U|²) increases ∝ |U|²
- τ(|U|²) decreases ∝ 1/|U|²

As |U|² increases:
1. Particle decays faster (shorter lifetime)
2. Production rate increases
3. BR_limit(τ) changes (depends on detector geometry)

**Question**: Is BR_limit(τ) monotonic? If not, the exclusion region might not be contiguous!

---

## 3. Input Data Analysis

### 3.1 BR vs Lifetime Exclusions

**Source**: `output/csv/hnlLL_m{mass}GeVLLP_exclusion_data.csv`

These files contain the detector-level exclusions calculated by `decayProbPerEvent.py`.

**Example data for m=15 GeV**:
```python
# Read the file
df = pd.read_csv('output/csv/hnlLL_m15GeVLLP_exclusion_data.csv')

# Columns:
# - lifetime_s: HNL lifetime in seconds
# - ctau_m: Decay length cτ in meters
# - BR_limit: Branching ratio exclusion limit
# - mean_event_decay_prob: Mean event-level decay probability
# - mean_single_particle_decay_prob: Mean single-particle decay probability
```

**Data summary for all masses**:
```
m15: BR_limit range = 9.43e-09 to 2.11e-03
m23: BR_limit range = 9.61e-09 to 2.74e-03
m31: BR_limit range = 1.08e-08 to 4.54e-02
m39: BR_limit range = 1.14e-08 to 5.71e-03
m47: BR_limit range = 1.35e-08 to 5.12e+00  ⚠️ BR > 1 (unphysical)
m55: BR_limit range = 1.58e-08 to 7.15e+02  ⚠️ BR > 100
m63: BR_limit range = 1.95e-08 to 1.88e+01  ⚠️ BR > 10
m71: BR_limit range = 2.55e-08 to 4.71e+02  ⚠️ BR > 100
```

**⚠️ ISSUE #4**: BR limits > 1 at high masses

For masses ≥47 GeV, the BR limits exceed 1 (some exceed 100). This means:
- The detector has extremely poor sensitivity (can only exclude BR>100%)
- These limits are physically meaningless
- **Root cause**: Likely very low acceptance at high masses

**Critical question**: Are the upstream BR_limit calculations correct?

The BR_limit formula in `decayProbPerEvent.py`:
```python
# N_signal = BR × N_W_produced × ε_geom × P_decay(τ)
# For 95% CL: N_signal = 3
# Therefore: BR_limit = 3 / (N_W × ε_geom × P_decay)

BR_limit = 3.0 / (n_events * mean_event_decay_prob)
```

If `mean_event_decay_prob` is very small (10⁻³ or 10⁻⁴), then BR_limit becomes large!

---

## 4. Numerical Results

### 4.1 Phase Space Factors

For muon-coupled HNLs (m_μ = 0.106 GeV, m_W = 80.4 GeV):

```
Mass (GeV) | x_N    | f(m_N) | BR(|U|²=1)
-----------|--------|--------|------------
15         | 0.187  | 0.931  | 0.103
23         | 0.286  | 0.836  | 0.093
31         | 0.386  | 0.698  | 0.078
39         | 0.485  | 0.529  | 0.059
47         | 0.585  | 0.346  | 0.038
55         | 0.684  | 0.173  | 0.019
63         | 0.784  | 0.034  | 0.004
71         | 0.883  | 0.003  | 0.0003
```

Note: Phase space suppression becomes severe above 55 GeV!

### 4.2 Lifetime Calculations

For coupling |U_μ|² = 10⁻⁵:

```
Mass (GeV) | m^5 (GeV^5) | Γ_N (GeV)   | τ (s)      | cτ (m)
-----------|-------------|-------------|------------|--------
15         | 7.59×10⁵    | 1.35×10⁻¹⁶  | 4.87×10⁻⁹  | 1.46
23         | 6.44×10⁶    | 1.14×10⁻¹⁵  | 5.76×10⁻¹⁰ | 0.17
31         | 2.86×10⁷    | 5.07×10⁻¹⁵  | 1.30×10⁻¹⁰ | 0.039
39         | 9.03×10⁷    | 1.60×10⁻¹⁴  | 4.11×10⁻¹¹ | 0.012
47         | 2.29×10⁸    | 4.06×10⁻¹⁴  | 1.62×10⁻¹¹ | 0.0049
55         | 5.03×10⁸    | 8.93×10⁻¹⁴  | 7.37×10⁻¹² | 0.0022
63         | 9.93×10⁸    | 1.76×10⁻¹³  | 3.74×10⁻¹² | 0.0011
71         | 1.80×10⁹    | 3.20×10⁻¹³  | 2.06×10⁻¹² | 0.0006
```

Strong mass dependence due to m⁵!

### 4.3 Calculated Coupling Limits

**Results**:
```
Mass (GeV) | |U_μ|² limit | Status
-----------|--------------|--------
15         | 2.08×10⁻²    | Valid
23         | 2.98×10⁻²    | Valid
31         | 5.75×10⁻¹    | Valid
39         | 9.02×10⁻²    | Valid
47         | NaN          | No exclusion (BR_limit > 1)
55         | NaN          | No exclusion (BR_limit >> 1)
63         | NaN          | No exclusion (BR_limit >> 1)
71         | NaN          | No exclusion (BR_limit >> 1)
```

**⚠️ MAJOR ISSUE**: These limits are **far too weak**!

Expected LLP detector sensitivity: |U|² ~ 10⁻⁵ to 10⁻⁷
Observed: |U|² ~ 10⁻² to 10⁻¹

**Factor of ~10³ to 10⁵ discrepancy!**

---

## 5. Diagnostic Checks

### 5.1 Check: Forward Calculation at m=15 GeV

Let's trace through for |U_μ|² = 0.02 (the claimed limit):

**Production:**
```
BR(|U|²=0.02, m=15) = 0.02 × 0.931 × (1/9)
                     = 0.02 × 0.103
                     = 2.07×10⁻³
```

**Lifetime:**
```
Γ_N = 0.02 × 1.7×10⁻³ × (1.166×10⁻⁵)² × (15)⁵
    = 0.02 × 1.7×10⁻³ × 1.36×10⁻¹⁰ × 7.59×10⁵
    = 3.49×10⁻¹⁵ GeV

τ = 6.582×10⁻²⁵ / 3.49×10⁻¹⁵ = 1.89×10⁻¹⁰ s

cτ = 0.057 m = 5.7 cm
```

**Question**: At τ = 1.89×10⁻¹⁰ s, what is BR_limit from the exclusion file?

Need to check: Does BR_limit(τ=1.89×10⁻¹⁰) ≈ 2.07×10⁻³?

If yes: Algorithm is correct, physics formulas might be wrong
If no: Algorithm has a bug

### 5.2 Check: Detector Acceptance

The BR_limit depends on detector acceptance. High BR_limit means low acceptance.

For m=47 GeV with BR_limit ~ 5:
```
BR_limit = 3 / (N_events × P_decay × ε_geom)

If BR_limit = 5 and N_events = 1M:
5 = 3 / (10⁶ × P_decay × ε_geom)
P_decay × ε_geom = 6×10⁻⁷
```

This is extremely small acceptance!

**Possible causes**:
1. High-mass HNLs don't reach detector (wrong kinematics)
2. Detector geometry not sensitive to high-mass regime
3. Bug in acceptance calculation

---

## 6. Code Locations

### 6.1 Main Scripts

**Coupling conversion**: `hnl_coupling_limit.py`
- Lines 29-49: `phase_space_factor()`
- Lines 51-84: `br_w_to_hnl()`
- Lines 86-124: `hnl_lifetime()`
- Lines 160-224: `compute_coupling_limit()`
- Lines 226-271: `create_coupling_mass_plot()`

**Detector analysis**: `decayProbPerEvent.py`
- Lines ~300-400: Lifetime scanning loop
- Lines ~450-500: BR_limit calculation
- Lines ~545-556: Exclusion data export

### 6.2 Input Files

**Simulation data**: `output/csv/hnlLL_m{mass}GeVLLP.csv`
- Format: `event,id,pt,eta,phi,momentum,mass`
- 1M events per mass point

**Exclusion data**: `output/csv/hnlLL_m{mass}GeVLLP_exclusion_data.csv`
- Columns: `lifetime_s, ctau_m, BR_limit, mean_event_decay_prob, mean_single_particle_decay_prob`
- 20 lifetime points per mass

### 6.3 External Data

**Experimental limits**: `external/{MATHUSLA,CODEX,ANUBIS}.csv`
- Format varies
- Used for comparison plots

---

## 7. Hypotheses for Bug Location

### Hypothesis A: Production BR Formula Wrong

**Current**:
```
BR(W → ℓ N) = |U_ℓ|² × f(m_N) × (1/9)
```

**Issue**: Factor of 1/9 might be wrong or missing additional terms

**Test**: Compare with published HNL production formulas (e.g., arXiv:1805.08567, arXiv:1912.07622)

---

### Hypothesis B: Lifetime Formula Wrong

**Current**:
```
Γ_N = C × |U_ℓ|² × G_F² × m_N⁵
C = 1.7×10⁻³ (calibrated at m=10 GeV)
```

**Issues**:
1. C might not be universal (mass-dependent)
2. Missing form factors or phase space corrections
3. Calibration might be wrong

**Test**:
- Check published lifetime formulas
- Verify calibration against multiple experimental points
- Calculate C(m) from first principles

---

### Hypothesis C: Interpolation/Algorithm Bug

**Current approach**:
1. Scan |U|² from 10⁻¹² to 1
2. For each |U|²: compute BR(|U|²) and τ(|U|²)
3. Check if BR(|U|²) > BR_limit(τ(|U|²))
4. Return minimum excluded |U|²

**Potential issues**:
1. BR_limit(τ) interpolation might extrapolate incorrectly
2. Exclusion region might not be contiguous
3. Scan resolution (1000 points) might miss fine structure

**Test**:
- Plot BR(|U|²) vs τ(|U|²) parametrically
- Overlay with BR_limit(τ) curve
- Visually check intersection

---

### Hypothesis D: Detector Sensitivity Wrong (Upstream)

**Issue**: BR_limit values might be calculated incorrectly in `decayProbPerEvent.py`

**Formula**:
```python
BR_limit = 3.0 / (n_events * mean_event_decay_prob)
```

**Questions**:
1. Is the factor of 3 correct for 95% CL?
2. Should this be particle-level or event-level probability?
3. Is geometric acceptance included properly?
4. Are event weights handled correctly?

**Test**:
- Manually calculate expected signal for known (BR, τ, m)
- Compare with detector simulation
- Check if BR_limit formula is correct

---

## 8. Recommended Debugging Steps

### Step 1: Validate Physics Formulas

**Action**: Compare with authoritative sources
- [ ] Check production BR formula against arXiv:1805.08567
- [ ] Verify lifetime formula against PDG or HNL reviews
- [ ] Recalculate phase space factors independently
- [ ] Test calibration constant against multiple experiments

**Files to check**: `hnl_coupling_limit.py` lines 29-124

---

### Step 2: Trace Single Example

**Action**: Full manual calculation for m=15 GeV, |U|²=0.02

```python
# Production
m = 15
coupling_sq = 0.02
BR_expected = br_w_to_hnl(coupling_sq, m, 'mu')
print(f"BR_expected = {BR_expected:.3e}")

# Lifetime
tau_expected = hnl_lifetime(coupling_sq, m)
print(f"tau_expected = {tau_expected:.3e} s")
print(f"ctau = {tau_expected * 3e8:.3e} m")

# Load exclusion data
df = pd.read_csv('output/csv/hnlLL_m15GeVLLP_exclusion_data.csv')

# Find BR_limit at this lifetime
from scipy.interpolate import interp1d
interp = interp1d(np.log10(df.lifetime_s), np.log10(df.BR_limit))
log_br_lim = interp(np.log10(tau_expected))
BR_limit = 10**log_br_lim
print(f"BR_limit(tau={tau_expected:.3e}) = {BR_limit:.3e}")

# Check exclusion
print(f"Excluded? {BR_expected > BR_limit}")
print(f"Ratio: BR_expected/BR_limit = {BR_expected/BR_limit:.3f}")
```

**Expected outcome**: If algorithm is correct, ratio should be ≈1 at the limit

---

### Step 3: Visualize Exclusion Contours

**Action**: Create diagnostic plot

```python
import matplotlib.pyplot as plt

coupling_range = np.logspace(-12, 0, 1000)
mass = 15

BR_expected = [br_w_to_hnl(c, mass, 'mu') for c in coupling_range]
tau_expected = [hnl_lifetime(c, mass) for c in coupling_range]

# Load BR_limit(tau)
df = pd.read_csv(f'output/csv/hnlLL_m{mass}GeVLLP_exclusion_data.csv')

# Plot in (tau, BR) space
plt.figure(figsize=(10, 8))
plt.loglog(df.lifetime_s, df.BR_limit, 'r-', linewidth=2, label='BR_limit(τ)')
plt.loglog(tau_expected, BR_expected, 'b-', linewidth=2, label='BR(|U|²) vs τ(|U|²)')
plt.xlabel('Lifetime (s)')
plt.ylabel('Branching Ratio')
plt.legend()
plt.grid(True, alpha=0.3)
plt.title(f'Exclusion Check: m={mass} GeV')
plt.savefig(f'debug_exclusion_m{mass}.png')
```

**Expected**: Curves should intersect at |U|² = limit

---

### Step 4: Check Upstream BR_limit Calculation

**Action**: Audit `decayProbPerEvent.py`

**Key questions**:
1. Is `mean_event_decay_prob` calculated correctly?
2. Does it include all geometric acceptance factors?
3. Is the normalization to W production rate correct?
4. Are event weights handled properly?

**Test**: For a known (BR, τ, m), inject signal and check if BR_limit agrees

---

### Step 5: Literature Comparison

**Action**: Find published HNL limits

Look for:
- MATHUSLA projections: arXiv:1901.04346
- CODEX-b projections: arXiv:1911.00762
- Recent HNL searches at LHC

**Compare**: Do our coupling limits match published projections?

---

## 9. Expected Corrections

Based on typical LLP detector sensitivities, we expect:

**milliQan-like detector at HL-LHC:**
- Luminosity: 3000 fb⁻¹
- Distance: ~20m from IP
- Expected sensitivity: |U|² ~ 10⁻⁶ to 10⁻⁵ for m ~ 10-50 GeV

**Our results** (|U|² ~ 0.02-0.6) are **10,000× too weak!**

This strongly suggests either:
1. Physics formulas are wrong by factor ~10⁴
2. Detector acceptance is underestimated by factor ~10⁴
3. Combination of multiple errors

---

## 10. Summary

**What works**:
- ✓ Calibration check passes for m=10 GeV reference point
- ✓ Phase space factors look reasonable
- ✓ Algorithm logic appears sound
- ✓ Data files exist and are properly formatted

**What's wrong**:
- ✗ Coupling limits are 10³-10⁵ times weaker than expected
- ✗ High-mass points (47-71 GeV) have no sensitivity (BR_limit > 1)
- ✗ Overall sensitivity curve doesn't match expectations

**Most likely culprits** (in order):
1. **BR_limit calculation in `decayProbPerEvent.py`** - Missing normalization or wrong formula
2. **Production BR formula** - Factor of 1/9 or phase space factor incorrect
3. **Lifetime formula calibration** - C constant wrong or mass-dependent
4. **Algorithm implementation** - Interpolation or scan issues

**Next steps**:
1. Validate all physics formulas against literature
2. Trace through single example by hand
3. Check upstream BR_limit calculation
4. Create diagnostic visualization plots
5. Compare with published experimental projections

---

## Appendix: Complete Physics Formulas from Literature

### A.1 HNL Production (arXiv:1805.08567)

For W → ℓ N in the Minimal Heavy Neutrino model:

```
Γ(W → ℓ N) = (G_F² m_W³)/(6π) × |U_ℓ|² × λ^(1/2)(1, x_N², x_ℓ²) × [...]

where λ(a,b,c) = a² + b² + c² - 2ab - 2ac - 2bc
```

This should be compared with our simplified formula!

### A.2 HNL Decay Width (PDG Review)

For N → ℓ + hadrons (dominant at high mass):

```
Γ_N = (G_F² m_N⁵)/(192π³) × |U_ℓ|² × [sum over channels with form factors]
```

The constant C in our code (1.7×10⁻³) should match:
```
C = 1/(192π³) ≈ 1.7×10⁻⁵ (?)
```

**⚠️ DISCREPANCY**: Our C = 1.7×10⁻³ vs theoretical 1.7×10⁻⁵

This is a **factor of 100 error** in the decay width!

---

**END OF DOCUMENT**

*This document should be continuously updated as bugs are found and fixed.*
