# Analysis Stage Validation

**Stage:** Limit Calculation with HNLCalc
**Language:** Python 3
**Input:** CSV files from production stage
**Output:** Exclusion limits |U|²_min, |U|²_max

---

## 📋 Purpose

Validate that the Python analysis correctly computes HNL sensitivity using:
- **Per-parent counting** (most critical decision)
- External cross-sections (no double-counting)
- HNLCalc for theory predictions
- Proper geometry and lifetime calculations

---

## 📁 Files in This Folder (7 files)

### 1. **`MULTI_HNL_METHODOLOGY.md`** ⭐ MOST CRITICAL

Explains per-parent vs per-event counting.

**Key concept:**
```python
# CORRECT: Per-parent
N_sig = Σ_parents [L × σ_parent × BR × ε_geom]

# WRONG: Per-event
N_sig = L × σ_??? × P_event
```

**Impact if wrong:** ~50% sensitivity loss

---

### 2. **`VALIDATION.md`**

Validation report confirming correctness.

**Expected benchmark:** 2.6 GeV muon → |U_μ|² ∈ [6.9×10⁻⁹, 2.4×10⁻⁵]

---

### 3. **`u2_limit_calculator.py`** ⭐ CORE

Main analysis driver.

**Lines ~70-228:** `expected_signal_events()` function
- Groups by `parent_id`
- Per-parent loop
- Uses HNLCalc for BR and cτ₀

---

### 4. **`per_parent_efficiency.py`**

Geometry: ray-tracing and boost factors.

**Outputs per HNL:**
- `hits_tube`: Intersects detector?
- `entry_distance`: Distance to entry [m]
- `path_length`: Path inside [m]
- `beta_gamma`: Boost factor p/m

---

### 5. **`hnl_model_hnlcalc.py`**

HNLCalc wrapper.

**Provides:**
- `ctau0_m`: Proper lifetime [m]
- `production_brs()`: {pdg: BR}

---

### 6. **`production_xsecs.py`** ⭐ CRITICAL

External cross-sections (NOT from Pythia).

**Key values:**
```python
σ(B⁰) = 4×10⁸ pb
σ(B⁺) = 4×10⁸ pb
σ(Bs) = 1×10⁸ pb
```

---

### 7. **`test_26gev_muon.py`** + **`test_closure.py`**

Benchmark and closure tests.

---

## 🔍 Critical Checks

### 1. Per-Parent Counting

**Location:** `u2_limit_calculator.py:~186-228`

```python
for pid in unique_parents:
    mask = (parent_id == pid)
    ε_parent = weighted_average(P_decay[mask])
    N_sig += L × σ_parent × BR × ε_parent
```

**Red flag:** `1 - np.prod(1-P_i)` → per-event (WRONG)

---

### 2. Cross-Section Source

**Must be:** External from literature (`production_xsecs.py`)
**Not:** From Pythia `.sigmaGen()` (double-counts!)

---

### 3. Boost Factors

```python
beta_gamma = momentum / mass
lam = beta_gamma * ctau0_m
```

---

### 4. Decay Probability

```python
P = exp(-entry_distance/lam) × (1 - exp(-path_length/lam))
```

Both terms required!

---

## 🧪 Tests

```bash
cd ../../analysis_pbc_test

# Closure test
conda run -n llpatcolliders python tests/closure_anubis/test_expected_signal_events_kernel.py

# Benchmark
conda run -n llpatcolliders python tests/test_26gev_muon.py
```

**Expected:** Tests pass, limits at [10⁻⁹, 10⁻⁵]

---

## ✅ Checklist

- [ ] `MULTI_HNL_METHODOLOGY.md` justifies per-parent
- [ ] `u2_limit_calculator.py` implements per-parent loop
- [ ] `production_xsecs.py` provides external σ
- [ ] Boost: `beta_gamma = p/m`
- [ ] Decay: Both survival + decay terms
- [ ] Tests pass

---

**Validation time:** 20-30 min
**Critical:** Per-parent counting + external cross-sections
