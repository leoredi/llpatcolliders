# Defensive Programming and Diagnostics Guide

This document describes the defensive programming improvements in the HNL coupling limit analysis to prevent silent failures and data loss.

## 1. Comprehensive NaN Filtering (Lines 309-328)

### Problem
If `weight`, `beta_gamma`, `entry_distance`, or `path_length` contained NaNs (from corrupt CSVs, partial writes, or manual edits), they would propagate through physics calculations:

1. `lam = beta_gamma * ctau0_m` → NaN
2. `P_decay` → NaN
3. `eff_parent = np.sum(w * P) / w_sum` → NaN
4. Result: "No sensitivity" for all mass points with no clear error message

### Solution
Check ALL critical columns before physics calculations:
```python
cols_to_check = ["parent_id", "weight", "beta_gamma", "entry_distance", "path_length"]
mask_valid = np.ones(len(geom_df), dtype=bool)

for col in cols_to_check:
    if col in geom_df.columns:
        mask_valid &= geom_df[col].notna().to_numpy() & np.isfinite(geom_df[col].to_numpy())
    else:
        print(f"[WARN] m={mass_str} {flavour}: Missing required column '{col}'.")
        return None
```

### Benefits
- **Early detection**: Corrupt data caught before physics calculations
- **Clear logging**: "Dropping N rows with NaNs in geometry/weights"
- **Prevents silent failures**: No more mysterious "No sensitivity" results
- **Future-proof**: Protects against edge cases

---

## 2. PDG Coverage Diagnostics (Lines 186-228)

### Problem
If a parent PDG code appears in simulation but is:
- **Not in HNLCalc**: `BR_parent = 0.0` → event silently discarded
- **Not in production_xsecs.py**: `sigma_parent_pb = 0.0` → event silently discarded

These events are lost without warning, potentially causing incorrect physics results.

**Example**: PDG 310 (KS0) appears in CSV but HNLCalc doesn't model K→HNL, so all KS0 events are discarded.

### Solution
Track and log missing PDG codes:
```python
# Track diagnostics for missing PDG codes
missing_br_pdgs = []
missing_xsec_pdgs = []

for pid in unique_parents:
    BR_parent = br_per_parent.get(int(pid), 0.0)
    if BR_parent <= 0.0:
        missing_br_pdgs.append(int(pid))
        continue

    sigma_parent_pb = get_parent_sigma_pb(int(pid))
    if sigma_parent_pb <= 0.0:
        missing_xsec_pdgs.append(int(pid))
        continue

# Log warnings (only at first eps2 point to avoid spam)
if missing_br_pdgs and eps2 == 1e-12:
    n_lost = np.sum(np.isin(parent_id, missing_br_pdgs))
    print(f"[WARN] Mass {mass_GeV:.2f} GeV: {len(missing_br_pdgs)} parent PDG(s) have no HNLCalc BR: {missing_br_pdgs}")
    print(f"       → Discarding {n_lost} events (silent data loss)")
```

**Smart logging**: Only warn once per mass point (at eps2=1e-12), avoiding spam during 100-point scans.

### Diagnostic Tool
Created `diagnostic_pdg_coverage.py` to audit coverage across simulation, HNLCalc, and cross-sections.

**Usage**:
```bash
conda run -n llpatcolliders python limits/diagnostic_pdg_coverage.py
```

**Output**:
```
⚠️  WARNING: 1 PDG codes in CSV but NOT in HNLCalc:
   [310]  # KS0
   → These events will have BR=0.0 and contribute nothing to signal

PDG      Name       In CSV   HNLCalc    Has σ    Status
--------------------------------------------------------------------------------
310      KS0        ✓        ✗          ✓        ⚠️  LOST
411      D+         ✓        ✓          ✓        OK
421      D0         ✓        ✓          ✓        OK
```

### Benefits
- **Visibility**: No more silent data loss
- **Actionable**: Identifies exactly which PDG codes are missing
- **Preventive**: Run before analysis to catch coverage gaps
- **Maintainable**: Catches mismatches when updating HNLCalc or Pythia

---

## 3. Current Coverage Status

### ✅ Fully Covered
- 411 (D+), 421 (D0), 431 (Ds+), 4122 (Λc+)
- 511 (B0), 521 (B+), 531 (Bs0), 5122 (Λb0)

### ⚠️ Silent Data Loss
- **310 (KS0)**: Has cross-section but no HNLCalc BR
  - **Impact**: All KS0→HNL events discarded
  - **Action**: Verify if KS0 fraction is significant

### 📝 Not Simulated
- 15 (τ), 541 (Bc+), 4132 (Ξc+), 5232 (Σb), 5332 (Ωb)
- Expected — Pythia focuses on common mesons

---

## 4. Multi-HNL Event Handling

### The Question
Pythia events can contain **multiple HNLs** from different parents:
```csv
event,parent_id,...
44,511,...       # B0 → N
44,-531,...      # Bs → N
44,411,...       # D+ → N
44,-431,...      # Ds → N
```

Should we count this as one event or four production channels?

### The Answer: Per-Parent Counting

We use **per-parent counting**:
```python
N_sig = Σ_parents [ L × σ_parent × BR(parent→ℓN) × ε_geom(parent) ]
```

Each HNL is treated independently according to its parent's production cross-section.

### Why This is Correct

1. **Different parents have different cross-sections**:
   - σ(pp → D0) ≈ 2.8 × 10¹⁰ pb
   - σ(pp → B0) ≈ 4.0 × 10⁸ pb
   - Cannot assign single σ to multi-parent events

2. **Matches PBC methodology**:
   - MATHUSLA (arXiv:1811.00927)
   - ANUBIS (arXiv:1909.13022)
   - CODEX-b (arXiv:1911.00481)

3. **Physically transparent**:
   - Separates production (theory) from decay (geometry)
   - Allows tracking D vs B contributions
   - Compatible with HNLCalc BRs

**See `MULTI_HNL_METHODOLOGY.md` for detailed explanation with examples.**

---

## 5. Files Modified

1. **limits/u2_limit_calculator.py**
   - Lines 309-328: Comprehensive NaN filtering
   - Lines 186-228: PDG coverage diagnostics
   - Lines 77-118: Multi-HNL methodology docstring

2. **limits/diagnostic_pdg_coverage.py** (New)
   - Standalone diagnostic tool
   - Audits PDG coverage across simulation/HNLCalc/cross-sections

3. **limits/MULTI_HNL_METHODOLOGY.md** (New)
   - Comprehensive explanation of per-parent counting
   - Example calculations and comparisons

---

## 6. Recommendations

### Before Analysis
```bash
# Run diagnostic to check PDG coverage
conda run -n llpatcolliders python limits/diagnostic_pdg_coverage.py
```

### During Analysis
- Monitor logs for "[WARN] parent PDG(s) have no..." messages
- Check if missing PDG codes are significant (>1% of events)

### After Adding New Simulations
- Re-run diagnostic to verify coverage
- Update production_xsecs.py if new parents appear

---

## 7. Testing

### Test NaN filtering:
```python
# Manually corrupt a geometry CSV
df = pd.read_csv("output/csv/geometry/HNL_mass_1.0_muon_geom.csv")
df.loc[0:10, "beta_gamma"] = np.nan
df.to_csv("output/csv/geometry/HNL_mass_1.0_muon_geom.csv", index=False)

# Run analysis - should see:
# [INFO] m=1.0 muon: Dropping 11 rows with NaNs in geometry/weights.
```

### Test PDG diagnostics:
```bash
conda run -n llpatcolliders python limits/diagnostic_pdg_coverage.py
# Should report KS0 (310) as LOST
```

---

## Summary

These fixes transform the analysis from "silently failing" to "loudly diagnosing":

| Failure Mode | Before | After |
|--------------|--------|-------|
| NaN in geometry CSV | Silent NaN propagation → "No sensitivity" | Early filtering + clear warning |
| Unknown parent PDG | Silent discard (BR=0 or σ=0) | Warning logged with event count |
| Missing cross-section | Silent discard | Diagnostic tool catches gaps |
| Corrupt CSV | Mysterious failures | Clear error messages |

**Key principles**:
1. If data is being discarded, the user should know about it
2. If the methodology is non-obvious, document it (see MULTI_HNL_METHODOLOGY.md)
3. If an assumption is fragile, validate it at runtime
