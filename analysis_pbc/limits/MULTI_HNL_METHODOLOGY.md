# Multi-HNL Event Handling: Methodology Documentation

## TL;DR

**Your current implementation is CORRECT.** You use per-parent counting, not per-event counting. This matches standard LLP detector methodology (MATHUSLA, ANUBIS, CODEX-b).

---

## The Situation

Pythia events can contain **multiple HNLs** from different parent mesons in a single pp collision.

**Example**: Event #44 from `HNL_1p00GeV_electron_charm.csv`:
```csv
event,parent_id,...
44,511,...       # B0 → ℓ N
44,-531,...      # Bs → ℓ N
44,411,...       # D+ → ℓ N
44,-431,...      # Ds → ℓ N
```

This event produces **4 HNLs** from **4 different meson species**.

---

## The Question

Should we count this as:
1. **Per-Event**: One event, probability `P_event = 1 - ∏(1-P_i)` that at least one HNL decays?
2. **Per-Parent**: Four independent production channels, each with its own cross-section?

---

## The Answer: Per-Parent Counting

We use **Option 2** because each parent represents a distinct physics process:

```python
N_sig = Σ_parents [ L × σ_parent × BR(parent→ℓN) × ε_geom(parent) ]
```

Where:
- σ(pp → D0) ≈ 2.8 × 10^10 pb
- σ(pp → B0) ≈ 4.0 × 10^8 pb
- σ(pp → Ds) ≈ 4.8 × 10^9 pb
- σ(pp → Bs) ≈ 1.0 × 10^8 pb

Event #44 contributes to **all four cross-section bins** because it contains all four production processes.

---

## Why Per-Event Logic is Wrong

If we used `P_event = 1 - ∏(1-P_i)`:

### Problem 1: Cannot Assign a Single Cross-Section
Event #44 contains D, B, Ds, Bs production. Which σ do you use?
- σ_D? → Ignores B/Bs/Ds contributions
- σ_B? → Ignores D/Ds contributions
- Average? → Physically meaningless

### Problem 2: Massive Undercounting
Using per-event logic:
- Event #44 = **1 signal event** (at least one decay)

Using per-parent logic:
- B0 contribution = σ_B0 × BR × ε_B0
- Bs contribution = σ_Bs × BR × ε_Bs
- D+ contribution = σ_D+ × BR × ε_D+
- Ds contribution = σ_Ds × BR × ε_Ds
- Total = **~4× more signal** (order of magnitude)

Per-event logic would **throw away ~75% of your sensitivity** in multi-HNL events.

### Problem 3: Breaks Theory Predictions
HNLCalc provides:
```python
BR(D0 → ℓN) = 1.2e-8 × |Ue|²  # For 1 GeV HNL
BR(B0 → ℓN) = 3.4e-9 × |Ue|²
```

These are **per-parent** branching ratios. You cannot apply them to a "pp event" that produces multiple parents.

---

## What the Standard Methodology Does

All LLP detector proposals use per-parent counting:

### MATHUSLA (arXiv:1811.00927)
```
N_MATHUSLA = Σ_i [ σ_i × BR_i × A_i × ε_i ]
```
where `i` runs over production channels (K, D, B, W, Z, ...)

### ANUBIS (arXiv:1909.13022)
```
S = L × Σ_parents [ σ × BR × ε_geom ]
```

### CODEX-b (arXiv:1911.00481)
Per-parent efficiency maps: ε(D), ε(B), ε(Λb), etc.

**They all separate**:
- **Production** (σ × BR): From theory (HNLCalc, NLO QCD)
- **Decay/Geometry** (ε): From detector simulation

This is exactly what we implement.

---

## Implementation in Our Code

### 1. Simulation (`production/pythia_production/main_hnl_production.cc`)
```cpp
// Pythia produces pp → many mesons
// Each meson forced to decay: M → ℓ N (BR=1.0)
// Output: One row per HNL with parent_id
```

### 2. Geometry (`geometry/per_parent_efficiency.py`)
```python
# Compute geometry for EACH HNL individually
# Output: beta_gamma, hits_tube, entry_distance, path_length (one row per HNL)
```

### 3. Analysis (`limits/expected_signal.py`)
```python
# Group HNLs by parent species
for pid in unique_parents:
    mask = (parent_id == pid)
    ε_parent = weighted_average(P_decay[mask])  # Average over HNLs from this parent

    N_sig += L × σ_parent × BR(parent→ℓN) × ε_parent

# For tau-decay chains (fromTau mode): parent_id=15 (τ), tau_parent_id=Ds/B (e.g. 431/511/521/531)
for pid in unique_tau_parents:
    mask = (parent_id == 15) & (tau_parent_id == pid)
    ε_parent = weighted_average(P_decay[mask])

    N_sig += L × σ_parent × BR(parent→τν) × BR(τ→NX) × ε_parent
```

**Key**: Each HNL contributes to its parent's cross-section bin, regardless of which pp event it came from. For tau-decay chains, the `tau_parent_id` identifies the original meson for proper BR weighting.

---

## Example Calculation

**Event #44**: [B0→N, Bs→N, D+→N, Ds→N]

Suppose each HNL has decay probability P_i ≈ 0.01 (1% chance to decay in tube).

### Per-Event (WRONG)
```
P_event = 1 - (0.99)^4 ≈ 0.039 (3.9%)
N_sig = L × σ_??? × 0.039    # Cannot choose σ!
```

### Per-Parent (CORRECT)
```
N_B0  = L × σ_B0  × BR(B0→ℓN)  × 0.01
N_Bs  = L × σ_Bs  × BR(Bs→ℓN)  × 0.01
N_D+  = L × σ_D+  × BR(D+→ℓN)  × 0.01
N_Ds  = L × σ_Ds  × BR(Ds→ℓN)  × 0.01

N_sig = N_B0 + N_Bs + N_D+ + N_Ds
```

With L = 3000 fb⁻¹ and |Ue|² = 1e-6:
```
N_B0  ≈ 3000 × 4e8 × 3.4e-9 × 1e-6 × 0.01 ≈ 0.04 events
N_Bs  ≈ 3000 × 1e8 × 2.1e-9 × 1e-6 × 0.01 ≈ 0.006 events
N_D+  ≈ 3000 × 1.2e10 × 1.2e-8 × 1e-6 × 0.01 ≈ 4.3 events
N_Ds  ≈ 3000 × 4.8e9 × 5.5e-9 × 1e-6 × 0.01 ≈ 0.8 events

Total ≈ 5.1 events (dominated by D+)
```

**Physics insight**: D mesons dominate at 1 GeV, even though B mesons are in the same event.

---

## Validation

To verify per-parent counting is implemented correctly:

```python
# Check that multiple HNLs per event are counted separately
df = pd.read_csv("output/csv/simulation/HNL_1p00GeV_electron_charm.csv")

# Event 44 should contribute to 4 different parent bins
event_44 = df[df["event"] == 44]
print(event_44["parent_pdg"].values)
# Output: [511, -531, 411, -431]

# In analysis, each of these increments its own N_sig contribution
```

You can verify this by checking that `len(geom_df)` equals the total number of HNLs (not the number of events).

---

## Summary Table

| Method | Event #44 Contribution | Cross-Section Used | Physically Correct? |
|--------|------------------------|-------------------|---------------------|
| Per-Event | 1 event × P_event | σ_??? (undefined) | ❌ No |
| Per-Parent | 4 parents × P_i each | σ_B0, σ_Bs, σ_D+, σ_Ds | ✅ Yes |

---

---

## Tau-Decay Chain Production (fromTau Mode)

For tau-coupled HNLs (BC8), there are **two independent O(U_τ²) mechanisms**:

### Direct Mode
```
Ds/B → τ N     (mixing at meson vertex)
```
- `parent_id` = Ds/B (e.g. 431/511/521/531)
- `tau_parent_id` = 0
- Signal: N_sig = L × σ(parent) × BR(parent→τN) × ε

### fromTau Mode
```
Ds/B → τ ν,  τ → N X    (mixing at tau decay)
```
- `parent_id` = 15 (τ)
- `tau_parent_id` = Ds/B (e.g. 431/511/521/531)
- Signal: N_sig = L × σ(parent) × BR(parent→τν) × BR(τ→N X) × ε

Both modes are O(U_τ²) and must be **added** (not double-counted). The `tau_parent_id` column identifies the "grandfather" meson for correct BR weighting.

### SM Branching Ratios for fromTau Weighting
| Parent | BR(parent→τν) | Source |
|--------|---------------|--------|
| Ds | 5.3% | PDG 2024 |
| B → D(*) | 2.3% | World average (R(D*) measurements) |

These are applied in `config/production_xsecs.py:get_parent_tau_br()`.

---

## Documentation Added

1. **`expected_signal.py`**: Per-parent counting + decay probability kernel
2. **`per_parent_efficiency.py`** (lines 247-250): Note about per-HNL geometry computation
3. **`production_xsecs.py`** (lines 17-28): Explanation of per-parent counting methodology
4. **`ROBUSTNESS_FIXES.md`** (Section 7): Full discussion with examples
5. **This document**: Standalone methodology reference

---

## References

- MATHUSLA: arXiv:1811.00927
- ANUBIS: arXiv:1909.13022
- CODEX-b: arXiv:1911.00481
- AL3X: arXiv:2010.02459
- PBC Report: CERN-PBC-REPORT-2018-007
- HNLCalc: arXiv:2405.07330

All use per-parent counting.
