# Meson-EW Production Channel Combination

## Issue Summary

In the transition region (4-8 GeV), HNLs can be produced from BOTH:
1. **Meson production** (Pythia): B/D → ℓN
   - Parent PDG codes: 511 (B⁰), 521 (B±), 421 (D⁰), 431 (Ds), etc.
   - Files: `HNL_5p0GeV_muon_beauty.csv`

2. **Electroweak production** (MadGraph): W/Z → ℓN
   - Parent PDG codes: 23 (Z), 24 (W±)
   - Files: `HNL_5p0GeV_muon_ew.csv`

## Is This Double-Counting?

Meson and EW production are different mechanisms and can be summed using per-parent counting:
- σ(pp → B) × BR(B → μN) ~ O(10⁸ pb) × O(|U|²) at 5 GeV
- σ(pp → W) × BR(W → μN) ~ O(10⁸ pb) × O(0.1|U|²) at 5 GeV

**Open point:** Tau cascade samples (`_fromTau`) are now included by default when combining. If you want to avoid overlap between direct W/Z → τN and τ → πN cascades, you may choose to drop `_fromTau` files before combining. Per-parent counting still separates parent PDGs, but the policy should be explicit.

## Current Status (Jan 2025)

- Pythia meson and MadGraph EW samples can both be present.
- `combine_production_channels.py` accepts `_direct` **and** `_fromTau` CSVs; no automatic filter is applied.
- Decide whether to keep or drop `_fromTau` when combining with EW tau samples; document the choice.

## Future Action Required

When Pythia production completes, follow these steps:

### 1. Check for Overlapping Files

```bash
cd output/csv/simulation
ls HNL_*GeV_muon_*.csv | grep -E "(4p|5p|6p|7p|8p)"
```

If you see BOTH beauty and ew files at the same mass:
```
HNL_5p0GeV_muon_beauty.csv
HNL_5p0GeV_muon_ew.csv
```

Then you need to combine them.

### 2. Combine Production Channels

```bash
cd analysis_pbc
python limits/combine_production_channels.py --dry-run  # Preview
python limits/combine_production_channels.py            # Execute
```

This will:
- Find all masses with multiple production files (including `_fromTau` unless you remove them first)
- Concatenate CSVs (preserving all parent PDG codes)
- Create unified files: `HNL_5p0GeV_muon_combined.csv`

### 3. Clean Up Original Files

```bash
# Archive regime-specific files
mkdir -p output/csv/simulation_backup
mv output/csv/simulation/HNL_*_beauty.csv simulation_backup/
mv output/csv/simulation/HNL_*_ew.csv simulation_backup/

# Move combined files to main directory
mv output/csv/simulation/combined/* output/csv/simulation/
```

### 4. Re-run Analysis

```bash
cd analysis_pbc
python limits/run_serial.py
```

The per-parent counting in `limits/expected_signal.py` will correctly sum:
```python
N_sig = Σ_parents [ L × σ(parent) × BR(parent→ℓN) × ε_geom(parent) ]
      = L × [σ(B) × BR(B→μN) × ε_B] + [σ(W) × BR(W→μN) × ε_W]
```

## Why This Approach is Correct

The analysis uses **per-parent counting** (not per-event):

```python
# From limits/expected_signal.py (per-parent counting)
for pid in unique_parents:
    BR_parent = br_per_parent.get(int(pid), 0.0)      # Different for B vs W
    sigma_parent_pb = get_parent_sigma_pb(int(pid))   # Different for B vs W

    mask_parent = np.abs(parent_id) == pid
    eff_parent = np.sum(weights[mask_parent] * P_decay[mask_parent]) / w_sum

    total_expected += lumi_fb * (sigma_parent_pb * 1e3) * BR_parent * eff_parent
```

Each parent species (B⁰, B±, W±, Z) is weighted by its OWN cross-section and branching ratio. No double-counting!

## Example at m = 5 GeV

Suppose at 5 GeV:
- Pythia generates 10k HNLs from B mesons
- MadGraph generates 5k HNLs from W bosons

**Combined CSV has 15k rows:**
```csv
event,weight,hnl_id,parent_pdg,pt,eta,phi,...
1,1.0,9900014,511,12.3,1.5,0.8,...     # From B⁰
2,1.0,9900014,521,8.7,-0.5,2.1,...     # From B±
...
10001,1.0,9900014,24,45.2,2.1,-1.3,... # From W⁺
10002,1.0,9900014,-24,38.5,-1.8,0.6,...# From W⁻
...
```

**Analysis correctly sums:**
```
N_sig(|U|²) = N_B(|U|²) + N_W(|U|²)
            = [σ_B × BR_B(|U|²) × ε_B] + [σ_W × BR_W(|U|²) × ε_W]
            ✓ Correct!
```

## References

- Review: `PHYSICS_REVIEW.md` Section 2.2
- Tool: `analysis_pbc/limits/combine_production_channels.py`
- Analysis: `analysis_pbc/limits/expected_signal.py`

---

**Status:** Documented and resolved for future production runs
**Date:** 2025-12-02
