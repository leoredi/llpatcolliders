# Mass Ranges for HNL Production Regimes

## Mass Grid Configuration

**Uniform Grid:** 0.1 GeV spacing from 0.2 to 10.0 GeV
- **Total points:** 99 per flavor
- **Range:** 0.2, 0.3, 0.4, ..., 9.8, 9.9, 10.0 GeV
- **Applies to:** All three flavors (electron, muon, tau)

## Production Regime Selection

**Source:** Lines 159-163 in `main_hnl_production.cc`

```cpp
std::string getProductionRegime(double mHNL) {
    if (mHNL < 0.5) return "kaon";       // Kaon-dominated regime
    if (mHNL < 2.0) return "charm";      // Charm-dominated regime
    return "beauty";                      // Beauty regime (2.0-10.0 GeV)
}
```

### Production Regimes

| Regime | Mass Range (GeV) | Production Channel | Pythia Card | Parent Mesons |
|--------|------------------|-------------------|-------------|---------------|
| **Kaon** | **m < 0.5** | K± → ℓ± N | `hnl_Kaon.cmnd` | K±, K⁰_L |
| **Charm** | **0.5 ≤ m < 2.0** | D mesons → ℓ± N | `hnl_Dmeson.cmnd` | D⁰, D±, D_s±, Λ_c |
| **Beauty** | **2.0 ≤ m ≤ 10.0** | B mesons → ℓ± N | `hnl_Bmeson.cmnd` | B⁰, B±, B_s, Λ_b |

**Note:** Masses above ~5 GeV will use the beauty card but may produce few/no HNLs due to kinematic constraints (B meson mass ~5.28 GeV).

## From Documentation (CLAUDE.md)

The documentation states:
- K-dominance: 0.2-0.5 GeV
- D-dominance: 0.5-2 GeV
- B-dominance: 2-5 GeV

### Comparison

✅ **Documentation matches code for supported regimes**

The documentation ranges are accurate for the meson production regimes (< 5 GeV).

## Production Cross-Sections by Regime

From `config/production_xsecs.py` and CLAUDE.md:

| Regime | Typical σ (pb) | Relative Importance |
|--------|---------------|-------------------|
| Kaon (m < 0.5 GeV) | σ(pp → K⁺) ≈ 5.0 × 10¹⁰ pb | Highest rate, lowest BR(K→ℓN) |
| Charm (0.5-2 GeV) | σ(pp → D⁰) ≈ 2.8 × 10¹⁰ pb | High rate, moderate BR |
| Beauty (2-5 GeV) | σ(pp → B⁰) ≈ 4.0 × 10⁸ pb | Lower rate, higher BR |

## Production Statistics

**Uniform 0.1 GeV grid from 0.2 to 10.0 GeV:**
- **Electron**: 99 mass points (all regimes)
- **Muon**: 99 mass points (all regimes)
- **Tau**: 99 mass points (dual mode for m < 1.64 GeV)

**Total simulations:** ~312 runs (99 + 99 + 99 + ~15 tau fromTau mode)

**Regime distribution:**
- **Kaon (0.2-0.4 GeV)**: 3 points per flavor
- **Charm (0.5-1.9 GeV)**: 15 points per flavor
- **Beauty (2.0-10.0 GeV)**: 81 points per flavor

**Note:** High-mass points (5-10 GeV) may produce few/no HNLs due to kinematic limits.

## Physics Rationale for Boundaries

### m = 0.5 GeV (Kaon → Charm transition)
- **Below 0.5 GeV**: Kaon decays (K → ℓN) dominate
  - m_K± = 0.494 GeV, m_K⁰ = 0.498 GeV
  - Phase space opens fully around m_HNL ~ 0.2-0.4 GeV
- **Above 0.5 GeV**: D meson decays start to dominate
  - m_D⁰ = 1.865 GeV, m_D± = 1.870 GeV
  - Larger BR(D→ℓN) compensates for lower σ(D production)

### m = 2.0 GeV (Charm → Beauty transition)
- **Below 2.0 GeV**: D meson decays dominate
- **Above 2.0 GeV**: B meson decays become important
  - m_B⁰ = 5.280 GeV, m_B± = 5.279 GeV
  - Excellent phase space for 2-5 GeV HNLs

### m = 5.0 GeV (Upper limit of meson production)
- **At 5.0 GeV**: B mesons (m_B ≈ 5.28 GeV) still kinematically allowed
- **Above 5.0 GeV**: B meson decays to HNL become increasingly suppressed
- **Production limit**: This system supports up to 5.0 GeV

## File Naming Convention

**Pattern:** `HNL_{mass}GeV_{flavor}_{regime}.csv`

**Mass format:** Decimal point → 'p' (e.g., 2.6 GeV → `2p60`)

**Examples:**
```
HNL_0p40GeV_electron_kaon.csv       # 0.4 GeV, kaon regime
HNL_1p20GeV_muon_charm.csv          # 1.2 GeV, charm regime
HNL_2p60GeV_electron_beauty.csv     # 2.6 GeV, beauty regime
HNL_4p80GeV_muon_beauty.csv         # 4.8 GeV, beauty regime
```

## Special Cases

### Tau Coupling
- **Cannot use kaon regime**: m_τ = 1.777 GeV >> m_K
- **Starts at charm regime**: Minimum mass ~0.5 GeV for tau
- **Two production modes:**
  - `direct`: B/D → τ N (mixing at vertex)
  - `fromTau`: B/D → τ ν, then τ → N X (mixing at tau decay)

## Summary

**Mass grid:** Uniform 0.1 GeV spacing from 0.2 to 10.0 GeV (99 points per flavor)

**Production regimes (from getProductionRegime):**
- Kaon: m < 0.5 GeV
- Charm: 0.5 ≤ m < 2.0 GeV
- Beauty: 2.0 ≤ m ≤ 10.0 GeV

**Kinematically viable range:** 0.2 - 5.0 GeV (limited by B meson mass ~5.28 GeV)

**Status:** ✅ Uniform grid configured for all three flavors
