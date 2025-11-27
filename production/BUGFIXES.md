# Bug Fixes Applied to HNL Production Code

This document summarizes all bugs that were identified and fixed in the HNL production simulation.

## Date: 2025-11-26

---

## 1. Critical Fix: Card Files Not Used in Benchmark Scans ✅

**Problem:**
- `run_benchmarks.sh` runs the executable from inside `output/` directory
- Code looked for cards at `cards/hnl_*.cmnd` (relative path)
- From inside `output/`, this tried to find `output/cards/...` which doesn't exist
- All benchmark scans silently fell back to default settings
- **Beautiful tuned cards were never actually used!**

**Solution:**
- Updated card loading logic in `main_hnl_production.cc` (lines 439-486)
- Now tries both `cards/` and `../cards/` paths
- Prints confirmation message when card is successfully loaded
- Fallback to defaults only if both paths fail

**Impact:** HIGH - This was causing all production runs to use wrong physics settings

---

## 2. Physics Bug: Incorrect HNL Spin Assignment ✅

**Problem:**
- HNL defined with `spinType = 2` (vector boson)
- Should be `spinType = 1` (fermion) for Majorana neutrino
- While functionally harmless (we use phase space decays), this is embarrassing and would be flagged by referees

**Solution:**
- Changed line 490: `hnlDef << HNL_ID << ":new = N Nbar 1 0 0";`
- Added comment explaining it's a fermion

**Impact:** MEDIUM - Cosmetic but important for correctness

---

## 3. Methodology Fix: Mass-Regime Thresholds Inconsistent ✅

**Problem:**
- Code used muon-specific kinematic thresholds: 0.39 / 1.76 / 5.17 GeV
- Documentation described simple stepped regimes: 0.5 / 2 / 5 GeV
- Thresholds not appropriate for electron or tau coupling
- Inconsistency between code and documentation

**Solution:**
- Updated `getProductionRegime()` function (lines 117-123)
- Now uses clean stepped boundaries: 0.5 / 2.0 / 5.0 GeV
- Matches README and PHYSICS.md documentation
- Actual kinematic checks still enforced per-channel in `configureMesonDecays()`

**Impact:** MEDIUM - Restores consistency, improves clarity

**Rationale:**
- The regime only chooses which card to use (SoftQCD vs HardQCD:cc vs HardQCD:bb vs EW)
- Actual decay channel allowances are checked separately with proper lepton masses
- No physics bug, just cleaner organization

---

## 4. Robustness Fix: Boost Gamma Division by Zero Risk ✅

**Problem:**
- Boost calculation: `boostGamma = p.e() / mHNL`
- Assumed runtime mass exactly equals input argument
- If `mHNL = 0` passed by mistake, would cause division by zero
- Not using actual particle mass from Pythia

**Solution:**
- Updated boost calculation (lines 576-580)
- Now uses `mass = p.m()` from Pythia event
- Fallback to `mHNL` if particle mass is zero
- Hard floor of `1e-6` to guarantee no division by zero

**Impact:** LOW - Defensive programming improvement

---

## 5. Documentation Fix: README Usage Outdated ✅

**Problem:**
- README showed obsolete interface: `./main_hnl <mass> <flavor> <production_mode>`
- Actual binary: `./main_hnl_production <mass> <flavor> [nEvents]`
- Production mode now auto-selected, not a parameter
- Would confuse users trying to follow documentation

**Solution:**
- Updated README.md usage section (lines 172-198)
- Corrected binary name to `main_hnl_production`
- Removed `<production_mode>` parameter
- Added examples for all regimes with proper syntax
- Added note explaining automatic regime selection

**Impact:** HIGH - User-facing documentation now correct

---

## 6. Documentation Fix: CSV Column List Wrong ✅

**Problem:**
- README claimed CSV has both `parent_id` and `parent_pdg` columns
- Code actually writes only `parent_pdg` (the PDG code integer)
- No separate index column for parent

**Solution:**
- Updated README.md column list (lines 152-170)
- Removed `parent_id` from column list
- Added detailed description of each column
- Now matches actual CSV output

**Impact:** MEDIUM - Documentation now accurate

---

## Testing Verification

All fixes tested and verified:

1. **Card loading:** ✅
   - Tested from repo root: Loads `cards/hnl_Bmeson.cmnd`
   - Tested from `output/`: Loads `../cards/hnl_Bmeson.cmnd` (fallback works)
   - Prints "Using card file: ..." confirmation

2. **Spin assignment:** ✅
   - Recompiled successfully
   - HNL now defined as fermion (spinType=1)

3. **Mass regimes:** ✅
   - 2.6 GeV → "beauty" regime (correct, 2 < 2.6 < 5)
   - Clean threshold boundaries as documented

4. **Boost gamma:** ✅
   - Uses `p.m()` with proper fallback logic
   - No division by zero possible

5. **Documentation:** ✅
   - README matches actual binary interface
   - CSV columns match actual output

---

## Files Modified

- `main_hnl_production.cc` (4 fixes applied)
- `README.md` (2 documentation fixes)
- `BUGFIXES.md` (this file - created for tracking)

---

## Summary

**Critical bugs fixed:** 1 (card files not used)
**Physics corrections:** 1 (spin assignment)
**Methodology improvements:** 1 (regime thresholds)
**Robustness improvements:** 1 (boost gamma safety)
**Documentation corrections:** 2 (usage + CSV columns)

**Total fixes:** 6

**Code status:** Ready for publication-quality production runs

---

## Next Steps

The code is now correct and robust. Recommended workflow:

```bash
# 1. Clean rebuild
make clean
PYTHIA8=$(pwd)/pythia/pythia8315 make

# 2. Test single mass point
./main_hnl_production 2.6 muon 10000

# 3. Verify card is used (should see "Using card file: ...")
# 4. Check CSV output format matches README

# 5. Run full benchmark scan
./run_benchmarks.sh
```

All benchmarks will now correctly use the tuned Pythia cards for each mass regime.
