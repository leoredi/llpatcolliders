# Session Notes: Tau Jaggedness Investigation (Jan 2026)

## Summary
The tau money plot jaggedness was traced to **inconsistent inputs** in the tau pipeline rather than a physics bug. Two main causes:

1) **Stale combined files missing EW contributions**  
   Some tau `_combined.csv` files were created *before* EW production finished. Those combined files therefore lacked EW rows even though the EW files existed later.

2) **Mixed production statistics for tau fromTau**  
   Older tau masses (e.g., 0.25/0.30/0.40/0.45 GeV) still used fromTau files produced with **100k–1M** events, while newer runs used **10M**. This caused point-to-point statistical fluctuations in tau limits.

3) **Decay file selection was not tau-aware**  
   Tau decay files use `lightfonly` / `lightfstau` / `lightfstauK` naming above 0.42 GeV, but the selection logic only recognized electron/muon categories. This was fixed so tau uses the correct decay-file tags.

## Key Observations
- **EW (MadGraph) always has ~100k rows** by construction (one HNL per event).
- **Pythia meson rows vary** widely and are much lower for tau, especially fromTau.
- Tau thus has larger statistical noise; e/μ are smoother because meson rows are O(1e5–2e5).
- Combined files are written to `output/csv/simulation/` and, prior to the patch, could delete originals.

## Changes Made
- `analysis_pbc/limits/combine_production_channels.py` now **keeps originals by default**.
  - New flag `--delete-originals` performs deletion explicitly.
- `analysis_pbc/decay/rhn_decay_library.py` now uses **flavour-aware decay file selection** and recognizes tau tags (`lightfonly`, `lightfstau`, `lightfstauK`).
- `config_mass_grid.py`: electron/muon meson grids now include the tau fine low-mass points (unified high granularity).
  - Added low-mass points for e/μ: 0.23, 0.26, 0.29, 0.32, 0.38, 0.41, 0.44, 0.47, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.35, 1.45, 1.55.
  - New e/μ meson grid size: 93 points (was 75). Regenerate production for these new masses.
- EW grids are now unified to the same 93-point grid for all flavours. Regenerate MadGraph for the added low-mass points if you want full EW coverage.

## Current Tau Mass Grid
From `config_mass_grid.py`:

```
0.20, 0.23, 0.26, 0.29, 0.32, 0.35, 0.38, 0.41, 0.44, 0.47, 0.50,
0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00,
1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40, 1.45, 1.50, 1.55, 1.60,
1.70, 1.80, 1.90, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00, 3.20, 3.40
```

Masses 0.25/0.30/0.40/0.45 are **not** in the current tau grid and are remnants from older runs.

## Commands Used (Clean Pipeline)

### 1) Re-generate tau Pythia (10M fromTau)
Run in bash from `production/pythia_production/`.

### 2) Re-combine tau channels
```
cd analysis_pbc
conda run -n llpatcolliders python limits/combine_production_channels.py --flavour tau
```

### 3) Re-run tau analysis
```
cd analysis_pbc
conda run -n llpatcolliders python run_tau_only.py
```

### 4) Re-plot money plot
```
cd money_plot
conda run -n llpatcolliders python plot_money_island.py
```

## Current Status (end of session)
- Tau combined files recreated with recent timestamps.
- Tau limits regenerated with `run_tau_only.py` (89 points).
- Money plot updated after re-analysis.

## Remaining Notes
- If residual wiggles remain in tau below ~1 GeV, they are likely MC statistics and discrete decay-file mass points; fix by raising fromTau stats further, adjusting decay seeds, or apply plotting-only smoothing.
- Additional small jaggedness sources:
  - Decay file mass snapping (nearest file selection) above the low-mass analytical threshold.
  - |U|² scan quantization (100 steps between 1e-12 and 1e-2; log step ≈ 0.101).
  - RNG seed dependence in decay sampling (different `--decay-seed` can shift limits by ~1 grid step).
- EW/meson row imbalance is expected and not a bug; normalization happens via σ × BR in analysis.



 Current Status: ✓ All 44 tau combined files have consistent timestamps (Jan 25 09:35-09:36), no stale files detected.

  Row Count by Mass Regime
  ┌────────────┬───────────────┬──────┬───────────┬───────────┐
  │   Regime   │  Mass Range   │  EW  │   Meson   │   Total   │
  ├────────────┼───────────────┼──────┼───────────┼───────────┤
  │ Kaon       │ 0.20-0.47 GeV │ 100k │ ~44k      │ ~144k     │
  ├────────────┼───────────────┼──────┼───────────┼───────────┤
  │ Charm      │ 0.50-1.60 GeV │ 100k │ ~110-132k │ ~210-232k │
  ├────────────┼───────────────┼──────┼───────────┼───────────┤
  │ Transition │ 1.70-1.90 GeV │ 100k │ 87        │ ~100k     │
  ├────────────┼───────────────┼──────┼───────────┼───────────┤
  │ Beauty     │ 2.00+ GeV     │ 100k │ ~86k      │ ~186k     │
  └────────────┴───────────────┴──────┴───────────┴───────────┘
  e/μ/τ Meson Comparison
  ┌──────────┬─────────────────┬──────────────────┐
  │ Flavour  │ Kaon (0.35 GeV) │ Charm (1.00 GeV) │
  ├──────────┼─────────────────┼──────────────────┤
  │ Electron │ 46,598          │ 199,163          │
  ├──────────┼─────────────────┼──────────────────┤
  │ Muon     │ 46,238          │ 198,869          │
  ├──────────┼─────────────────┼──────────────────┤
  │ Tau      │ 44,572          │ 110,031          │
  └──────────┴─────────────────┴──────────────────┘
  Key Observations

  1. Transition region (1.70-1.90 GeV): Only 87 charm rows — this is expected physics. The D_s→τν→HNL phase space closes as m_HNL approaches m_τ (1.777 GeV). These
   masses are dominated by EW production.
  2. Tau has ~2× fewer charm rows than e/μ at 1.00 GeV (110k vs 199k). This is because tau production from D_s→τν→HNL has tighter kinematics than direct D→ℓν→HNL.
  3. Residual jaggedness below ~1 GeV in the plot likely reflects the lower tau meson statistics (~44k kaon, ~110k charm) compared to e/μ (~46k kaon, ~199k charm).
   This is intrinsic to the tau production mechanism, not a bug.
  4. No stale files: Old mass points (0.25, 0.30, 0.40, 0.45 GeV) confirmed absent.




  
