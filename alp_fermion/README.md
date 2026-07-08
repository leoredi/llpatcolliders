# `alp_fermion/` — PBC benchmark BC10 (fermiophilic ALP) for GRENDEL

A pseudoscalar axion-like particle `a` with **universal coupling to SM
fermions** (single parameter `1/f`, Wilson coefficients `c_f = 1`), produced in
**B → K a** (the flavour-violating `b → s a` top penguin) and decaying to the
heaviest kinematically open charged fermions (`a → μμ, ττ, cc/ss …`). This is
the genuine charged-track "ALP" benchmark.

BC10 is **coupling-controlled**: the single `1/f` fixes the production yield, the
lifetime cτ, and the visible branching ratios *simultaneously*. The limit is
therefore set **model-completely** (like the HNL `|U|²`, not like the
model-agnostic `(BR, cτ)` scan of `higgs/`): for each `(m_a, 1/f)` the actual
cτ, production, and visible BR are evaluated and `N_signal ≥ 3` (background-free,
3000 fb⁻¹) is required, producing a **closed island** in the `(m_a, 1/f)` plane.

## What is reused (imported, never copied)
- `higgs/grendel_geometry.py`, `higgs/reco_common.py` — the shared GRENDEL mesh
  and the 4-hit two-track vertex reconstruction (`wall_inner_outer`,
  `reconstruct_3d`, `timing_chi2_4hit`) + PR#13 selection.
- `hnl/analysis/decay_reco_acceptance.py` — `build_event_mc`, `best_two_directions`,
  `selection_mask`, `scan_u2` (the acceptance MC + lifetime reweighting).
- `hnl/analysis/_engine.py`, `format_bridge.py`, `exclusion.py` — geometry
  ray-cast, CSV→kinematics bridge, island extraction.
- `hnl/production/fonll` + `decay_engine/kinematics.py` — FONLL bottom sampler
  (14 TeV) and the two-body decay sampler for B → K a.

The key trick that makes the reuse exact: map the coupling to the HNL scan
variable `u2 ≡ (1/f / 1/f_ref)²`. Off the single reference point
(`model.INV_F_REF`), production ∝ `u2` (since BR(B→K a) ∝ (1/f)²) and
cτ ∝ `1/u2` (since Γ_tot ∝ (1/f)²) — exactly the structure `scan_u2` already
implements. The only new analysis code is the `u2 ↔ 1/f` remap and the island
extraction.

## Layout
```
alp_fermion/
  model.py          model layer: a→ff widths, Γ_tot, cτ, BRs, BR(B→K a) (b→s a penguin)
  alp_production.py B (FONLL) → B→K a → a four-vector CSVs (weight,E,px,py,pz), weighted by BR(B→K a)
  templates.py      per-mass ALP rest-frame decay templates (channels ∝ BR), FairShip-compatible npz
  sensitivity.py    (m_a, 1/f) closed-island scan; N_signal ≥ 3, 3000 fb⁻¹
  plot.py           BC10 island in the (m_a, 1/f) plane (+ optional overlay curves)
  paths.py          output-dir policy (ALP_TMP_DIR override) + hnl/ import shim
  tests/test_model.py
  EXTERNAL_INPUTS_NEEDED.md   (B→K a normalization, data-driven hadronic width, overlay curves)
```

## Run (use the `hnl` conda env — has trimesh + networkx)
```
PY=/Volumes/sandbox/conda/envs/hnl/bin/python
$PY -m alp_fermion.alp_production --n-pool 120000      # a four-vector CSVs
$PY -m alp_fermion.templates       --n-templates 20000 # decay templates
$PY -m alp_fermion.sensitivity                         # island CSV + plot
$PY -m pytest alp_fermion/tests -q
```
Outputs go to `alp_fermion/tmp/` (gitignored): `analysis/bc10_sensitivity.csv`
and `analysis/bc10_island.{png,pdf}`.

## Model layer references
- a→f f̄ widths and L_int: Bauer–Neubert–Thamm, JHEP 12 (2017) 044
  (arXiv:1708.00443); EFT review arXiv:2012.12272.
- b→s a top penguin (B→K a): BNT and "Flavour probes of ALPs"
  (arXiv:2110.10698); ALPINIST (arXiv:2105.10806).
- Benchmark definition / curves: 2025 PBC report (arXiv:2505.00947), using the
  unified FIP calculation (arXiv:2311.00507).
See `EXTERNAL_INPUTS_NEEDED.md` for the two analytic placeholders (B→K a
absolute normalization; data-driven hadronic width above 2 m_π) and the overlay
curves.
