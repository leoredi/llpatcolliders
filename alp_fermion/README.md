# `alp_fermion/` — PBC benchmark BC10 (fermiophilic ALP) for GRENDEL

A pseudoscalar axion-like particle `a` with **universal coupling to SM
fermions** (single parameter `1/f`, Wilson coefficients `c_f = 1`), produced in
**B → K⁽ⁱ⁾ a** (the flavour-violating `b → s a` penguin, summed over the full
kaon tower as in GKOZ arXiv:2310.03524) and decaying to leptons and — via the
data-driven GKOZ spectral tables — to hadrons. This is the genuine
charged-track "ALP" benchmark.

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
  model.py          model layer: a→ff/γγ/hadrons widths (GKOZ tables), Γ_tot, cτ, BRs,
                    BR(B→K⁽ⁱ⁾ a) over the kaon tower (b→s a penguin, one-loop RG coefficient)
  alp_production.py B (FONLL) → B→K⁽ⁱ⁾ a → a four-vector CSVs (weight,E,px,py,pz)
  templates.py      per-mass ALP rest-frame decay templates (channels ∝ visible BR)
  sensitivity.py    (m_a, 1/f) closed-island scan; N_signal ≥ 3, 3000 fb⁻¹
  plot.py           BC10 island in the (m_a, 1/f) plane (+ optional overlay curves)
  paths.py          output-dir policy (ALP_TMP_DIR override) + hnl/ import shim
  data/alpinist/    digitized GKOZ decay-width tables via ALPINIST (PROVENANCE.md, pinned SHA)
  data/senscalc_2501/ provenance and production audit for the pinned arXiv:2501.04525 upgrade
  data/published/   canonical published sensitivity curve (CSV + MANIFEST; see its README)
  tools/compute_cbs_alpinist.py  regenerates the b→s a RG coefficient (BSD-3 ALPINIST port)
  tools/export_senscalc_2501.{py,wls}  verifies and exports the three SensCalc decay inputs
  tests/test_model.py
  EXTERNAL_INPUTS_NEEDED.md   (input provenance + residual systematics; overlays dropped)
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
- Coupling convention + a→f f̄ widths: Bauer–Neubert–Thamm, JHEP 12 (2017) 044
  (arXiv:1708.00443). Our 1/f axis is the BNT `g_aff = c_f m_f/f` convention
  (= 2/f in the GKOZ normalisation).
- Production + decay phenomenology: GKOZ, "ALPs with universal fermion
  couplings — revisited" (arXiv:2310.03524), via the ALPINIST implementation
  (arXiv:2105.10806, BSD-3, pinned SHA in `data/alpinist/COMMIT_SHA.txt`):
  one-loop RG b→s a coefficient (`tools/compute_cbs_alpinist.py`), kaon-tower
  form factors (Boiarska et al., arXiv:1904.10447), digitized hadronic + γγ
  width tables (`data/alpinist/PROVENANCE.md` — incl. the normalisation
  cross-checks).
- Benchmark definition: 2025 PBC report (arXiv:2505.00947); unified FIP
  calculation arXiv:2311.00507; state-of-the-art hadronic treatment
  arXiv:2501.04525.
See `EXTERNAL_INPUTS_NEEDED.md` for residual systematics (production
normalisation <~20% in 1/f; >3 GeV width from the matched perturbative
continuation) and the dropped overlay curves.
