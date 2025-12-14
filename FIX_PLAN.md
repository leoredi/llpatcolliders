# Remaining Improvements / TODO (Actionable)

This file tracks **only what is left to do** in `llpatcolliders`. Completed items are intentionally removed to keep this list short and current.

---

## Remaining TODOs (short list)

1. **Geometry performance + correctness**: vectorize ray-tracing and handle multi-intersection rays in `analysis_pbc/geometry/per_parent_efficiency.py`.
2. **Cross-section documentation clarity**: add explicit literature refs + assumptions in `analysis_pbc/config/production_xsecs.py`.
3. **Optional security hardening**: replace `eval()` in `analysis_pbc/models/hnl_model_hnlcalc.py` with an AST whitelist evaluator.

---

## 1) Geometry performance + multi-intersection correctness

**Problem:** `analysis_pbc/geometry/per_parent_efficiency.py` currently:
- loops row-by-row (`df.iterrows()` around `analysis_pbc/geometry/per_parent_efficiency.py:323`) → slow for O(1e5–1e6) rows,
- calls `mesh.ray.intersects_location` one-ray-at-a-time (`analysis_pbc/geometry/per_parent_efficiency.py:345`),
- assumes exactly two intersections and uses only the first pair for `path_length` (`analysis_pbc/geometry/per_parent_efficiency.py:355`), which can be wrong for complex geometry.

### Concrete implementation steps

File: `analysis_pbc/geometry/per_parent_efficiency.py`

1. **Vectorize direction calculation**
   - Build `theta = 2 * arctan(exp(-eta))` for all rows (already implemented in `eta_phi_to_direction` but currently called per-row).
   - Compute `ray_directions` as an `(N, 3)` float array.
   - Build a `valid_mask` for finite `eta/phi` and finite, non-zero direction norm; drop invalid rows once up-front.

2. **Batch ray intersections**
   - Use `mesh.ray.intersects_location(ray_origins=origins, ray_directions=directions)` with batches (e.g. 10k rays per batch).
   - Use `origins = np.repeat(origin_arr[None, :], N_batch, axis=0)`.

3. **Compute distances along ray robustly**
   - For each returned intersection point `loc`, compute parametric distance `t = dot(loc - origin, direction_unit)` (not Euclidean `||loc-origin||`, which can be wrong if direction isn’t perfectly normalized).

4. **Handle multi-intersection rays**
   - Intersections come with `index_ray`; group `t` values per ray.
   - Sort `t` per ray and compute:
     - `entry_distance = t_sorted[0]`
     - `path_length = Σ_k (t_sorted[2k+1] - t_sorted[2k])` over all full pairs
   - Mark `hits_tube = True` only if a ray has at least one full pair (`len(t_sorted) >= 2`).

5. **Chunk-level robustness**
   - If trimesh raises `RTreeError` for a batch, split the batch and retry (binary split), then fall back to per-ray only for the smallest failing chunk.
   - Keep counters for skipped rays (invalid dirs, RTreeError) and print a final summary.

6. **Preserve the output contract**
   - Keep output columns: `beta_gamma`, `hits_tube`, `entry_distance`, `path_length`.
   - Preserve the current “old/new CSV compatibility” mapping (`parent_pdg→parent_id`, `p→momentum`).

### Acceptance / validation

- Time a representative run (e.g. ~100k rows) before/after and record wall time.
- Compare `hits_tube` fraction and `entry_distance/path_length` distributions before/after:
  - allow small numerical differences,
  - but large shifts indicate a bug in `t` computation or multi-intersection handling.

---

## 2) Cross-section references + assumptions (docs/comments only)

**Problem:** the numbers in `analysis_pbc/config/production_xsecs.py` are “PBC-ish”, but the file does not clearly document:
- which measurements/theory inputs the constants correspond to,
- whether `σ(W)`/`σ(Z)` are inclusive and at what perturbative order,
- why charm fragmentation fractions don’t sum exactly to 1.

### Concrete changes

File: `analysis_pbc/config/production_xsecs.py`

1. Add explicit references (comment-only is fine) for:
   - `SIGMA_CCBAR_PB`, `SIGMA_BBBAR_PB`
   - `SIGMA_W_PB`, `SIGMA_Z_PB`
   - `SIGMA_KAON_PB` (and warn that it’s the dominant systematic for kaon regime)
2. Clarify `σ(W)` / `σ(Z)` definitions (inclusive W±? Z/γ* mass window?).
3. Add a one-line comment explaining why `FRAG_C_*` sum ≠ 1 (missing excited states/baryons) and why that’s acceptable at this level.

---

## 3) Optional hardening: remove `eval()` from HNLCalc wrapper

**Problem:** `analysis_pbc/models/hnl_model_hnlcalc.py` uses `eval()` (even with restricted builtins) to evaluate HNLCalc-provided BR expressions. This is a security/code-review smell.

### Concrete changes

File: `analysis_pbc/models/hnl_model_hnlcalc.py`

1. Replace `_safe_eval()` with an AST-based whitelist evaluator:
   - Allowed nodes: literals, names (`hnl`, `mass`, `coupling`, `np`), arithmetic ops, tuples/lists, and calls/attributes only on `hnl`/`np`.
   - Reject everything else: comprehensions, subscripting, imports, lambdas, unknown attribute access, etc.
2. Add a tiny unit-style check (no framework needed) that ensures a malicious expression is rejected.

---

## Verification checklist (for remaining items)

After **geometry refactor**:
- `cd analysis_pbc && conda run -n llpatcolliders python tests/closure_anubis/test_expected_signal_events_kernel.py`
- Run a representative geometry preprocessing and compare hit fractions/distributions vs a known baseline.

After **any physics/kernel change**:
- `cd analysis_pbc && conda run -n llpatcolliders python tests/test_26gev_muon.py`

