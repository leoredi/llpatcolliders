## MadGraph Form-Factor Semileptonic Samples (Optional Third Production Step)

Goal: replace the phase-space Pythia semileptonic decays (B→DℓN, D→KℓN, Λb→ΛcℓN, …) with MadGraph samples that include realistic form factors. These are **shape-only** samples; normalization still comes from `production_xsecs.py` (σ) and `production_brs()` (BR) in the analysis.

### What to generate
- Channels: charm (`D→KℓN`) and beauty (`B→DℓN`, `Bs→DsℓN`, `Λb→ΛcℓN`, `Λc→ΛℓN` as available).
- Flavours: electron, muon, tau (one coupling set to 1 at generation time).
- Mass grid: match `config_mass_grid.py` masses for the corresponding flavour.

### File naming (required)
Place the converted CSVs in `output/csv/simulation/` with:
```
HNL_<mass>GeV_<flavour>_charm_ff.csv
HNL_<mass>GeV_<flavour>_beauty_ff.csv
```
Use 2-decimal mass strings (e.g., `2p60`), same columns as Pythia/MG EW (`event,weight,hnl_id,parent_pdg,pt,eta,phi,p,E,mass,prod_x_mm,prod_y_mm,prod_z_mm,beta_gamma`), and set `parent_pdg` to the true parent (411/421/431 for charm; 511/521/531/5122 for beauty baryon; 4122 for Λc).

### Normalization rules
- `weight` should be a relative MC weight (typically 1.0). **Do not** bake in cross-sections.
- Cross-sections and BRs are applied later by the analysis; these samples only provide kinematics.

### Workflow (after standard Pythia meson + MG EW steps)
1) Generate MadGraph semileptonic LHE with form factors for the desired channels/masses/flavours.
2) Convert LHE → CSV (match the column schema and naming above).
3) Drop the CSVs into `output/csv/simulation/`. You may keep or remove the Pythia charm/beauty CSVs; the pipeline now **prefers** `*_ff` files and will ignore the base versions if both exist.
4) Run:
   ```
   cd analysis_pbc
   conda run -n llpatcolliders python limits/combine_production_channels.py
   conda run -n llpatcolliders python limits/run.py --parallel
   ```

### Why this works
- `combine_production_channels.py` and `run.py` now match `_ff` regimes and automatically pick them over the phase-space versions, so no double counting.
- Per-parent normalization still uses `production_xsecs.get_parent_sigma_pb` and `production_brs()`, keeping physics consistent.

### TODO (not automated)
A dedicated MG driver and cards for these decay-level processes are not included here. Until then, treat this as an expert/manual step for power users who want improved semileptonic kinematics.
