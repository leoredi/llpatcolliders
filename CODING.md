# Code Architecture Documentation

## 1. Language and environment

- **Analysis pipeline:** Python 3.11 (conda env `llpatcolliders`)
- **Production simulation:** C++17 linked against Pythia 8.315; MadGraph5 via Docker
- **Shell orchestration:** Bash (not fish, despite macOS default)
- **Key Python dependencies:** numpy, pandas, matplotlib, scipy, trimesh, shapely, rtree, tqdm

Install: `conda env create -f environment.yml && conda activate llpatcolliders`

Build C++: `cd production/pythia_production && PYTHIA8=$(pwd)/pythia8315 make`

## 2. File-by-file reference

### 2.1 `config_mass_grid.py`

Central mass grid definition. 146 points from 0.2 to 17.0 GeV.

- **MASS_GRID**: sorted list of floats. Fine spacing (15 MeV) below 2 GeV, coarser above.
- **format_mass_for_filename(mass)**: `2.60 → "2p60"`. Used everywhere for CSV naming.
- Imported by: `run_parallel_production.sh` (via `load_mass_grid.sh`), `run_hnl_scan.py`, `run.py`

### 2.2 `production/pythia_production/main_hnl_production.cc`

C++ event generator. Registers HNL as particle 9900012, configures meson decay channels with 100% BR, runs Pythia event loop, writes CSV.

**Key functions:**
- `getLeptonInfo(flavour)` → PDG IDs, masses, labels
- `getProductionRegime(mass, flavour)` → "kaon" / "charm" / "beauty"
- `configureMesonDecays(pythia, mass, ...)` → enables 2-body and 3-body channels per regime
- `configureMesonDecaysToTauNu(pythia, mass)` → enables parent → τν for fromTau mode
- `configureTauDecays(pythia, mass)` → configures τ → HNL X with adaptive branching ratios
- `findPhysicalParent(event, i)` → traces event record to find true parent PDG
- `main()` → argument parsing, Pythia setup, event loop, CSV write

**CLI:** `./main_hnl_production <mass_GeV> <flavour> [nEvents] [fromTau]`

**Output CSV columns:** event, weight, hnl_id, parent_pdg, tau_parent_id, pt, eta, phi, p, E, mass, prod_x_mm, prod_y_mm, prod_z_mm, beta_gamma

### 2.3 `production/pythia_production/run_parallel_production.sh`

Orchestrates parallel Pythia runs across all mass points and flavours. Sources `load_mass_grid.sh` for the mass array.

**Key variables:**
- `MAX_PARALLEL=12` — concurrent job limit
- `NEVENTS=100000` — events per mass point
- `FROMTAU_MASS_THRESHOLD=1.77` — max mass for fromTau mode
- `FLAVOUR` / `MODE` — from CLI args

**Key functions:**
- `count_jobs()` — counts background jobs
- `wait_for_slot()` — blocks until a slot opens
- `run_production_job(mass, flavour, mode)` — launches one Pythia run in background, moves output CSV to `output/csv/simulation/`
- `count_tau_runs()` — accounts for dual-mode tau runs

**Output structure:** `output/csv/simulation/HNL_<mass>GeV_<flavour>_<regime>[_fromTau].csv`

### 2.4 `production/pythia_production/load_mass_grid.sh`

Loads MASS_GRID from Python into a bash array. Auto-detects conda environment. Called via `source ./load_mass_grid.sh`.

### 2.5 `production/pythia_production/monitor_production.sh`

Status dashboard: shows running job count, CSV file count, recent completions, latest log tail.

### 2.6 `production/madgraph_production/scripts/run_hnl_scan.py`

MadGraph5 driver for EW production. Runs inside Docker.

**Pipeline per (flavour, mass):**
1. `generate_process()` — creates MG5 process directory from proc_card template
2. `write_cards_to_process()` — writes run_card.dat and param_card.dat with mass/mixing substitution
3. `generate_events()` — runs `bin/generate_events` to produce LHE
4. `extract_cross_section()` — parses banner/log for σ
5. `convert_lhe_to_csv()` — calls `LHEParser` to produce pipeline-compatible CSV
6. `cleanup_workdir()` — removes MG5 working directory

**Configuration:** mixing set to |U|² = 1 at generation time. Physical coupling applied analytically in analysis.

### 2.7 `production/madgraph_production/scripts/lhe_to_csv.py`

Parses LHE files and extracts HNL (PDG 9900012) 4-vectors. Determines parent W/Z via:
1. Direct mother chain traversal
2. Process ID mapping from LHE header `<init>` block
3. Default parent fallback if only one boson type in file

Output matches Pythia CSV format exactly.

### 2.8 `production/madgraph_production/scripts/validate_xsec.py`

Sanity-checks MadGraph cross-sections against literature values (arXiv:1805.08567). Log-space interpolation of expected ranges per mass.

### 2.9 `analysis_pbc/config/production_xsecs.py`

Stores all production cross-sections, fragmentation fractions, and tau branching ratios.

**Key functions:**
- `get_parent_sigma_pb(pdg_id)` → cross-section in pb for a parent particle, including fragmentation
- `get_parent_tau_br(pdg_id)` → BR(parent → τν) for tau-chain production
- `get_sigma_summary()` → formatted table of all cross-sections

**Particle coverage:** K±(321), K_L(130), D⁰(421), D±(411), D_s±(431), Λ_c(4122), B⁰(511), B±(521), B_s(531), Λ_b(5122), W(24), Z(23)

### 2.10 `analysis_pbc/models/hnl_model_hnlcalc.py`

Safe wrapper around the vendored HNLCalc.

**Class `HNLModel(mass_GeV, Ue2, Umu2, Utau2)`:**
- `ctau0_m` property → proper decay length in metres
- `production_brs()` → dict {abs(parent_pdg): branching_ratio} for all kinematically allowed channels
  - Iterates 2-body and 3-body channels from HNLCalc
  - Adds W and Z contributions with explicit phase-space and helicity factors:
    - BR(W → ℓN) = |U|² × BR(W → ℓν)_SM × (1 − r²)² × (1 + r²) where r = m_N/m_W
    - BR(Z → νN) = |U|² × BR(Z → νν̄)_SM/3 × same structure

**Safe expression evaluator:** HNLCalc channel definitions contain string expressions like `"hnl.get_2body_br(411, 13)"`. These are evaluated via a restricted AST walker (`_SafeExprEvaluator`) that only allows:
- Names: `hnl`, `mass`, `coupling`, `np`
- HNLCalc methods: `get_2body_br`, `get_3body_dbr_*`, `get_2body_br_tau`
- NumPy: `sqrt`, `pi`
- Basic arithmetic: +, -, *, /, **

### 2.11 `analysis_pbc/geometry/per_parent_efficiency.py`

Builds detector mesh and computes ray-mesh intersections.

**`build_drainage_gallery_mesh()`:**
- 46 hard-coded vertices in x-z plane at y = 22 m (global coordinates)
- Tube radius 1.54 m (1.4 m × 1.1 envelope)
- Returns trimesh.Trimesh object (32-segment rings)

**`preprocess_hnl_csv(csv_path, output_path)`:**
- Reads production CSV
- Converts (eta, phi) to 3D unit directions
- Batch ray-mesh intersection (batch_size=10000)
- Handles trimesh RTreeError via binary subdivision of batches
- Computes entry_distance, path_length per event from intersection hit data
- Saves augmented DataFrame with hits_tube, entry_distance, path_length columns

**`eta_phi_to_directions(eta, phi)`:** vectorised conversion via θ = 2·arctan(e⁻η)

### 2.12 `analysis_pbc/decay/decay_detector.py`

Computes decay acceptance from pre-computed decay event libraries.

**`build_decay_cache(geom_df, mass, flavour, ...)`:**
- Loads decay events from RHN library (or generates with MG5+Pythia)
- For each HNL that hits the tube:
  - Picks a random decay event
  - Boosts all decay products from rest frame to lab frame (Lorentz boost along HNL direction)
  - Stores charged product directions as unit vectors
- Returns `DecayCache` with pre-computed charged directions and random decay uniforms

**`compute_decay_acceptance(geom_df, ctau0_m, mesh, decay_cache, separation_m)`:**
- Decay probability: P = exp(-d_entry/λ) × [1 − exp(-d_path/λ)] where λ = βγ·cτ₀
- For each event with P > 0:
  - Projects charged product directions from a decay point (sampled from exponential) to the mesh
  - Computes pairwise track separations at the mesh surface
  - Requires ≥ 2 charged tracks separated by ≥ separation_m

**`compute_separation_pass_static(geom_df, mesh, decay_cache, separation_m)`:**
- Faster version using fixed decay point at midpoint of path through tube
- Returns boolean array (one per event)

**Lorentz boost implementation:**
- `_rotation_matrix_from_z(direction)` — Rodrigues formula rotation
- `_boost_along_direction(p4, beta_gamma, direction)` — full relativistic boost

### 2.13 `analysis_pbc/decay/rhn_decay_library.py`

Manages external MATHUSLA RHN decay event databases.

**File selection logic in `select_decay_file(mass, flavour)`:**
1. Lists all `.txt` files in the flavour's decay directory
2. Categorises by name pattern (inclDs, inclDD, inclD, nocharm, nocharmnoss, lightfonly, analytical)
3. Below low_mass_threshold (0.42–0.53 GeV): prefers analytical
4. Above threshold: follows priority order (inclDs > inclDD > ...)
5. Selects closest mass within each priority category
6. Warns if |Δm| > 0.5 GeV

**`load_decay_events(path)`:** parses text files where events are separated by blank lines, each line is `E, px, py, pz, mass, pid`.

### 2.14 `analysis_pbc/decay/generate_hnl_decay_events.py`

Generates HNL decay samples using MadGraph5 + Pythia8 hadronisation.

**Process:** n1 → ℓ q q̄, n1 → ν q q̄, n1 → ℓ ℓ ν, etc.

**Output CSV columns:** event, pid, E, px, py, pz, mass, mass_GeV, ctau0_m, gamma_tot_GeV

### 2.15 `analysis_pbc/limits/expected_signal.py`

Core signal calculation.

**`couplings_from_eps2(eps2, benchmark)`:** maps eps2 + benchmark code ("100"/"010"/"001") to (Ue2, Umu2, Utau2).

**`expected_signal_events(geom_df, mass_GeV, eps2, benchmark, lumi_fb, ...)`:**
1. Gets HNL model (ctau₀, production BRs) — either fresh from HNLCalc or via scaling
2. Filters events that hit the tube
3. For each parent species:
   - Looks up σ_parent from production_xsecs
   - Gets BR(parent → HNL) from model (or scaled)
   - Computes per-event decay probability P_i = exp(-d_entry/λ) × [1 − exp(-d_path/λ)]
   - Applies separation cut
   - Sums: N += L × σ × BR × Σ(w_i × P_i × sep_i)
4. For tau-chain events (tau_parent_id > 0):
   - Uses parent tau BR × parent → τν BR
5. Returns total N_sig

**`scan_eps2_for_mass(geom_df, mass, benchmark, lumi_fb, ...)`:**
- Scans 100 log-spaced eps2 values from 10⁻¹² to 10⁻²
- Computes HNLCalc once at eps2_ref, scales for all others
- Finds eps2_min, eps2_max where N_sig crosses 2.996
- Returns dict with eps2_min, eps2_max, peak_events

### 2.16 `analysis_pbc/limits/run.py`

Main orchestrator. Manages file discovery, geometry caching, parallelisation, and output aggregation.

**`run_flavour(flavour, benchmark, lumi_fb, ...)`:**
1. Discovers production CSVs in `output/csv/simulation/`
2. Groups by mass point
3. For each mass: preprocesses geometry (cached), builds decay cache, runs eps2 scan
4. Collects results into DataFrame
5. Saves to `output/csv/analysis/HNL_U2_limits_summary.csv`

**File selection logic:**
- Prefers `*_combined.csv` (output of combine_production_channels.py)
- Falls back to best individual regime: beauty > charm > kaon
- EW files (`*_ew.csv`) added as supplementary

**Parallelisation:** ProcessPoolExecutor with configurable worker count.

**Geometry caching:** saves preprocessed CSVs to `output/csv/geometry/`. Invalidates if source CSV is newer.

**CLI entry point:** `python limits/run.py --parallel --workers 12 [--flavour X] [--mass Y]`

### 2.17 `analysis_pbc/limits/combine_production_channels.py`

Merges production CSVs at each mass point across regimes.

**File pattern:** `HNL_<mass>GeV_<flavour>_<regime>[_ff][_direct|_fromTau].csv`

**Selection rules:**
- Form-factor files (`_ff`) preferred over phase-space when both exist
- One file per (regime, mode) pair
- Output: `HNL_<mass>GeV_<flavour>_combined.csv` with added columns: source_regime, source_mode, source_is_ff

### 2.18 `analysis_pbc/limits/timing_utils.py`

Single context manager `_time_block(name, timing_dict)` for performance profiling. No-op if timing_dict is None.

### 2.19 `analysis_pbc/run_tau_only.py`

Convenience script: reruns tau analysis and merges into existing results CSV without recomputing electron/muon.

### 2.20 `analysis_pbc/scripts/check_hnlcalc_scaling.py`

Validation script: verifies ctau₀ ∝ 1/|U|² and BR ∝ |U|² across multiple mass/flavour/eps2 combinations. Reports relative errors against tolerance (default 5×10⁻⁴).

### 2.21 `analysis_pbc/tests/test_scaling_vs_per_eps2.py`

Unit test: creates a mock geometry DataFrame and verifies that the scaled computation matches the full per-eps2 HNLCalc computation to < 10⁻¹⁰ relative error.

### 2.22 `money_plot/plot_money_island.py`

Reads `HNL_U2_limits_summary.csv` and produces 3-panel exclusion plot (electron, muon, tau).

**`append_tip_point_if_needed()`:** closes the island at the high-mass edge by interpolating a synthetic point where peak_events crosses the 3σ threshold. Prevents the island from ending with an open contour.

## 3. Data flow and caching

```
[Production CSVs]
    → preprocess_hnl_csv() → [Geometry cache CSVs]
    → build_decay_cache() → [in-memory DecayCache]
    → expected_signal_events() → N_sig(mass, eps2)
    → scan_eps2_for_mass() → (eps2_min, eps2_max)
    → run_flavour() → [results CSV]
    → plot_money_island() → [PNG]
```

**Geometry cache:** `output/csv/geometry/geom_HNL_<mass>GeV_<flavour>_<regime>.csv`. Invalidated when source is newer. Delete the directory to force full recomputation.

**Decay cache:** in-memory only (not persisted). Built once per mass point and reused across the eps2 scan.

## 4. Naming conventions

- Masses in filenames: `2p60` (two decimals, dot → p)
- Flavours: `electron`, `muon`, `tau` (lowercase)
- Regimes: `kaon`, `charm`, `beauty`, `ew`
- Modes: `direct`, `fromTau`
- Benchmarks: `"100"` (electron), `"010"` (muon), `"001"` (tau)
- Parent PDGs use absolute values in BR dictionaries (e.g. 521 for B±)

## 5. Error handling patterns

- Missing production files: skipped with warning, mass point gets NaN limits
- Missing cross-sections: warning printed, parent contribution zeroed
- Missing tau BRs: warning, tau-chain contribution zeroed
- HNLCalc failures: caught and logged, mass point gets NaN
- trimesh RTreeError: batch subdivided and retried
- Empty CSVs (< 1000 bytes): filtered out by combine_production_channels.py
