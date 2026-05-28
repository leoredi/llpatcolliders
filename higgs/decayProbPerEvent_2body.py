import os
import sys
import numpy as np
import pandas as pd
import matplotlib
# Batch mode: use the non-interactive Agg backend unless the user has set
# MPLBACKEND or passes --interactive. Default is batch so the script runs
# on headless machines without a DISPLAY.
INTERACTIVE = ('--interactive' in sys.argv)
if not INTERACTIVE and not os.environ.get('MPLBACKEND'):
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.integrate import quad


def show_or_close():
    """plt.show() in interactive mode, plt.close('all') in batch mode."""
    if INTERACTIVE:
        plt.show()
    else:
        plt.close('all')

from grendel_geometry import (
    SPEED_OF_LIGHT, DETECTOR_THICKNESS,
    calculate_decay_length, cache_geometry, mesh_fiducial,
)

M_ELECTRON = 0.000511  # GeV/c²

# Analysis cuts
P_CUT   = 0.600    # GeV/c — minimum electron momentum
SEP_MIN = 0.01    # m — minimum separation at detector (1 cm)
SEP_MAX = 10.0      # m — maximum separation at detector
DCA_CUT = 0.1  # m (10 cm) — maximum DCA between reconstructed tracks
# Conditional max-sep_outer cut, motivated by the geometric statement that
# a real vertex satisfies sep_outer = open_angle × (D + L) with D bounded
# by the fiducial path length. For near-parallel tracks (open_angle at the
# angular-resolution noise floor), the largest sep_outer a real vertex can
# produce is ~θ_noise × max_path ≈ 30 mrad × 50 m ≈ 1.5 m. Above that,
# we are looking at two physically separated parallel particles
# (e.g., beam-halo or cosmic coincidences), not a vertex.
# When open_angle is clearly above the noise floor (diverging tracks),
# no upper bound on sep_outer is applied.
THETA_PARALLEL      = 0.100  # rad (100 mrad — covers parallel + moderate-opening regime)
SEP_OUT_MAX_PARALLEL = 0.30  # m — max sep_outer applied only when open_angle < THETA_PARALLEL

# IP-muon-transit veto. A single muon shooting outward from the CMS IP makes
# four hits that lie on a single straight line. The "collinearity" is the RMS
# perpendicular residual of all 4 hits to their best-fit line (SVD of the 4×3
# centered hit matrix; sqrt((S[1]^2 + S[2]^2) / 4)). For a real two-body decay
# the 4 hits trace a "V", so collinearity ≈ (sep_in + sep_out)/4 — large. For
# a muon transit it sits at the hit-noise floor (~σ_hit). To protect collimated
# signal (where the V is narrow and naturally looks like a line), the cut only
# fires when sep_outer > SEP_OUT_COLLIN_GATE, i.e. the topology is not in the
# collimated regime where signal and muon are geometrically indistinguishable.
COLLIN_MIN          = 0.030  # m (30 mm) — min collinearity for events above the gate
SEP_OUT_COLLIN_GATE = 0.30   # m — collinearity cut applies only when sep_outer > this

outString = "0p5GeV"
sample_csv = "LLP0p5GeVSmall.csv"

# Tracking resolution
HIT_RESOLUTION = 0.003  # m (3 mm per layer)
N_LAYERS       = 2       # number of tracking layers (stations)

# ============================================================
# Two-body decay acceptance (analytical)
# ============================================================
#
# LLP (mass M, momentum p_LLP) → e+ e-
#
# Rest frame: E* = M/2,  p* = sqrt(M²/4 - m_e²)
# Isotropic in cosθ*, φ*
#
# Lab frame (boost along LLP direction, m_e → 0 limit):
#   E_{1,2}  = γ M/2 (1 ± β cosθ*)
#   |p|_{1,2} ≈ E_{1,2}  (since E >> m_e)
#
# Opening angle:
#   cos(θ_12) = 1 - 2 / [γ²(1 - β² cos²θ*)]
#
# Accepted region of |cosθ*|:
#
#   c_S(d)      < |cosθ*| < min(c_P, c_max_sep(d), β)
#   ↑ min sep                 ↑ momentum cut  ↑ max sep  ↑ forward
#
# Accepted fraction: A(d) = max(0, c_upper(d) - c_S(d))
#
# Detectable decay probability:
#   P = ∫_{d_entry}^{d_exit} (1/λ) e^{-d/λ} · A(d) dd
#

def compute_c_upper(gamma, beta, mass, p_cut=P_CUT):
    """
    Upper limit on |cosθ*| from momentum cut and forward requirement.
    
    Momentum cut on softer electron:
      |p_2| > p_cut  →  E_2 > sqrt(p_cut² + m_e²)
      γ M/2 (1 - β|cosθ*|) > sqrt(p_cut² + m_e²)
      |cosθ*| < (1 - 2·sqrt(p_cut² + m_e²) / (γM)) / β
    
    Forward requirement: both electrons forward → |cosθ*| < β
    """
    E_min = np.sqrt(p_cut**2 + M_ELECTRON**2)
    c_P = (1.0 - 2.0 * E_min / (gamma * mass)) / beta
    c_P = min(c_P, 1.0)
    return min(c_P, beta)


def compute_c_S(d_remaining, gamma, beta, sep_min=SEP_MIN):
    """
    Minimum |cosθ*| from minimum separation requirement.
    
    Separation ≈ θ_12 × d_remaining > sep_min
    Need θ_12 > θ_min → need |cosθ*| > c_S
    
    From cos(θ_12) = 1 - 2/(γ²(1 - β²cos²θ*)):
      θ_12 > θ_min  ↔  cos(θ_12) < cos(θ_min)
      β²cos²θ* > 1 - 2/(γ²(1 - cos(θ_min)))
    """
    if d_remaining <= 0:
        return 1.0
    
    theta_min = sep_min / d_remaining
    if theta_min <= 0:
        return 0.0
    
    cos_theta_min = np.cos(min(theta_min, np.pi))
    denom = gamma**2 * (1.0 - cos_theta_min)
    if denom <= 0:
        return 0.0
    
    val = (1.0 - 2.0 / denom) / beta**2
    if val <= 0:
        return 0.0
    return np.sqrt(val)


def compute_c_max_sep(d_remaining, gamma, beta, sep_max=SEP_MAX):
    """
    Maximum |cosθ*| from maximum separation requirement.
    
    Separation ≈ θ_12 × d_remaining < sep_max
    Need θ_12 < θ_max → need |cosθ*| < c_max_sep
    
    From cos(θ_12) = 1 - 2/(γ²(1 - β²cos²θ*)):
      θ_12 < θ_max  ↔  cos(θ_12) > cos(θ_max)
      β²cos²θ* < 1 - 2/(γ²(1 - cos(θ_max)))
    """
    if d_remaining <= 0:
        return 1.0  # No constraint at exit point
    
    theta_max = sep_max / d_remaining
    if theta_max >= np.pi:
        return 1.0  # No constraint
    
    cos_theta_max = np.cos(theta_max)
    denom = gamma**2 * (1.0 - cos_theta_max)
    if denom <= 0:
        return 1.0
    
    val = (1.0 - 2.0 / denom) / beta**2
    if val <= 0:
        return 0.0  # Even cosθ*=0 exceeds max separation
    if val >= 1.0:
        return 1.0  # No constraint
    return np.sqrt(val)


def acceptance_weighted_decay_prob(entry_d, exit_d, gamma, beta, mass,
                                    decay_length, p_cut=P_CUT,
                                    sep_min=SEP_MIN, sep_max=SEP_MAX):
    """
    Compute decay probability weighted by two-body decay acceptance.
    
    P = ∫_{entry}^{exit} (1/λ) e^{-d/λ} · A(d) dd
    
    A(d) = max(0, c_upper(d) - c_lower(d))
    where:
      c_lower = c_S(d)                          [min separation]
      c_upper = min(c_P, c_max_sep(d), β)       [momentum + max sep + forward]
    """
    c_P_upper = compute_c_upper(gamma, beta, mass, p_cut)
    
    if c_P_upper <= 0:
        return 0.0
    
    def integrand(d):
        d_remaining = exit_d - d
        c_lower = compute_c_S(d_remaining, gamma, beta, sep_min)
        c_upper_sep = compute_c_max_sep(d_remaining, gamma, beta, sep_max)
        c_upper = min(c_P_upper, c_upper_sep)
        A = max(0.0, c_upper - c_lower)
        return (1.0 / decay_length) * np.exp(-d / decay_length) * A
    
    result, _ = quad(integrand, entry_d, exit_d, limit=100,
                     epsabs=1e-15, epsrel=1e-10)
    return result


def unweighted_decay_prob(entry_d, exit_d, decay_length):
    """Decay probability with no acceptance cuts."""
    path_length = exit_d - entry_d
    return np.exp(-entry_d / decay_length) * (1.0 - np.exp(-path_length / decay_length))


# ============================================================
# Processing functions
# ============================================================

def process_with_acceptance(csv_file, lifetime_seconds, geo_cache,
                             p_cut=P_CUT, sep_min=SEP_MIN, sep_max=SEP_MAX):
    df = pd.read_csv(csv_file)
    n = len(df)
    hits = geo_cache['hits']
    
    df['hits_tube'] = hits
    df['decay_probability'] = 0.0
    df['decay_probability_no_cuts'] = 0.0
    df['acceptance'] = 0.0
    
    for idx in np.where(hits)[0]:
        p = geo_cache['momentum'][idx]
        m = geo_cache['mass'][idx]
        gamma = geo_cache['gamma'][idx]
        beta = geo_cache['beta'][idx]
        entry = geo_cache['entry_d'][idx]
        exit_ = geo_cache['exit_d'][idx]
        decay_length = calculate_decay_length(p, m, lifetime_seconds)
        
        p_with = acceptance_weighted_decay_prob(
            entry, exit_, gamma, beta, m, decay_length, p_cut, sep_min, sep_max)
        p_no = unweighted_decay_prob(entry, exit_, decay_length)
        
        df.iat[idx, df.columns.get_loc('decay_probability')] = p_with
        df.iat[idx, df.columns.get_loc('decay_probability_no_cuts')] = p_no
        if p_no > 0:
            df.iat[idx, df.columns.get_loc('acceptance')] = p_with / p_no
    
    event_stats = []
    for event_idx, group in df.groupby('event'):
        n_part = len(group)
        p1 = group.iloc[0]['decay_probability'] if n_part >= 1 else 0
        p2 = group.iloc[1]['decay_probability'] if n_part >= 2 else 0
        p1_nc = group.iloc[0]['decay_probability_no_cuts'] if n_part >= 1 else 0
        p2_nc = group.iloc[1]['decay_probability_no_cuts'] if n_part >= 2 else 0
        
        p_at_least_one = 1 - (1 - p1) * (1 - p2)
        p_at_least_one_nc = 1 - (1 - p1_nc) * (1 - p2_nc)
        
        event_stats.append({
            'event': event_idx,
            'n_particles_hitting_tube': group['hits_tube'].sum(),
            'prob_at_least_one_decays': p_at_least_one,
            'prob_at_least_one_no_cuts': p_at_least_one_nc,
            'prob_both_decay': p1 * p2,
            'particle1_decay_prob': p1,
            'particle2_decay_prob': p2,
        })
    
    event_df = pd.DataFrame(event_stats)
    return df, event_df


def analyze_decay_vs_lifetime(csv_file, geo_cache, lifetime_range,
                               p_cut=P_CUT, sep_min=SEP_MIN, sep_max=SEP_MAX):
    df_base = pd.read_csv(csv_file)
    n_events = df_base['event'].nunique()
    
    results = {
        'lifetimes': lifetime_range,
        'mean_single_particle_decay_prob': [],
        'mean_single_no_cuts': [],
        'mean_at_least_one_decay_prob': [],
        'mean_at_least_one_no_cuts': [],
        'mean_both_decay_prob': [],
        'exclusion': [],
        'exclusion_no_cuts': [],
        'mean_acceptance': [],
        'total_events': n_events
    }
    
    for lifetime in tqdm(lifetime_range, desc="Scanning lifetimes"):
        df, event_df = process_with_acceptance(
            csv_file, lifetime, geo_cache, p_cut, sep_min, sep_max)
        
        hits = df[df['hits_tube']]
        mean_single = hits['decay_probability'].mean() if len(hits) > 0 else 0
        mean_single_nc = hits['decay_probability_no_cuts'].mean() if len(hits) > 0 else 0
        mean_acc = hits['acceptance'].mean() if len(hits) > 0 else 0
        
        mean_p1 = event_df['prob_at_least_one_decays'].mean()
        mean_p1_nc = event_df['prob_at_least_one_no_cuts'].mean()
        mean_both = event_df['prob_both_decay'].mean()
        
        results['mean_single_particle_decay_prob'].append(mean_single)
        results['mean_single_no_cuts'].append(mean_single_nc)
        results['mean_at_least_one_decay_prob'].append(mean_p1)
        results['mean_at_least_one_no_cuts'].append(mean_p1_nc)
        results['mean_both_decay_prob'].append(mean_both)
        results['exclusion'].append(3 / (mean_p1 * 3000 * 52E3))
        results['exclusion_no_cuts'].append(3 / (mean_p1_nc * 3000 * 52E3))
        results['mean_acceptance'].append(mean_acc)
    
    return results


def sample_separations(geo_cache, lifetime_seconds, n_samples_per_particle=100,
                       rng_seed=42, hit_resolution=HIT_RESOLUTION,
                       n_layers=N_LAYERS):
    """
    Monte Carlo sample decay positions and rest-frame angles to build
    distributions of electron-pair separations, pointing angles, and DCA.

    For each particle that hits the fiducial volume:
      1. Sample decay position d from (1/λ) exp(-d/λ) within [entry, exit]
      2. Sample |cosθ*| uniformly in [0, 1]
      3. Compute lab-frame electron angles from boost kinematics
      4. Place true hits at inner/outer tracking layers (separated by
         DETECTOR_THICKNESS), then smear by hit_resolution
      5. Reconstruct tracks from smeared hits → separation, pointing, DCA
      6. Check that reconstructed vertex (PCA) is inside the fiducial volume

    Set hit_resolution=0 to recover truth-level distributions (DCA=0).

    Returns:
        separations, weights, momenta, pointing_angles, p_soft, dca,
        vtx_in_fiducial, open_angles, sep_outers, d_implieds, collinearities
        (all 1-D arrays of the same length)
    """
    rng = np.random.default_rng(rng_seed)

    hits = geo_cache['hits']
    hit_idx = np.where(hits)[0]

    all_seps = []
    all_weights = []
    all_momenta = []
    all_pointing = []
    all_p_soft = []
    all_dca = []
    all_vtx_in = []
    all_open = []
    all_sep_out = []
    all_dimp = []
    all_collin = []

    for idx in hit_idx:
        entry = geo_cache['entry_d'][idx]
        exit_ = geo_cache['exit_d'][idx]
        gamma = geo_cache['gamma'][idx]
        beta = geo_cache['beta'][idx]
        mass = geo_cache['mass'][idx]
        p_llp = geo_cache['momentum'][idx]
        N = n_samples_per_particle

        decay_length = calculate_decay_length(p_llp, mass, lifetime_seconds)
        path_length = exit_ - entry

        # Inverse CDF sampling of decay position within [entry, exit]
        u = rng.uniform(0, 1, N)
        exp_entry = np.exp(-entry / decay_length)
        exp_exit = np.exp(-exit_ / decay_length)
        denom = exp_entry - exp_exit
        if denom < 1e-300:
            continue

        d_samples = -decay_length * np.log(exp_entry - u * denom)
        D = exit_ - d_samples  # distance from vertex to inner tracking layer

        p_decay = exp_entry * (1 - np.exp(-path_length / decay_length))
        w = p_decay / N

        # Sample rest-frame decay angle
        cos_theta_star = rng.uniform(0, 1, N)
        sin_theta_star = np.sqrt(1 - cos_theta_star**2)

        # Softer electron momentum (truth, not affected by resolution)
        p_soft = gamma * mass / 2 * (1 - beta * cos_theta_star)

        # Lab-frame angles from LLP direction (small-angle regime)
        # e1 (harder): tanα₁ = sinθ* / γ(cosθ*+β)
        # e2 (softer): tanα₂ = sinθ* / γ(β-cosθ*)  — opposite side of LLP axis
        tan_a1 = sin_theta_star / (gamma * (cos_theta_star + beta))
        tan_a2 = sin_theta_star / (gamma * (beta - cos_theta_star))

        L = DETECTOR_THICKNESS

        # ---- 3D reconstruction in local frame ----
        # Origin = true vertex, ẑ = LLP direction, x̂ = decay plane
        # True hit positions at inner (z=D) and outer (z=D+L) tracking layers
        x1_in_true  = D * tan_a1
        x2_in_true  = -D * tan_a2
        x1_out_true = (D + L) * tan_a1
        x2_out_true = -(D + L) * tan_a2
        # y = 0 for both (in the decay plane)

        if hit_resolution > 0:
            noise = rng.normal(0, hit_resolution, (8, N))
            x1_in  = x1_in_true  + noise[0]
            y1_in  = noise[1]
            x1_out = x1_out_true + noise[2]
            y1_out = noise[3]
            x2_in  = x2_in_true  + noise[4]
            y2_in  = noise[5]
            x2_out = x2_out_true + noise[6]
            y2_out = noise[7]
        else:
            x1_in  = x1_in_true;  y1_in  = np.zeros(N)
            x1_out = x1_out_true; y1_out = np.zeros(N)
            x2_in  = x2_in_true;  y2_in  = np.zeros(N)
            x2_out = x2_out_true; y2_out = np.zeros(N)

        # Reconstructed track directions (unnormalised)
        dx1 = x1_out - x1_in;  dy1 = y1_out - y1_in   # dz = L
        dx2 = x2_out - x2_in;  dy2 = y2_out - y2_in

        # --- Separation at inner and outer tracking layers ---
        sep     = np.sqrt((x1_in  - x2_in )**2 + (y1_in  - y2_in )**2)
        sep_out = np.sqrt((x1_out - x2_out)**2 + (y1_out - y2_out)**2)

        # --- Pointing angle (bisector vs ẑ = LLP direction) ---
        mag1 = np.sqrt(dx1**2 + dy1**2 + L**2)
        mag2 = np.sqrt(dx2**2 + dy2**2 + L**2)
        bx = dx1 / mag1 + dx2 / mag2
        by = dy1 / mag1 + dy2 / mag2
        bz = L / mag1 + L / mag2
        bmag = np.sqrt(bx**2 + by**2 + bz**2)
        pointing = np.arccos(np.clip(bz / bmag, -1, 1))

        # --- Opening angle between the two reconstructed tracks ---
        cos_open = (dx1 * dx2 + dy1 * dy2 + L * L) / (mag1 * mag2)
        open_angle = np.arccos(np.clip(cos_open, -1, 1))

        # --- Implied vertex distance: d_implied = sep_inner / open_angle ---
        # Geometric identity for two tracks emerging from a common point:
        # sep_inner = open_angle × (vertex → inner-layer distance).
        # Robust against the PCA-degeneracy that hits near-parallel tracks.
        d_implied = np.where(open_angle > 1e-9, sep / open_angle, np.inf)

        # --- Collinearity: RMS perpendicular residual of the 4 hits to their
        # best-fit line. ≈ σ_hit for a single-line topology (IP-muon transit),
        # ≈ (sep_inner + sep_outer)/4 for a true V-shaped decay vertex.
        pts = np.stack([
            np.stack([x1_in,  y1_in,  np.zeros(N)],   axis=-1),
            np.stack([x1_out, y1_out, np.full(N, L)], axis=-1),
            np.stack([x2_in,  y2_in,  np.zeros(N)],   axis=-1),
            np.stack([x2_out, y2_out, np.full(N, L)], axis=-1),
        ], axis=1)  # (N, 4, 3)
        ctr = pts.mean(axis=1, keepdims=True)
        _, S_svd, _ = np.linalg.svd(pts - ctr, full_matrices=False)
        collin = np.sqrt((S_svd[:, 1]**2 + S_svd[:, 2]**2) / 4)

        # --- DCA between reconstructed tracks ---
        # Lines: P_i + t_i * d_i, with P_i at inner-layer hits (z = D)
        # Cross product n = d1 × d2
        nx = dy1 * L - L * dy2
        ny = L * dx2 - dx1 * L
        nz = dx1 * dy2 - dy1 * dx2
        n_mag = np.sqrt(nx**2 + ny**2 + nz**2)

        wx = x1_in - x2_in
        wy = y1_in - y2_in
        # wz = 0 (both at z = D)

        dca = np.where(n_mag > 1e-30,
                       np.abs(wx * nx + wy * ny) / n_mag,
                       np.sqrt(wx**2 + wy**2))

        # --- PCA (point of closest approach) → reconstructed vertex ---
        a_coeff = dx1**2 + dy1**2 + L**2
        b_coeff = dx1 * dx2 + dy1 * dy2 + L**2
        c_coeff = dx2**2 + dy2**2 + L**2
        d_coeff = dx1 * wx + dy1 * wy
        e_coeff = dx2 * wx + dy2 * wy
        det = a_coeff * c_coeff - b_coeff**2
        det_safe = np.where(np.abs(det) < 1e-30, 1e-30, det)
        t1 = (b_coeff * e_coeff - c_coeff * d_coeff) / det_safe
        t2 = (a_coeff * e_coeff - b_coeff * d_coeff) / det_safe

        # PCA midpoint z in local frame (vertex at z=0)
        pca_z = D + 0.5 * (t1 + t2) * L

        # Reconstructed vertex distance from IP along LLP direction
        d_vtx = d_samples + pca_z
        vtx_in = (d_vtx > entry) & (d_vtx < exit_)

        all_seps.append(sep)
        all_weights.append(np.full(N, w))
        all_momenta.append(np.full(N, p_llp))
        all_pointing.append(pointing)
        all_p_soft.append(p_soft)
        all_dca.append(dca)
        all_vtx_in.append(vtx_in)
        all_open.append(open_angle)
        all_sep_out.append(sep_out)
        all_dimp.append(d_implied)
        all_collin.append(collin)

    if not all_seps:
        empty = np.array([])
        return (empty, empty, empty, empty, empty, empty,
                np.array([], dtype=bool), empty, empty, empty, empty)

    return (np.concatenate(all_seps),
            np.concatenate(all_weights),
            np.concatenate(all_momenta),
            np.concatenate(all_pointing),
            np.concatenate(all_p_soft),
            np.concatenate(all_dca),
            np.concatenate(all_vtx_in),
            np.concatenate(all_open),
            np.concatenate(all_sep_out),
            np.concatenate(all_dimp),
            np.concatenate(all_collin))




def build_cutflow(seps, pointing, weights, momenta, p_soft, dca, vtx_in,
                  open_angle, sep_outer, collinearity,
                  p_cut=P_CUT, sep_min=SEP_MIN, sep_max=SEP_MAX,
                  dca_cut=DCA_CUT, theta_parallel=THETA_PARALLEL,
                  sep_out_max_parallel=SEP_OUT_MAX_PARALLEL,
                  collin_min=COLLIN_MIN,
                  sep_out_collin_gate=SEP_OUT_COLLIN_GATE,
                  pointing_cut=None):
    """
    Apply signal selection cuts sequentially and return a cutflow table.

    Parameters
    ----------
    dca_cut : float
        Maximum DCA between reconstructed tracks (m).
    theta_parallel, sep_out_max_parallel : float
        Conditional max-sep_outer cut: when open_angle < theta_parallel,
        require sep_outer < sep_out_max_parallel. For larger open_angle
        no upper bound on sep_outer is applied.
    collin_min, sep_out_collin_gate : float
        IP-muon-transit veto. When sep_outer > sep_out_collin_gate, require
        collinearity > collin_min (rejects events where all 4 hits sit on
        a single straight line, i.e. an IP-shot muon traversing the tunnel).
        The gate protects collimated signal (small sep_outer) from a cut
        that cannot discriminate it from muons at this resolution.
    pointing_cut : float or None
        Maximum pointing angle (rad). If None, no pointing cut is applied.

    Returns
    -------
    cutflow : list of dicts with keys 'cut', 'efficiency',
              'marginal_efficiency', 'weighted_yield'
    """
    total_w = weights.sum()
    if total_w == 0:
        return []

    mask = np.ones(len(seps), dtype=bool)
    rows = []

    def add_row(name, m):
        w = weights[m].sum()
        prev_w = weights[prev_mask].sum()
        rows.append({
            'cut': name,
            'weighted_yield': w,
            'efficiency': w / total_w,
            'marginal_efficiency': w / prev_w if prev_w > 0 else 0,
        })

    prev_mask = mask.copy()
    add_row('Decay in fiducial', mask)

    prev_mask = mask.copy()
    mask = mask & (p_soft >= p_cut)
    add_row(f'p_soft > {p_cut*1000:.0f} MeV/c', mask)

    prev_mask = mask.copy()
    mask = mask & (seps >= sep_min) & (sep_outer >= sep_min)
    add_row(f'sep_in & sep_out > {sep_min*1000:.0f} mm', mask)

    prev_mask = mask.copy()
    mask = mask & (seps <= sep_max)
    add_row(f'sep_in < {sep_max*100:.0f} cm', mask)

    prev_mask = mask.copy()
    mask = mask & (dca <= dca_cut)
    add_row(f'DCA < {dca_cut*100:.1f} cm', mask)

    # Conditional max sep_outer: tight bound only when tracks are
    # nearly parallel (open_angle below the angular-resolution floor)
    prev_mask = mask.copy()
    is_parallel = open_angle < theta_parallel
    parallel_ok = (~is_parallel) | (sep_outer < sep_out_max_parallel)
    mask = mask & parallel_ok
    add_row(f'sep_out<{sep_out_max_parallel:.1f}m if θ<{theta_parallel*1000:.0f}mrad',
            mask)

    # Conditional collinearity cut: kills IP-muon transits (4 hits on a
    # single line). Applied only when sep_outer is above the gate so that
    # collimated signal is unaffected.
    prev_mask = mask.copy()
    gated = sep_outer > sep_out_collin_gate
    collin_ok = (~gated) | (collinearity > collin_min)
    mask = mask & collin_ok
    add_row(f'collin>{collin_min*1000:.0f}mm if sep_out>{sep_out_collin_gate*100:.0f}cm',
            mask)

    prev_mask = mask.copy()
    mask = mask & vtx_in
    add_row('Vertex in fiducial (PCA)', mask)

    if pointing_cut is not None:
        prev_mask = mask.copy()
        mask = mask & (pointing <= pointing_cut)
        add_row(f'pointing < {pointing_cut*1000:.1f} mrad', mask)

    return rows


# ============================================================
# Main
# ============================================================
if __name__ == "__main__":
    origin = [0, 0, 0]
    
    print(f"  Detector thickness: {DETECTOR_THICKNESS*100:.0f} cm")
    print(f"  Cuts: p_e > {P_CUT*1000:.0f} MeV/c, "
          f"{SEP_MIN*1000:.0f} mm < separation < {SEP_MAX*100:.0f} cm")
    
    geo_cache = cache_geometry(sample_csv, mesh_fiducial, origin)
    
    # --- Diagnostics ---
    print("\n" + "="*50)
    print("ACCEPTANCE DIAGNOSTICS")
    print("="*50)
    
    test_mass = 15.0
    print(f"\nM = {test_mass} GeV → e+e-, p_cut = {P_CUT*1000:.0f} MeV/c, "
          f"sep: {SEP_MIN*1000:.0f} mm – {SEP_MAX*100:.0f} cm")
    print(f"\n{'p_LLP':>8} {'γ':>6} {'c_P':>7} {'θ_min':>10} "
          f"{'d(1mm)':>8} {'d(10cm)':>9} {'θ_12=10cm/1m':>14}")
    
    for p_test in [20, 50, 100, 200, 500]:
        E = np.sqrt(p_test**2 + test_mass**2)
        g = E / test_mass
        b = p_test / E
        c_P = compute_c_upper(g, b, test_mass)
        theta_min = 2.0 / g
        d_min_sep = SEP_MIN / theta_min
        d_max_sep = SEP_MAX / theta_min
        sep_at_1m = theta_min * 1.0
        print(f"{p_test:>8.0f} {g:>6.1f} {c_P:>7.4f} {theta_min*1000:>8.1f} mrad "
              f"{d_min_sep*100:>7.2f}cm {d_max_sep:>7.1f}m {sep_at_1m*100:>12.1f}cm")
    
    print(f"\nMax sep cut effect: at d_remaining = 1m, θ_max = {SEP_MAX}m / 1m = 100 mrad")
    print("  → kills decays far from detector for low-γ (large opening angle)")
    print(f"  → for γ=1.67 (p=20 GeV): θ_min=1200 mrad >> 100 mrad → heavily constrained")
    print(f"  → for γ=33.4 (p=500 GeV): θ_min=60 mrad < 100 mrad → mostly unconstrained")
    
    # Single lifetime
    print("\n" + "="*50)
    print("SINGLE LIFETIME ANALYSIS")
    print("="*50)
    
    lifetime = 100e-8
    df_results, event_df = process_with_acceptance(
        sample_csv, lifetime, geo_cache)
    
    hits = df_results[df_results['hits_tube']]
    print(f"\nτ = {lifetime*1e9:.1f} ns:")
    print(f"  Particles hitting fiducial: {len(hits)}")
    if len(hits) > 0:
        print(f"  Mean acceptance: {hits['acceptance'].mean():.4f}")
        print(f"  Mean P_decay (with cuts):    {hits['decay_probability'].mean():.6f}")
        print(f"  Mean P_decay (without cuts): {hits['decay_probability_no_cuts'].mean():.6f}")
    print(f"  Mean P(≥1) with cuts:    {event_df['prob_at_least_one_decays'].mean():.6e}")
    print(f"  Mean P(≥1) without cuts: {event_df['prob_at_least_one_no_cuts'].mean():.6e}")
    
    # --- Separation histogram ---
    print("\n" + "="*50)
    print("SEPARATION DISTRIBUTION")
    print("="*50)
    
    (seps, weights, momenta, pointing, p_soft, dca, vtx_in,
     open_angle, sep_outer, d_implied, collinearity) = sample_separations(
        geo_cache, lifetime, n_samples_per_particle=200)

    print(f"  Hit resolution: {HIT_RESOLUTION*1000:.1f} mm, "
          f"N layers: {N_LAYERS}, DCA cut: {DCA_CUT*100:.1f} cm")
    print(f"  Conditional max sep_outer: < {SEP_OUT_MAX_PARALLEL:.1f} m "
          f"when open_angle < {THETA_PARALLEL*1000:.0f} mrad")
    print(f"  IP-muon-transit veto: collinearity > {COLLIN_MIN*1000:.0f} mm "
          f"when sep_outer > {SEP_OUT_COLLIN_GATE*100:.0f} cm")
    
    if len(seps) > 0:
        fig_sep, axes_sep = plt.subplots(1, 3, figsize=(18, 5))
        
        ax = axes_sep[0]
        bins = np.linspace(0, 0.5, 100)
        ax.hist(seps, bins=bins, weights=weights, color='steelblue',
                edgecolor='black', linewidth=0.3, alpha=0.8)
        ax.axvline(SEP_MIN, color='red', linestyle='--', linewidth=2,
                   label=f'min sep = {SEP_MIN*1000:.0f} mm')
        ax.axvline(SEP_MAX, color='red', linestyle='-', linewidth=2,
                   label=f'max sep = {SEP_MAX*100:.0f} cm')
        ax.axvspan(SEP_MIN, SEP_MAX, color='green', alpha=0.1, label='Accepted')
        ax.set_xlabel('Separation at detector (m)')
        ax.set_ylabel('Weighted counts (decay prob.)')
        ax.set_title(f'Separation distribution (τ = {lifetime*1e9:.0f} ns)\n'
                     f'All decays, weighted by P(decay)')
        ax.legend(fontsize=9)
        ax.set_xlim(0, 0.5)
        
        ax2 = axes_sep[1]
        bins_log = np.logspace(-4, 1, 100)
        ax2.hist(seps, bins=bins_log, weights=weights, color='steelblue',
                 edgecolor='black', linewidth=0.3, alpha=0.8)
        ax2.axvline(SEP_MIN, color='red', linestyle='--', linewidth=2,
                    label=f'min sep = {SEP_MIN*1000:.0f} mm')
        ax2.axvline(SEP_MAX, color='red', linestyle='-', linewidth=2,
                    label=f'max sep = {SEP_MAX*100:.0f} cm')
        ax2.axvspan(SEP_MIN, SEP_MAX, color='green', alpha=0.1, label='Accepted')
        ax2.set_xscale('log')
        ax2.set_xlabel('Separation at detector (m)')
        ax2.set_ylabel('Weighted counts (decay prob.)')
        ax2.set_title('Log-scale separation\n(showing full range)')
        ax2.legend(fontsize=9)
        
        ax3 = axes_sep[2]
        mask_finite = np.isfinite(seps) & (seps > 0)
        h = ax3.hist2d(momenta[mask_finite], seps[mask_finite] * 100,
                       bins=[np.linspace(0, 500, 50), np.linspace(0, 500, 50)],
                       weights=weights[mask_finite],
                       cmap='viridis', cmin=1e-20)
        ax3.axhline(SEP_MIN * 100, color='red', linestyle='--', linewidth=2,
                    label=f'min = {SEP_MIN*1000:.0f} mm')
        ax3.axhline(SEP_MAX * 100, color='red', linestyle='-', linewidth=2,
                    label=f'max = {SEP_MAX*100:.0f} cm')
        ax3.set_xlabel('LLP momentum (GeV/c)')
        ax3.set_ylabel('Separation at detector (cm)')
        ax3.set_title('Separation vs LLP momentum')
        ax3.legend(fontsize=9, loc='upper right')
        plt.colorbar(h[3], ax=ax3, label='Weighted counts')
        
        plt.tight_layout()
        plt.savefig('separation_histogram'+outString+'.png', dpi=150)
        show_or_close()
        
        in_window = (seps >= SEP_MIN) & (seps <= SEP_MAX)
        frac_accepted = weights[in_window].sum() / weights.sum()
        print(f"  Fraction of decays in separation window: {frac_accepted:.3f}")
        print(f"  Median separation (all decays): {np.median(seps)*100:.1f} cm")
        accepted_seps = seps[in_window]
        if len(accepted_seps) > 0:
            print(f"  Median separation (accepted):  {np.median(accepted_seps)*100:.1f} cm")

    # --- Pointing angle distribution ---
    if len(pointing) > 0:
        print("\n" + "="*50)
        print("POINTING ANGLE DISTRIBUTION")
        print("="*50)

        # Apply separation acceptance to pointing angles
        in_window = (seps >= SEP_MIN) & (seps <= SEP_MAX)
        pointing_acc = pointing[in_window]
        weights_acc = weights[in_window]
        momenta_acc = momenta[in_window]

        pointing_mrad = pointing_acc * 1000

        fig_pt, axes_pt = plt.subplots(1, 3, figsize=(18, 5))

        ax = axes_pt[0]
        bins_pt = np.linspace(0, np.percentile(pointing_mrad, 99.5), 80)
        ax.hist(pointing_mrad, bins=bins_pt, weights=weights_acc,
                color='darkorange', edgecolor='black', linewidth=0.3, alpha=0.8)
        ax.set_xlabel('Pointing angle (mrad)')
        ax.set_ylabel('Weighted counts (decay prob.)')
        ax.set_title(f'Pointing angle: bisector vs LLP direction\n'
                     f'(τ = {lifetime*1e9:.0f} ns, accepted decays)')
        median_pt = np.median(pointing_mrad)
        ax.axvline(median_pt, color='red', linestyle='--', linewidth=2,
                   label=f'Median = {median_pt:.2f} mrad')
        ax.legend(fontsize=9)

        ax2 = axes_pt[1]
        bins_log_pt = np.logspace(np.log10(max(pointing_mrad.min(), 1e-3)),
                                  np.log10(pointing_mrad.max()), 80)
        ax2.hist(pointing_mrad, bins=bins_log_pt, weights=weights_acc,
                 color='darkorange', edgecolor='black', linewidth=0.3, alpha=0.8)
        ax2.set_xscale('log')
        ax2.set_xlabel('Pointing angle (mrad)')
        ax2.set_ylabel('Weighted counts (decay prob.)')
        ax2.set_title('Log-scale pointing angle')

        ax3 = axes_pt[2]
        mask_fin = np.isfinite(pointing_mrad) & (pointing_mrad > 0)
        h = ax3.hist2d(momenta_acc[mask_fin], pointing_mrad[mask_fin],
                       bins=[np.linspace(0, 500, 50),
                             np.linspace(0, np.percentile(pointing_mrad, 99), 50)],
                       weights=weights_acc[mask_fin],
                       cmap='viridis', cmin=1e-20)
        ax3.set_xlabel('LLP momentum (GeV/c)')
        ax3.set_ylabel('Pointing angle (mrad)')
        ax3.set_title('Pointing angle vs LLP momentum')
        plt.colorbar(h[3], ax=ax3, label='Weighted counts')

        plt.tight_layout()
        plt.savefig('pointing_angle'+outString+'.png', dpi=150)
        show_or_close()

        print(f"  Median pointing angle (accepted): {median_pt:.3f} mrad")
        p90 = np.percentile(pointing_mrad, 90)
        p99 = np.percentile(pointing_mrad, 99)
        print(f"  90th percentile: {p90:.3f} mrad")
        print(f"  99th percentile: {p99:.3f} mrad")

        # Convert to spatial resolution at the IP distance
        mean_entry = np.nanmean(geo_cache['entry_d'][geo_cache['hits']])
        print(f"  At mean entry distance ({mean_entry:.0f} m):")
        print(f"    Median miss distance: {median_pt * mean_entry:.1f} mm")
        print(f"    90th pct miss distance: {p90 * mean_entry:.1f} mm")

    # --- DCA distribution ---
    if len(dca) > 0:
        print("\n" + "="*50)
        print("DCA DISTRIBUTION")
        print("="*50)

        in_window = (seps >= SEP_MIN) & (seps <= SEP_MAX)
        dca_acc = dca[in_window]
        weights_dca = weights[in_window]
        momenta_dca = momenta[in_window]

        dca_cm = dca_acc * 100

        fig_dca, axes_dca = plt.subplots(1, 3, figsize=(18, 5))

        ax = axes_dca[0]
        pct99 = np.percentile(dca_cm, 99.5) if len(dca_cm) > 0 else 10
        bins_dca = np.linspace(0, max(pct99, DCA_CUT * 100 * 2), 80)
        ax.hist(dca_cm, bins=bins_dca, weights=weights_dca,
                color='seagreen', edgecolor='black', linewidth=0.3, alpha=0.8)
        ax.axvline(DCA_CUT * 100, color='red', linestyle='--', linewidth=2,
                   label=f'DCA cut = {DCA_CUT*100:.1f} cm')
        ax.set_xlabel('DCA (cm)')
        ax.set_ylabel('Weighted counts (decay prob.)')
        ax.set_title(f'Track DCA (τ = {lifetime*1e9:.0f} ns, '
                     f'σ_hit = {HIT_RESOLUTION*1000:.1f} mm)')
        ax.legend(fontsize=9)

        ax2 = axes_dca[1]
        bins_log_dca = np.logspace(np.log10(max(dca_cm[dca_cm > 0].min(), 1e-4))
                                   if np.any(dca_cm > 0) else -4,
                                   np.log10(max(dca_cm.max(), 1)), 80)
        ax2.hist(dca_cm, bins=bins_log_dca, weights=weights_dca,
                 color='seagreen', edgecolor='black', linewidth=0.3, alpha=0.8)
        ax2.axvline(DCA_CUT * 100, color='red', linestyle='--', linewidth=2,
                    label=f'DCA cut = {DCA_CUT*100:.1f} cm')
        ax2.set_xscale('log')
        ax2.set_xlabel('DCA (cm)')
        ax2.set_ylabel('Weighted counts (decay prob.)')
        ax2.set_title('DCA (log scale)')
        ax2.legend(fontsize=9)

        ax3 = axes_dca[2]
        mask_fin = np.isfinite(dca_cm) & (dca_cm > 0)
        if np.any(mask_fin):
            h = ax3.hist2d(momenta_dca[mask_fin], dca_cm[mask_fin],
                           bins=[np.linspace(0, 500, 50),
                                 np.linspace(0, max(pct99, DCA_CUT * 100 * 2), 50)],
                           weights=weights_dca[mask_fin],
                           cmap='viridis', cmin=1e-20)
            ax3.axhline(DCA_CUT * 100, color='red', linestyle='--', linewidth=2,
                        label=f'DCA cut = {DCA_CUT*100:.1f} cm')
            ax3.set_xlabel('LLP momentum (GeV/c)')
            ax3.set_ylabel('DCA (cm)')
            ax3.set_title('DCA vs LLP momentum')
            ax3.legend(fontsize=9, loc='upper right')
            plt.colorbar(h[3], ax=ax3, label='Weighted counts')

        plt.tight_layout()
        plt.savefig('dca_' + outString + '.png', dpi=150)
        show_or_close()

        median_dca = np.median(dca_cm)
        w_total = weights_dca.sum()
        frac_dca = weights_dca[dca_acc <= DCA_CUT].sum() / w_total \
            if w_total > 0 else 0
        vtx_in_acc = vtx_in[in_window]
        frac_vtx = weights_dca[vtx_in_acc].sum() / w_total \
            if w_total > 0 else 0
        print(f"  Median DCA: {median_dca:.4f} cm")
        print(f"  Fraction passing DCA < {DCA_CUT*100:.1f} cm: {frac_dca:.4f}")
        print(f"  Fraction with vertex in fiducial: {frac_vtx:.4f}")
        print(f"  (fractions above computed after separation cuts)")

    # --- sep_outer distribution, split by open-angle band ---
    # Motivates the conditional sep_outer<X-when-parallel cut:
    # nearly-parallel signal stays at small sep_outer; clearly diverging
    # signal can occupy arbitrary sep_outer (legitimately).
    if len(sep_outer) > 0:
        print("\n" + "="*50)
        print("SEP_OUTER BY OPEN-ANGLE BAND")
        print("="*50)

        in_window = (seps >= SEP_MIN) & (seps <= SEP_MAX)
        so   = sep_outer[in_window]
        op   = open_angle[in_window]
        ww   = weights[in_window]

        bands = [
            (f'open_angle < {THETA_PARALLEL*1000:.0f} mrad',
                 (None, THETA_PARALLEL), 'steelblue'),
            (f'{THETA_PARALLEL*1000:.0f} ≤ open_angle < 100 mrad',
                 (THETA_PARALLEL, 0.100), 'darkorange'),
            ('open_angle ≥ 100 mrad',
                 (0.100, None), 'crimson'),
        ]

        fig_so, axes_so = plt.subplots(1, 2, figsize=(14, 5))

        bins_lin = np.linspace(0, max(2.0, SEP_OUT_MAX_PARALLEL * 1.5), 80)
        positive = so[so > 0]
        lo = np.log10(positive.min()) if len(positive) else -3
        hi = np.log10(max(so.max(), 10.0))
        bins_log = np.logspace(lo, hi, 80)

        for ax, bins, scale in [(axes_so[0], bins_lin, 'linear'),
                                (axes_so[1], bins_log, 'log')]:
            for label, (lo_b, hi_b), color in bands:
                mask = np.ones(len(op), dtype=bool)
                if lo_b is not None: mask &= (op >= lo_b)
                if hi_b is not None: mask &= (op <  hi_b)
                if mask.sum() == 0:
                    continue
                w_sum = ww[mask].sum()
                ax.hist(so[mask], bins=bins, weights=ww[mask],
                        histtype='step', linewidth=2, color=color,
                        label=f'{label}  (w={w_sum:.2e})')
            ax.axvline(SEP_OUT_MAX_PARALLEL, color='gray', linestyle='--',
                       linewidth=1.5,
                       label=f'parallel max = {SEP_OUT_MAX_PARALLEL:.1f} m')
            ax.set_xlabel('sep_outer (m)')
            ax.set_ylabel('Weighted counts (decay prob.)')
            ax.set_xscale(scale)
            if scale == 'log':
                ax.set_yscale('log')
            ax.legend(fontsize=8, loc='upper right')
            ax.grid(True, which='both', alpha=0.3)

        axes_so[0].set_title(f'sep_outer by open-angle band  (τ = {lifetime*1e9:.0f} ns)')
        axes_so[1].set_title('Log-log view')
        plt.tight_layout()
        plt.savefig('sep_outer_bands_' + outString + '.png', dpi=150)
        show_or_close()

        # Numerical summary
        def wpct(x, w, q):
            if len(x) == 0: return 0.0
            idx = np.argsort(x); xs, ws = x[idx], w[idx]
            c = np.cumsum(ws); return xs[np.searchsorted(c, q*c[-1])]
        print(f"  {'band':<32} {'med':>10} {'90%':>9} {'frac>parallel_max':>20}")
        for label, (lo_b, hi_b), _ in bands:
            mask = np.ones(len(op), dtype=bool)
            if lo_b is not None: mask &= (op >= lo_b)
            if hi_b is not None: mask &= (op <  hi_b)
            if mask.sum() == 0:
                continue
            x_, w_ = so[mask], ww[mask]
            med = wpct(x_, w_, 0.50); p90 = wpct(x_, w_, 0.90)
            frac = w_[x_ > SEP_OUT_MAX_PARALLEL].sum() / w_.sum()
            print(f"  {label:<32} {med*100:>8.2f}cm {p90*100:>7.2f}cm "
                  f"{frac:>20.4f}")

    # --- Collinearity distribution (IP-muon-transit veto diagnostic) ---
    if len(collinearity) > 0:
        print("\n" + "="*50)
        print("COLLINEARITY (IP-muon-transit veto)")
        print("="*50)

        in_window = (seps >= SEP_MIN) & (seps <= SEP_MAX)
        co = collinearity[in_window]
        so = sep_outer[in_window]
        ww = weights[in_window]

        below = so <= SEP_OUT_COLLIN_GATE
        above = so >  SEP_OUT_COLLIN_GATE

        fig_co, axes_co = plt.subplots(1, 3, figsize=(18, 5))

        co_mm = co * 1000
        ax = axes_co[0]
        bins = np.linspace(0, max(np.percentile(co_mm, 99.5), COLLIN_MIN*1000*3), 80)
        for lbl, m, color in [
            (f'sep_out ≤ {SEP_OUT_COLLIN_GATE*100:.0f} cm (gate closed)',
                 below, 'steelblue'),
            (f'sep_out > {SEP_OUT_COLLIN_GATE*100:.0f} cm (gate open)',
                 above, 'crimson'),
        ]:
            if m.sum() == 0:
                continue
            ax.hist(co_mm[m], bins=bins, weights=ww[m],
                    histtype='step', linewidth=2, color=color,
                    label=f'{lbl}  (w={ww[m].sum():.2e})')
        ax.axvline(COLLIN_MIN*1000, color='gray', linestyle='--', linewidth=1.5,
                   label=f'cut = {COLLIN_MIN*1000:.0f} mm')
        # IP-muon reference (from realistic ray-cast MC, σ_hit=3mm, L=24cm)
        ax.axvspan(0, 8.2, color='gray', alpha=0.15,
                   label='IP-muon transit range (≤ 8 mm)')
        ax.set_xlabel('collinearity (mm)')
        ax.set_ylabel('Weighted counts (decay prob.)')
        ax.set_title(f'Collinearity by sep_outer gate (τ = {lifetime*1e9:.0f} ns)')
        ax.legend(fontsize=8, loc='upper right')
        ax.grid(True, alpha=0.3)

        ax2 = axes_co[1]
        positive = co_mm[co_mm > 0]
        lo = np.log10(max(positive.min(), 0.1)) if len(positive) else -1
        hi = np.log10(max(co_mm.max(), 100.0))
        bins_log = np.logspace(lo, hi, 80)
        for lbl, m, color in [
            ('gate closed', below, 'steelblue'),
            ('gate open',   above, 'crimson'),
        ]:
            if m.sum() == 0:
                continue
            ax2.hist(co_mm[m], bins=bins_log, weights=ww[m],
                     histtype='step', linewidth=2, color=color, label=lbl)
        ax2.axvline(COLLIN_MIN*1000, color='gray', linestyle='--', linewidth=1.5)
        ax2.axvspan(0.1, 8.2, color='gray', alpha=0.15,
                    label='IP-muon range')
        ax2.set_xscale('log'); ax2.set_yscale('log')
        ax2.set_xlabel('collinearity (mm)')
        ax2.set_ylabel('Weighted counts')
        ax2.set_title('Log-log view')
        ax2.legend(fontsize=8, loc='upper left')
        ax2.grid(True, which='both', alpha=0.3)

        ax3 = axes_co[2]
        # 2D heat map: collinearity vs sep_outer with rejection region shaded
        mask_fin = np.isfinite(co) & np.isfinite(so) & (co > 0) & (so > 0)
        h = ax3.hist2d(so[mask_fin]*100, co_mm[mask_fin],
                       bins=[np.logspace(0, 3, 50), np.logspace(-1, 3, 50)],
                       weights=ww[mask_fin], cmap='viridis', cmin=1e-30)
        ax3.axvline(SEP_OUT_COLLIN_GATE*100, color='red', linestyle='--',
                    linewidth=1.5,
                    label=f'gate = {SEP_OUT_COLLIN_GATE*100:.0f} cm')
        ax3.axhline(COLLIN_MIN*1000, color='red', linestyle='--', linewidth=1.5,
                    label=f'cut = {COLLIN_MIN*1000:.0f} mm')
        # Shade rejection region (sep_out > gate AND collin < cut)
        ax3.fill_between([SEP_OUT_COLLIN_GATE*100, 1e3], 1e-1, COLLIN_MIN*1000,
                         color='red', alpha=0.15, label='reject (muon-like)')
        ax3.set_xscale('log'); ax3.set_yscale('log')
        ax3.set_xlabel('sep_outer (cm)')
        ax3.set_ylabel('collinearity (mm)')
        ax3.set_title('Conditional cut region')
        ax3.legend(fontsize=8, loc='lower right')
        plt.colorbar(h[3], ax=ax3, label='Weighted counts')

        plt.tight_layout()
        plt.savefig('collinearity_' + outString + '.png', dpi=150)
        show_or_close()

        # Numerical summary
        gated = above
        if gated.sum() > 0:
            x_, w_ = co[gated], ww[gated]
            def wpct(x, w, q):
                idx = np.argsort(x); xs, ws = x[idx], w[idx]
                c = np.cumsum(ws); return xs[np.searchsorted(c, q*c[-1])]
            med = wpct(x_, w_, 0.50); p1 = wpct(x_, w_, 0.01)
            frac_fail = w_[x_ < COLLIN_MIN].sum() / w_.sum()
            print(f"  Above gate (sep_out > {SEP_OUT_COLLIN_GATE*100:.0f} cm): "
                  f"median = {med*1000:.1f} mm, 1% = {p1*1000:.1f} mm")
            print(f"  Signal failing collinearity > {COLLIN_MIN*1000:.0f} mm "
                  f"in this band: {frac_fail*100:.3f}%")
        else:
            print(f"  No events above the gate (sep_out > "
                  f"{SEP_OUT_COLLIN_GATE*100:.0f} cm) for this sample.")

    # --- Cutflow table ---
    if len(seps) > 0:
        print("\n" + "="*50)
        print("CUTFLOW TABLE")
        print("="*50)

        cutflow = build_cutflow(seps, pointing, weights, momenta, p_soft,
                                dca, vtx_in, open_angle, sep_outer,
                                collinearity)
        cutflow_df = pd.DataFrame(cutflow)
        print(f"\n{'Cut':<25} {'Eff (cumul.)':<15} {'Eff (margin.)':<15}")
        print("-" * 55)
        for row in cutflow:
            print(f"{row['cut']:<25} {row['efficiency']:<15.4f} "
                  f"{row['marginal_efficiency']:<15.4f}")

        cutflow_df.to_csv('cutflow_' + outString + '.csv', index=False)
        print(f"\nCutflow saved to cutflow_{outString}.csv")

    # --- Save MC distributions for overlay plotting ---
    if len(seps) > 0:
        np.savez('mc_distributions_' + outString + '.npz',
                 seps=seps, pointing=pointing, weights=weights,
                 momenta=momenta, p_soft=p_soft,
                 dca=dca, vtx_in=vtx_in, open_angle=open_angle,
                 sep_outer=sep_outer, d_implied=d_implied,
                 collinearity=collinearity,
                 lifetime=lifetime, label=outString,
                 p_cut=P_CUT, sep_min=SEP_MIN, sep_max=SEP_MAX,
                 dca_cut=DCA_CUT,
                 theta_parallel=THETA_PARALLEL,
                 sep_out_max_parallel=SEP_OUT_MAX_PARALLEL,
                 collin_min=COLLIN_MIN,
                 sep_out_collin_gate=SEP_OUT_COLLIN_GATE,
                 hit_resolution=HIT_RESOLUTION, n_layers=N_LAYERS)
        print(f"MC distributions saved to mc_distributions_{outString}.npz")

    # Lifetime scan
    print("\n" + "="*50)
    print("LIFETIME SCAN")
    print("="*50)
    
    lifetimes = np.logspace(-10.5, -3.5, 20)
    scan = analyze_decay_vs_lifetime(sample_csv, geo_cache, lifetimes)
    
    # === Plotting ===
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    ax1 = axes[0, 0]
    ax1.semilogx(lifetimes * 1e9, scan['mean_acceptance'],
                 'b-', linewidth=2)
    ax1.set_xlabel('Lifetime (ns)')
    ax1.set_ylabel('Mean Acceptance')
    ax1.set_title(f'Two-body acceptance\n(p > {P_CUT*1000:.0f} MeV/c, '
                  f'{SEP_MIN*1000:.0f} mm < sep < {SEP_MAX*100:.0f} cm)')
    ax1.grid(True, which="both", ls="-", alpha=0.2)
    ax1.set_ylim(0, 1.05)
    
    ax2 = axes[0, 1]
    ax2.loglog(lifetimes * 1e9, scan['mean_at_least_one_decay_prob'],
               'r-', linewidth=2, label='With cuts')
    ax2.loglog(lifetimes * 1e9, scan['mean_at_least_one_no_cuts'],
               'r--', linewidth=2, alpha=0.5, label='Without cuts')
    ax2.set_xlabel('Lifetime (ns)')
    ax2.set_ylabel('Mean P(≥1 decay)')
    ax2.set_title('Effect of analysis cuts')
    ax2.grid(True, which="both", ls="-", alpha=0.2)
    ax2.legend()
    
    ax3 = axes[1, 0]
    ax3.loglog(lifetimes * 1e9, scan['mean_single_particle_decay_prob'],
               'b-', linewidth=2, label='With cuts')
    ax3.loglog(lifetimes * 1e9, scan['mean_single_no_cuts'],
               'b--', linewidth=2, alpha=0.5, label='Without cuts')
    ax3.set_xlabel('Lifetime (ns)')
    ax3.set_ylabel('Mean Decay Probability')
    ax3.set_title('Single Particle: With vs Without Cuts')
    ax3.grid(True, which="both", ls="-", alpha=0.2)
    ax3.legend()
    
    ax4 = axes[1, 1]
    ax4.loglog(lifetimes * SPEED_OF_LIGHT, scan['exclusion'],
               color='blue', linewidth=2, label="GRENDEL")
    # ax4.loglog(lifetimes * SPEED_OF_LIGHT, scan['exclusion_no_cuts'],
    #            color='blue', linewidth=2, linestyle='--', alpha=0.5,
    #            label="milliQan (no cuts)")
    ax4.set_xlabel(r'$c\tau$ (m)')
    ax4.set_ylabel('BR')
    ax4.grid(True, which="both", ls="-", alpha=0.2)
    
    ext = {}
    if df_results['mass'].iloc[0] == 15:
        ext["MATHUSLA"] = np.loadtxt("external/MATHUSLA.csv", delimiter=",")
        ext["CODEX"] = np.loadtxt("external/CODEX.csv", delimiter=",")
        ext["CMS"] = np.loadtxt("external/CMS.csv", delimiter=",")
        # ext["ANUBISOpt"] = np.loadtxt("external/ANUBISOpt.csv", delimiter=",")
        ext["ANUBISCons"] = np.loadtxt("external/ANUBISUpdateCons.csv", delimiter=",")
    
        ax4.loglog(ext["MATHUSLA"][:, 0], ext["MATHUSLA"][:, 1],
                   color="green", linewidth=2, label="MATHUSLA")
        ax4.loglog(ext["CODEX"][:, 0], ext["CODEX"][:, 1],
                   color="cyan", linewidth=2, label="CODEX-b")
        ax4.loglog(ext["CMS"][:, 0], ext["CMS"][:, 1],
                   color="purple", linewidth=2, label="CMS")
        # ax4.loglog(ext["ANUBISOpt"][:, 0], ext["ANUBISOpt"][:, 1],
        #            color="purple", linewidth=2, linestyle="--", label="ANUBIS Opt")
        ax4.loglog(ext["ANUBISCons"][:, 0], ext["ANUBISCons"][:, 1],
                   color="magenta", linewidth=2, linestyle="--", label="ANUBIS Cons")
        ax4.legend(fontsize=8, loc='lower right')
        ax4.set_ylim([1E-5,1])
    elif df_results['mass'].iloc[0] == 0.5:
        ext["CODEX"] = np.loadtxt("external/CODEX0p5.csv", delimiter=",")
    
        ax4.loglog(ext["CODEX"][:, 0], ext["CODEX"][:, 1],
                   color="cyan", linewidth=2, label="CODEX-b")
        ax4.legend(fontsize=8, loc='lower right')
        
    plt.tight_layout()
    plt.savefig('exclusion_2body'+outString+'.png', dpi=150)
    show_or_close()
    
    df_results.to_csv("particle_decay_results_2body.csv", index=False)
    event_df.to_csv("event_decay_statistics_2body.csv", index=False)
    print("\nResults saved.")
    print("Plots: exclusion_2body"+outString+".png, separation_histogram"+outString+".png")
