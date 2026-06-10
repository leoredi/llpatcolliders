"""
Cosmic-muon decay-in-flight background for GRENDEL.
====================================================

Models the background where an atmospheric muon enters the tunnel from
above (standard cos^2 theta_zenith angular distribution), decays in
flight inside the fiducial air volume (mu -> e nu nu), and the incoming
muon track + outgoing electron track form a two-track vertex that can
fake the LLP -> e+e- signal.

This complements:
  - cosmic_muon_check.py    : muon *transit* (straight-through, 4 collinear hits)
  - background_trident_*.py : muon *trident* (mu N -> mu N e+e- in air)

Physics / modelling choices (chosen with the user):
  - Direction: zenith from f(theta) = 3 cos^2(theta) sin(theta) on
    [0, theta_max], azimuth uniform; downgoing (CMS Y up, tunnel at +22 m).
  - Momentum: pluggable input spectrum, default exponential in |p| with
    mean 1 GeV/c, defined *at the cavern* (no overburden / energy loss).
  - Decay: full Michel electron-energy spectrum n(x) ~ x^2 (3 - 2x),
    isotropic rest-frame direction (unpolarised muon), Lorentz-boosted to
    the lab to get the electron lab momentum and the kink direction.
  - Reconstruction: the muon is measured at the wall it *enters* (two
    layers spanning DETECTOR_THICKNESS), the electron at the wall it
    *exits*; both segments are extrapolated to a common vertex and the
    full signal selection from decayProbPerEvent_2body is applied.
  - Normalisation: absolute rate (Hz and per livetime) from a configurable
    cosmic vertical intensity I_v [m^-2 s^-1 sr^-1]. The toy momentum
    spectrum sets only the shape; I_v sets the absolute scale.

Usage:
    python cosmic_decay_check.py [--n-rays 2000000] [--p-mean 1.0]
                                 [--spectrum exp] [--i-vertical 70.0]
                                 [--livetime-years 10] [--seed 42]
"""

import os
import sys
import argparse
import numpy as np
import matplotlib

from grendel_geometry import (
    mesh_fiducial, DETECTOR_THICKNESS, SPEED_OF_LIGHT, points_on_tracker,
    points_in_fiducial, classify_points_with_basis,
)
import reco_common
from decayProbPerEvent_2body import (
    P_CUT, SEP_MIN, SEP_MAX, DCA_CUT,
    THETA_PARALLEL, SEP_OUT_MAX_PARALLEL,
    COLLIN_MIN, SEP_OUT_COLLIN_GATE, COLLIN_FRAC, SEP_OUT_GATE,
    POINT_TIGHT_SEP_IN, SEP_IN_POINT_GATE, POINT_GLOBAL,
    HIT_RESOLUTION, _first_forward_hit,
)

# ------------------------------------------------------------------
# Physical constants
# ------------------------------------------------------------------
M_MUON     = 0.10566        # GeV/c^2
M_ELECTRON = 0.000511       # GeV/c^2
TAU_MUON   = 2.1969811e-6   # s
CTAU_MUON  = SPEED_OF_LIGHT * TAU_MUON   # ~658.6 m

# Default cosmic vertical intensity (muons above ~1 GeV at sea level).
# dI/dOmega = I_VERTICAL * cos^2(theta).  Units: m^-2 s^-1 sr^-1.
# Only sets the relative cos(theta) weighting; the absolute scale is fixed by
# anchoring the total muon-through-tunnel rate (see --muon-rate-hz).
I_VERTICAL_DEFAULT = 70.0

# Total cosmic-muon rate through the fiducial volume (Hz). The user supplies
# this from their separate overburden treatment; all weights are rescaled so
# the "muons through fiducial" stage equals this value.
MUON_RATE_HZ_DEFAULT = 600.0

# Minimum reconstructed-track momentum (GeV/c). Overrides the signal's
# 600 MeV electron cut for this cosmic-decay study.
P_SOFT_CUT = 0.010

# Targeted (L-scaled) collinearity veto: when sep_out > L require
# collinearity > COLLIN_FRAC * L (L = DETECTOR_THICKNESS). COLLIN_FRAC is
# imported from decayProbPerEvent_2body (single source, shared with the signal);
# --collin-frac overrides it for a given cosmic run.

# Per-hit timing resolution (s) and the timing-consistency chi^2 cut.
# The selection asks that the 4 hit times are consistent with particles
# leaving the reconstructed vertex and travelling at c (one free parameter,
# the vertex time t0; ndof = 3).
SIGMA_T_DEFAULT = 1.5e-9      # 1.5 ns
CHI2_TIMING_MAX = 9.0         # ndof = 3 (chi2(3) CDF ~ 0.97); the APPLIED final cut

# Per-track outgoing-velocity timing test (reported as a side-by-side
# comparison in the cutflow printout; NOT applied as a cut). For each
# track the two hits must be consistent with motion AWAY from the vertex at c:
#   pull = [(t_far - t_near) - (R_far - R_near)/c] / (sqrt2 * sigma_t)
# An incoming muon track gives a large NEGATIVE pull (its far hit is earlier).
# Require BOTH tracks pull > -VELO_NCUT. Signal daughters are outgoing
# (pull ~ N(0,1)); n_cut = 2.17 -> per-event signal eff ~0.97 (matches chi2<9).
VELO_NCUT = 2.17

# The tracker-layer spacing is NOT a separate knob: it equals DETECTOR_THICKNESS
# (inner layer at the fiducial face, outer at the tunnel wall). This is the same
# single constant the fiducial mesh is built from, so the mesh and the
# reconstruction can never disagree. To study a different spacing, change
# DETECTOR_THICKNESS in grendel_geometry.py (one place).

IP = np.array([0.0, 0.0, 0.0])


# ------------------------------------------------------------------
# Input momentum spectrum (pluggable)
# ------------------------------------------------------------------

def sample_momentum(n, rng, mean=1.0, kind='exp'):
    """
    Sample cosmic-muon momenta |p| (GeV/c) at the cavern.

    kind='exp'    : f(p) ~ exp(-p / mean)               (default, mean = `mean`)
    kind='powerlaw': f(p) ~ p^-2.7 above p_min=`mean`   (toy hard spectrum)
    kind='mono'   : delta function at p = `mean`
    """
    if kind == 'exp':
        return rng.exponential(mean, n)
    if kind == 'mono':
        return np.full(n, mean)
    if kind == 'powerlaw':
        # inverse-CDF for p^-2.7 on [mean, 200 GeV]
        a = 2.7
        lo, hi = mean, 200.0
        u = rng.uniform(0, 1, n)
        return ((1 - u) * lo**(1 - a) + u * hi**(1 - a)) ** (1.0 / (1 - a))
    raise ValueError(f"unknown spectrum kind: {kind}")


def sample_depth_momentum(cos_z, theta_max, rng, e_min=0.3):
    """
    Sample per-muon momenta |p| (GeV/c) from the realistic *depth* spectrum
    (Gaisser surface flux propagated through the rock overburden), using each
    muon's own zenith angle. Uses CosmicMuonFlux from background_trident_update.

    A (cos_theta, E) inverse-CDF table is built once and inverted per muon.
    """
    from background_trident_update import CosmicMuonFlux
    flux = CosmicMuonFlux()

    E_grid = np.geomspace(e_min, 2000.0, 400)
    cz_grid = np.linspace(np.cos(theta_max), 1.0, 60)
    cdf_table = np.empty((len(cz_grid), len(E_grid)))
    for i, cz in enumerate(cz_grid):
        pdf = np.array([flux.flux_at_depth(E, cz) for E in E_grid])
        pdf = np.maximum(pdf, 0.0)
        c = np.cumsum(pdf)
        cdf_table[i] = (c / c[-1]) if c[-1] > 0 else np.linspace(0, 1, len(E_grid))

    idx = np.clip(np.searchsorted(cz_grid, cos_z), 0, len(cz_grid) - 1)
    u = rng.uniform(0, 1, len(cos_z))
    E = np.empty(len(cos_z))
    for i in range(len(cz_grid)):
        m = idx == i
        if m.any():
            E[m] = np.interp(u[m], cdf_table[i], E_grid)
    return np.sqrt(np.maximum(E**2 - M_MUON**2, 1e-6))


# Path to a ROOT file holding the measured cavern-entry muon |p| spectrum
# (TH1D 'h_muon_momentum_at_cavern_entry', MeV). Set by main() from --momentum-root.
EMPIRICAL_ROOT_PATH = 'muon_entry_distributions.root'
# Optional generation window [GeV]: restrict muon |p| sampling to this range
# (importance sampling) and scale the rate by the window's spectral fraction.
# None => full spectrum. Set by main() from --gen-pmin/--gen-pmax.
GEN_PMIN = None
GEN_PMAX = None
_EMP_CACHE = {}


def _load_empirical(path):
    """Load (and cache) the measured cavern-entry momentum histogram.
    Returns dict with resolved-bin edges_GeV, per-bin counts, the overflow
    count (>10 GeV, shapeless) and the total (resolved + overflow)."""
    if path in _EMP_CACHE:
        return _EMP_CACHE[path]
    import uproot
    h = uproot.open(path)['h_muon_momentum_at_cavern_entry']
    c = h.values(flow=True)                  # [underflow, ...bins..., overflow]
    edges = h.axis().edges() / 1000.0        # MeV -> GeV, len = nbins+1
    counts = c[1:-1].astype(float)
    overflow = float(c[-1])
    res = dict(edges=edges, counts=counts, overflow=overflow,
               total=float(counts.sum() + overflow))
    _EMP_CACHE[path] = res
    return res


def empirical_window_fraction(p_min, p_max, path=None):
    """Fraction of the full measured spectrum within [p_min, p_max] GeV
    (resolved range only; assumes p_max <= the histogram's upper edge)."""
    d = _load_empirical(path or EMPIRICAL_ROOT_PATH)
    lo, hi = d['edges'][:-1], d['edges'][1:]
    overlap = np.clip(np.minimum(hi, p_max) - np.maximum(lo, p_min), 0.0, None)
    w = float((d['counts'] * overlap / (hi - lo)).sum())
    return w / d['total']


def sample_empirical_momentum(n, rng, path=None, tail_index=2.7, tail_max=2000.0):
    """Sample muon |p| (GeV/c) from the measured cavern-entry spectrum.

    Full-spectrum mode (GEN_PMIN/GEN_PMAX both None): the resolved range
    (<= 10 GeV, ~4.5% of muons) is drawn from the histogram by inverse-CDF;
    the unresolved overflow (>10 GeV, ~95.5%) from a p^-tail_index tail. Those
    high-betagamma muons give forward/collinear electrons removed by the
    collinearity cut, so the FINAL selected background is insensitive to the
    tail shape; only the low-p (large-kink) tail, which is resolved, matters.

    Restricted mode (GEN_PMIN/GEN_PMAX set): draw only within the window from
    the resolved histogram (importance sampling); the caller scales the rate by
    empirical_window_fraction() to recover the absolute normalisation.
    """
    d = _load_empirical(path or EMPIRICAL_ROOT_PATH)
    edges, counts = d['edges'], d['counts']
    lo, hi = edges[:-1], edges[1:]

    if GEN_PMIN is not None or GEN_PMAX is not None:
        pmin = 0.0 if GEN_PMIN is None else GEN_PMIN
        pmax = edges[-1] if GEN_PMAX is None else GEN_PMAX
        win = np.clip(np.minimum(hi, pmax) - np.maximum(lo, pmin), 0.0, None)
        wbin = counts * win / (hi - lo)              # expected counts in window
        if wbin.sum() <= 0:
            return np.full(n, 0.5 * (pmin + pmax))
        cdf = np.cumsum(wbin) / wbin.sum()
        bi = np.clip(np.searchsorted(cdf, rng.uniform(0, 1, n)), 0, len(lo) - 1)
        blo = np.maximum(lo[bi], pmin)
        bhi = np.minimum(hi[bi], pmax)
        return np.maximum(blo + rng.uniform(0, 1, n) * (bhi - blo), 1e-3)

    # full spectrum: resolved histogram + power-law overflow tail
    f_inrange = counts.sum() / d['total']
    p = np.empty(n)
    is_in = rng.uniform(0, 1, n) < f_inrange
    nin = int(is_in.sum())
    if nin:
        cdf = np.cumsum(counts) / counts.sum()
        bi = np.clip(np.searchsorted(cdf, rng.uniform(0, 1, nin)), 0, len(lo) - 1)
        p[is_in] = lo[bi] + rng.uniform(0, 1, nin) * (hi - lo)[bi]
    nt = n - nin
    if nt:
        a = tail_index
        ut = rng.uniform(0, 1, nt)
        p[~is_in] = ((1 - ut) * edges[-1]**(1 - a)
                     + ut * tail_max**(1 - a)) ** (1.0 / (1 - a))
    return np.maximum(p, 1e-3)


# ------------------------------------------------------------------
# Michel decay (mu -> e nu nu) electron kinematics
# ------------------------------------------------------------------

def sample_michel_x(n, rng):
    """
    Sample reduced electron energy x = E_e / E_max from the Michel
    spectrum n(x) ~ x^2 (3 - 2x) on [0, 1] (massless-e approx, unpolarised).
    Rejection sampling; the envelope max on [0,1] is n(1) = 1.
    """
    out = np.empty(n)
    filled = 0
    while filled < n:
        m = n - filled
        x = rng.uniform(0, 1, m)
        y = rng.uniform(0, 1, m)
        acc = y < x**2 * (3 - 2 * x)
        k = int(acc.sum())
        out[filled:filled + k] = x[acc]
        filled += k
    return out


def transverse_basis(d):
    """Two orthonormal vectors spanning the plane perpendicular to unit
    rows of d (shape (N,3)). Returns (e1, e2), each (N,3)."""
    helper = np.tile(np.array([0., 1., 0.]), (len(d), 1))
    par = np.abs((d * helper).sum(axis=1)) > 0.95
    helper[par] = np.array([1., 0., 0.])
    e1 = np.cross(d, helper)
    e1 /= np.linalg.norm(e1, axis=1, keepdims=True)
    e2 = np.cross(d, e1)
    return e1, e2


def michel_electron(p_mu, mu_dir, rng):
    """
    Given muon momentum p_mu (N,) and unit direction mu_dir (N,3), sample
    the lab-frame electron from a full Michel decay + Lorentz boost.

    Returns (p_e, dir_e): electron momentum magnitude (N,) and unit
    direction (N,3) in the lab.
    """
    n = len(p_mu)
    E_max = (M_MUON**2 + M_ELECTRON**2) / (2 * M_MUON)   # ~m_mu/2
    x = sample_michel_x(n, rng)
    E_star = np.maximum(x * E_max, M_ELECTRON)
    p_star = np.sqrt(np.maximum(E_star**2 - M_ELECTRON**2, 0.0))

    cos_ts = rng.uniform(-1, 1, n)
    sin_ts = np.sqrt(np.maximum(1 - cos_ts**2, 0.0))
    phi_s = rng.uniform(0, 2 * np.pi, n)

    E_mu = np.sqrt(p_mu**2 + M_MUON**2)
    gamma = E_mu / M_MUON
    beta = p_mu / E_mu

    # Boost along the muon direction (z* = mu_dir)
    p_long = gamma * (p_star * cos_ts + beta * E_star)
    p_perp = p_star * sin_ts

    e1, e2 = transverse_basis(mu_dir)
    p_vec = (p_long[:, None] * mu_dir
             + (p_perp * np.cos(phi_s))[:, None] * e1
             + (p_perp * np.sin(phi_s))[:, None] * e2)
    p_e = np.linalg.norm(p_vec, axis=1)
    dir_e = p_vec / p_e[:, None]
    return p_e, dir_e


# ------------------------------------------------------------------
# Cosmic ray sampling (cos^2 theta), with flux normalisation
# ------------------------------------------------------------------

def sample_cosmics(n_rays, y_throw, theta_max, rng):
    """
    Throw cosmic muons from a horizontal plane at Y = y_throw.

    Direction PDF (intensity measure): g(theta) ~ cos^2(theta) sin(theta)
    on [0, theta_max], azimuth uniform.  CDF -> cos theta sampling.

    Returns
    -------
    origins : (N,3)  throw positions on the plane
    dirs    : (N,3)  downgoing unit directions
    cos_z   : (N,)   cos(zenith)
    area    : float  throw-plane area (m^2)
    """
    u = rng.uniform(0, 1, n_rays)
    cos_min = np.cos(theta_max)
    # truncated inverse CDF of g ~ cos^2 sin on [0, theta_max]
    cos_z = ((1 - u) * (1 - cos_min**3) + cos_min**3) ** (1.0 / 3.0)
    sin_z = np.sqrt(np.maximum(1 - cos_z**2, 0.0))
    phi = rng.uniform(0, 2 * np.pi, n_rays)
    dx = sin_z * np.cos(phi)
    dy = -cos_z                       # downgoing (CMS Y up)
    dz = sin_z * np.sin(phi)
    dirs = np.column_stack([dx, dy, dz])

    bbox = mesh_fiducial.bounds
    pad = (y_throw - bbox[1, 1]) * np.tan(theta_max) + 5.0
    x_min, x_max = bbox[0, 0] - pad, bbox[1, 0] + pad
    z_min, z_max = bbox[0, 2] - pad, bbox[1, 2] + pad
    origins = np.column_stack([
        rng.uniform(x_min, x_max, n_rays),
        np.full(n_rays, y_throw),
        rng.uniform(z_min, z_max, n_rays),
    ])
    area = (x_max - x_min) * (z_max - z_min)
    return origins, dirs, cos_z, area


def traverse_fiducial(origins, dirs):
    """
    Ray-cast all cosmics against the fiducial mesh; return entry/exit
    distances for those crossing it (>=2 intersections).

    Returns boolean mask `hit` (N,) and arrays d_in, d_out (N,) with NaN
    where there is no crossing.
    """
    n = len(origins)
    d_in = np.full(n, np.nan)
    d_out = np.full(n, np.nan)
    locs, ray_idx, _ = mesh_fiducial.ray.intersects_location(
        ray_origins=origins, ray_directions=dirs, multiple_hits=True)
    if len(locs) == 0:
        return np.zeros(n, bool), d_in, d_out
    signed = np.einsum('ij,ij->i', locs - origins[ray_idx], dirs[ray_idx])
    fwd = signed > 1e-9
    locs, ray_idx, signed = locs[fwd], ray_idx[fwd], signed[fwd]
    # min/max signed distance per ray
    order = np.argsort(ray_idx)
    ray_idx, signed = ray_idx[order], signed[order]
    uniq, start = np.unique(ray_idx, return_index=True)
    for k, r in enumerate(uniq):
        end = start[k + 1] if k + 1 < len(start) else len(signed)
        seg = signed[start[k]:end]
        if len(seg) >= 2:
            d_in[r] = seg.min()
            d_out[r] = seg.max()
    hit = ~np.isnan(d_in)
    return hit, d_in, d_out


# ------------------------------------------------------------------
# Reconstruction of the muon+electron vertex
# ------------------------------------------------------------------

def line_closest_approach(a1, d1, a2, d2):
    """Per-row closest-approach midpoint of two lines a + t d (all (N,3))."""
    w0 = a1 - a2
    a = np.einsum('ij,ij->i', d1, d1)
    b = np.einsum('ij,ij->i', d1, d2)
    c = np.einsum('ij,ij->i', d2, d2)
    d = np.einsum('ij,ij->i', d1, w0)
    e = np.einsum('ij,ij->i', d2, w0)
    denom = a * c - b * b
    denom = np.where(np.abs(denom) < 1e-30, 1e-30, denom)
    t = (b * e - c * d) / denom
    s = (a * e - b * d) / denom
    p1 = a1 + t[:, None] * d1
    p2 = a2 + s[:, None] * d2
    return 0.5 * (p1 + p2)


def reconstruct(mu_ori, mu_dir, d_in, d_decay, dir_e, d_out_e,
                p_mu, p_e, sigma_hit, sigma_t, rng):
    """
    Build the 4 smeared hits (muon entry pair + electron exit pair),
    reconstruct all signal observables. Returns a dict of (N,) arrays.
    """
    n = len(d_in)
    L = DETECTOR_THICKNESS        # layer spacing = detector shell (single source)
    V = mu_ori + d_decay[:, None] * mu_dir

    # True hit positions. The layers are separated by L *radially* (inner at
    # the fiducial face, outer at the tunnel wall), so the along-track pair
    # spacing is L / |d . n_hat| -- the same convention as
    # reco_common.wall_inner_outer used by the signal. Plain-L spacing would
    # systematically shorten the stubs of oblique crossings, shrinking the
    # 4-hit collinearity of kinked vertices and underestimating the
    # background surviving the collinearity veto. Grazing tracks
    # (|d . n_hat| < 1e-6) are flagged unreconstructable, mirroring the NaN
    # behaviour of wall_inner_outer.
    P_ei = mu_ori + (d_in)[:, None] * mu_dir       # entry-inner (fiducial face)
    P_xi = V + (d_out_e)[:, None] * dir_e          # exit-inner  (fiducial face)

    def _along_track_spacing(P, d):
        theta, _, _, right, up = classify_points_with_basis(P)
        n_hat = np.cos(theta)[:, None] * right + np.sin(theta)[:, None] * up
        dn = np.abs(np.einsum('ij,ij->i', d, n_hat))
        return np.where(dn > 1e-6, L / np.maximum(dn, 1e-6), np.nan)

    s_mu = _along_track_spacing(P_ei, mu_dir)
    s_e = _along_track_spacing(P_xi, dir_e)
    ok_pair = np.isfinite(s_mu) & np.isfinite(s_e)
    s_mu = np.where(ok_pair, s_mu, L)              # placeholder; excluded below
    s_e = np.where(ok_pair, s_e, L)
    P_eo = P_ei - s_mu[:, None] * mu_dir           # entry-outer (tunnel wall)
    P_xo = P_xi + s_e[:, None] * dir_e             # exit-outer  (tunnel wall)

    # Geometric reconstruction via the shared single-source routine (same code
    # the signal uses). *_out = farther from the vertex: muon entry-outer and
    # electron exit-outer. Returns geometry + the smeared hits used downstream.
    geo = reco_common.reconstruct_3d(P_eo, P_ei, P_xi, P_xo, sigma_hit, rng)
    sep_in, sep_out = geo['sep'], geo['sep_outer']
    open_a, dca, collin = geo['open_angle'], geo['dca'], geo['collin']
    vtx_in, pointing, V_reco = geo['vtx_in'], geo['pointing'], geo['V_reco']
    H_eo, H_ei = geo['H_out1'], geo['H_in1']
    H_xi, H_xo = geo['H_in2'], geo['H_out2']

    # On-tracker: muon must enter through a tracker wall AND electron must
    # exit through a tracker wall (scintillator walls give no tracks), and
    # neither crossing may be grazing (ill-defined layer pair).
    on_tracker = (points_on_tracker(P_ei) & points_on_tracker(P_xi) & ok_pair)

    p_soft = np.minimum(p_mu, p_e)

    # ---- Timing selection ----
    # True hit times relative to the decay (t=0 at the vertex): the muon
    # reaches the vertex from outside (entry hits at NEGATIVE time), the
    # electron leaves the vertex outward (exit hits at POSITIVE time). Each
    # particle travels at its own beta. Times are then smeared by sigma_t.
    c = SPEED_OF_LIGHT
    beta_mu = p_mu / np.sqrt(p_mu**2 + M_MUON**2)
    beta_e = p_e / np.sqrt(p_e**2 + M_ELECTRON**2)
    r_eo = np.linalg.norm(P_eo - V, axis=1)
    r_ei = np.linalg.norm(P_ei - V, axis=1)
    r_xi = np.linalg.norm(P_xi - V, axis=1)
    r_xo = np.linalg.norm(P_xo - V, axis=1)
    t_true = np.stack([-r_eo / (beta_mu * c), -r_ei / (beta_mu * c),
                       r_xi / (beta_e * c), r_xo / (beta_e * c)], axis=1)
    t_meas = t_true + rng.normal(0, sigma_t, t_true.shape)

    # Hypothesis: all 4 hits come from particles leaving the reconstructed
    # vertex at c. Predicted time = t0 + R_i/c with R_i the measured-hit to
    # reco-vertex distance. Fit the single t0 (= mean residual) and form chi2.
    H = np.stack([H_eo, H_ei, H_xi, H_xo], axis=1)
    R = np.linalg.norm(H - V_reco[:, None, :], axis=2)
    x = t_meas - R / c
    resid = x - x.mean(axis=1, keepdims=True)
    timing_chi2 = (resid**2).sum(axis=1) / sigma_t**2          # ndof = 3
    timing_rms = np.sqrt((resid**2).mean(axis=1))

    # ---- Per-track outgoing-velocity test (t0-free, per track) ----
    # For each track, the far-from-vertex hit must arrive later than the near
    # one by exactly (R_far - R_near)/c if the particle leaves the vertex at c.
    def _track_pull(c0, c1):
        R0, R1 = R[:, c0], R[:, c1]
        far0 = R0 >= R1
        t_far = np.where(far0, t_meas[:, c0], t_meas[:, c1])
        t_near = np.where(far0, t_meas[:, c1], t_meas[:, c0])
        dR = np.abs(R0 - R1)
        resid_tr = (t_far - t_near) - dR / c
        return resid_tr / (np.sqrt(2.0) * sigma_t)

    velo_pull_mu = _track_pull(0, 1)   # muon hits  (entry-outer, entry-inner)
    velo_pull_e = _track_pull(2, 3)    # electron hits (exit-inner, exit-outer)
    velo_pull_min = np.minimum(velo_pull_mu, velo_pull_e)

    return dict(sep_in=sep_in, sep_out=sep_out, open_a=open_a, dca=dca,
                collin=collin, vtx_in=vtx_in, pointing=pointing,
                on_tracker=on_tracker, p_soft=p_soft, p_e=p_e, p_mu=p_mu,
                V_reco=V_reco, timing_chi2=timing_chi2, timing_rms=timing_rms,
                velo_pull_min=velo_pull_min, velo_pull_mu=velo_pull_mu,
                velo_pull_e=velo_pull_e,
                # geometry for event displays
                V_true=V, mu_dir=mu_dir, dir_e=dir_e, t_meas=t_meas,
                P_eo=P_eo, P_ei=P_ei, P_xi=P_xi, P_xo=P_xo)


# ------------------------------------------------------------------
# Volume sampler (efficient: every sample is an in-cavern decay)
# ------------------------------------------------------------------

def sample_fiducial_volume(n, seed):
    """Sample n points uniformly inside the fiducial volume by bbox rejection,
    using the fast points_in_fiducial test (vs trimesh.sample.volume_mesh, which
    calls the slow mesh.contains internally)."""
    rng = np.random.default_rng(seed)
    b = mesh_fiducial.bounds
    pts, got = [], 0
    # acceptance ~ V_fid / V_bbox; oversample generously to limit iterations
    accept = max(mesh_fiducial.volume /
                 np.prod(b[1] - b[0]), 0.02)
    while got < n:
        m = int((n - got) / accept * 1.25) + 1000
        cand = np.column_stack([rng.uniform(b[0, i], b[1, i], m) for i in range(3)])
        inside = cand[points_in_fiducial(cand)]
        if len(inside):
            pts.append(inside); got += len(inside)
    return np.vstack(pts)[:n]


def _volume_chunk(m, seed_chunk, rng, theta_max, p_mean, spectrum, i_vertical,
                  sigma_hit, sigma_t, depth_e_min, cmin):
    """One volume-sampling chunk: sample m decays, reconstruct the
    reconstructable ones, and return (reco_dict, sum(1/chord) over geom-valid,
    n_geom). The reco_dict adds inv_lam/cos_z/chord and a scalar n_reco."""
    Vc = sample_fiducial_volume(m, seed_chunk)
    u = rng.uniform(0, 1, m)
    cos_z = ((1 - u) * (1 - cmin**3) + cmin**3) ** (1.0 / 3.0)
    sin_z = np.sqrt(np.maximum(1 - cos_z**2, 0.0))
    phi = rng.uniform(0, 2 * np.pi, m)
    d = np.column_stack([sin_z * np.cos(phi), -cos_z, sin_z * np.sin(phi)])

    if spectrum == 'depth':
        p_mu = sample_depth_momentum(cos_z, theta_max, rng, e_min=depth_e_min)
    elif spectrum == 'empirical':
        p_mu = sample_empirical_momentum(m, rng)
    else:
        p_mu = sample_momentum(m, rng, mean=p_mean, kind=spectrum)
    lam = (p_mu / M_MUON) * CTAU_MUON

    entry = _first_forward_hit(mesh_fiducial, Vc, -d)
    fwd = _first_forward_hit(mesh_fiducial, Vc, d)
    a = np.linalg.norm(entry - Vc, axis=1)
    f_out = np.linalg.norm(fwd - Vc, axis=1)
    chord = a + f_out

    p_e, dir_e = michel_electron(p_mu, d, rng)
    exit_e = _first_forward_hit(mesh_fiducial, Vc, dir_e)
    d_out_e = np.linalg.norm(exit_e - Vc, axis=1)

    geom_ok = (np.isfinite(a) & np.isfinite(f_out) & (a > 1e-3) & (f_out > 1e-3))
    inv_chord_sum = float(np.sum(1.0 / chord[geom_ok]))
    n_geom = int(geom_ok.sum())
    good = (geom_ok & np.isfinite(d_out_e) & (d_out_e > 1e-3)
            & np.isfinite(exit_e[:, 0]))
    ng = int(good.sum())

    rc = reconstruct(Vc[good], d[good], -a[good], np.zeros(ng), dir_e[good],
                     d_out_e[good], p_mu[good], p_e[good], sigma_hit, sigma_t, rng)
    rc['inv_lam'] = 1.0 / lam[good]
    rc['cos_z'] = cos_z[good]
    rc['chord'] = chord[good]
    rc['n_reco'] = ng
    return rc, inv_chord_sum, n_geom


def run_volume(n_events, theta_max, p_mean, spectrum, i_vertical,
               sigma_hit, seed, sigma_t=SIGMA_T_DEFAULT,
               target_muon_rate_hz=MUON_RATE_HZ_DEFAULT, chunk=200_000,
               depth_e_min=0.3, spectrum_weight=1.0):
    """
    Volume importance sampler. Sample the decay vertex uniformly in the
    fiducial volume + a cos^2(theta) downgoing direction + the momentum
    spectrum; every sample is a genuine in-cavern decay (no throw-box waste).

    Normalisation uses track-length density = flux:
        decay rate = V_fid * Integral I(theta)/lambda dOmega
    With directions drawn from the cos^2 intensity pdf the cos^2 cancels, so
    per-sample decay-rate weight = pref / lambda, pref = V_fid*I_v*2pi*Z_g.
    The muon-through rate (for the anchor) is pref*<1/chord> (full muon chord),
    using <chord>_flux = 1/<1/chord>_volume.

    Processed in chunks (default 200k) so arbitrarily large N runs without a
    memory spike; the global normalisation is accumulated across chunks.
    """
    Vfid = mesh_fiducial.volume
    cmin = np.cos(theta_max)
    pref = Vfid * i_vertical * 2 * np.pi * (1 - cmin**3) / 3.0

    rng = np.random.default_rng(seed)
    reco_chunks = []
    inv_chord_sum, n_geom = 0.0, 0
    n_done, ci = 0, 0
    print(f"Volume sampling {n_events:,} in-cavern decays "
          f"(chunks of {chunk:,})...")
    while n_done < n_events:
        m = min(chunk, n_events - n_done)
        rc, ics, ngeo = _volume_chunk(m, seed + 1 + ci, rng, theta_max, p_mean,
                                      spectrum, i_vertical, sigma_hit, sigma_t,
                                      depth_e_min, cmin)
        reco_chunks.append(rc)
        inv_chord_sum += ics
        n_geom += ngeo
        n_done += m
        ci += 1
        print(f"  chunk {ci}: {n_done:,}/{n_events:,} ({rc['n_reco']}/{m} reco)")

    reco = {k: np.concatenate([rc[k] for rc in reco_chunks], axis=0)
            for k in reco_chunks[0] if k != 'n_reco'}
    mean_inv_chord = inv_chord_sum / n_geom
    muon_rate = pref * mean_inv_chord
    # decay-rate per sample [Hz]; spectrum_weight = generated window's spectral
    # fraction when importance-sampling a restricted momentum range (else 1).
    w_all = pref * reco['inv_lam'] / n_events * spectrum_weight
    if target_muon_rate_hz:
        w_all *= target_muon_rate_hz / muon_rate
        muon_rate = target_muon_rate_hz

    reco['w'] = w_all
    reco['nv'] = n_events
    reco['n_reco'] = len(w_all)
    reco['sigma_t_ns'] = sigma_t * 1e9
    reco['muon_rate_hz'] = muon_rate
    return reco


def _passers_to_csv(d, out_csv):
    """Write the collinearity-cut passers (dict of aligned arrays incl. 'w') to
    CSV with the same columns as dump_collin_passers."""
    import pandas as pd
    from grendel_geometry import local_transverse_xy, classify_points_with_basis
    if len(d.get('sep_in', [])) == 0:
        pd.DataFrame().to_csv(out_csv, index=False)
        return 0
    V = d['V_true']
    bx, by = local_transverse_xy(V)
    th, s, *_ = classify_points_with_basis(V)
    gpt = d['sep_in'] < SEP_IN_POINT_GATE
    pass_point = ((~gpt) | (d['pointing'] < POINT_TIGHT_SEP_IN)) \
        & (d['pointing'] < POINT_GLOBAL)
    pass_time = d['timing_chi2'] < CHI2_TIMING_MAX
    pd.DataFrame({
        'x_m': V[:, 0], 'y_m': V[:, 1], 'z_m': V[:, 2],
        'x_local_m': bx, 'y_local_m': by, 's_m': s, 'theta_deg': np.degrees(th),
        'p_mu_GeV': d['p_mu'], 'p_e_GeV': d['p_e'], 'p_soft_GeV': d['p_soft'],
        'sep_in_cm': d['sep_in'] * 100, 'sep_out_cm': d['sep_out'] * 100,
        'collin_mm': d['collin'] * 1000, 'kink_deg': np.degrees(d['open_a']),
        'dca_cm': d['dca'] * 100, 'pointing_mrad': d['pointing'] * 1000,
        'vtx_in': d['vtx_in'], 'timing_chi2': d['timing_chi2'],
        'velo_pull_min': d['velo_pull_min'], 'rate_hz': d['w'],
        'pass_vtx': d['vtx_in'], 'pass_pointing': pass_point,
        'pass_timing': pass_time,
        'pass_all': d['vtx_in'] & pass_point & pass_time,
    }).to_csv(out_csv, index=False)
    return len(V)


def run_volume_streaming(n_events, theta_max, p_mean, spectrum, i_vertical,
                         sigma_hit, seed, livetime_s, collin_csv,
                         sigma_t=SIGMA_T_DEFAULT,
                         target_muon_rate_hz=MUON_RATE_HZ_DEFAULT,
                         chunk=200_000, depth_e_min=0.3, spectrum_weight=1.0):
    """Memory-bounded volume run: accumulate the weighted cutflow per chunk and
    keep only the (rare) collinearity-cut passers, so arbitrarily large N runs
    in fixed memory. Prints the merged cutflow and writes the passers CSV."""
    Vfid = mesh_fiducial.volume
    cmin = np.cos(theta_max)
    pref = Vfid * i_vertical * 2 * np.pi * (1 - cmin**3) / 3.0
    rng = np.random.default_rng(seed)

    names, wsum, nmc = None, None, None
    velo_w, velo_n = 0.0, 0
    inv_chord_sum, n_geom = 0.0, 0
    passers = []
    n_done, ci = 0, 0
    print(f"Volume STREAMING {n_events:,} in-cavern decays "
          f"(chunks of {chunk:,})...")
    while n_done < n_events:
        mm = min(chunk, n_events - n_done)
        rc, ics, ngeo = _volume_chunk(mm, seed + 1 + ci, rng, theta_max, p_mean,
                                      spectrum, i_vertical, sigma_hit, sigma_t,
                                      depth_e_min, cmin)
        inv_chord_sum += ics
        n_geom += ngeo
        w_un = pref * rc['inv_lam'] / n_events * spectrum_weight  # unanchored rate/event
        stages = cutflow_stages(rc)
        if wsum is None:
            names = [s[0] for s in stages]
            wsum = [0.0] * len(stages)
            nmc = [0] * len(stages)
        for k, (_, mask) in enumerate(stages):
            wsum[k] += float(w_un[mask].sum())
            nmc[k] += int(mask.sum())
        pre_t = stages[PRE_TIMING_STAGE][1]
        vmask = pre_t & (rc['velo_pull_min'] > -VELO_NCUT)
        velo_w += float(w_un[vmask].sum())
        velo_n += int(vmask.sum())
        cmask = stages[COLLIN_STAGE][1]
        if cmask.any():
            sub = {k: v[cmask] for k, v in rc.items()
                   if isinstance(v, np.ndarray) and v.shape[0] == len(cmask)}
            sub['w_un'] = w_un[cmask]
            passers.append(sub)
        n_done += mm
        ci += 1
        print(f"  chunk {ci}: {n_done:,}/{n_events:,} ({rc['n_reco']}/{mm} reco)")

    mean_inv_chord = inv_chord_sum / n_geom
    muon_rate = pref * mean_inv_chord
    anchor = (target_muon_rate_hz / muon_rate) if target_muon_rate_hz else 1.0
    if target_muon_rate_hz:
        muon_rate = target_muon_rate_hz
    rates = [x * anchor for x in wsum]

    print(f"\n  Total muon rate through fiducial: {muon_rate:.3e} Hz")
    print(f"  Decay-in-fiducial rate (reco):    {rates[0]:.3e} Hz")
    final, ul = print_cutflow(names, rates, nmc, livetime_s,
                              velo_w * anchor, velo_n)
    if passers:
        merged = {k: np.concatenate([p[k] for p in passers])
                  for k in passers[0]}
        merged['w'] = merged['w_un'] * anchor
        _passers_to_csv(merged, collin_csv)
    else:
        _passers_to_csv({}, collin_csv)
    print(f"  {nmc[COLLIN_STAGE]} events pass the collinearity cut -> {collin_csv}")
    return dict(final=final, ul=ul, muon_rate=muon_rate)


# ------------------------------------------------------------------
# Main MC (throw-plane sampler, kept for cross-checks)
# ------------------------------------------------------------------

def run(n_rays, y_throw, theta_max, p_mean, spectrum, i_vertical,
        sigma_hit, seed, sigma_t=SIGMA_T_DEFAULT,
        target_muon_rate_hz=MUON_RATE_HZ_DEFAULT):
    rng = np.random.default_rng(seed)

    origins, dirs, cos_z, area = sample_cosmics(n_rays, y_throw, theta_max, rng)

    print(f"Ray-casting {n_rays:,} cosmics against the fiducial mesh...")
    hit, d_in, d_out = traverse_fiducial(origins, dirs)
    nv = int(hit.sum())
    print(f"  {nv:,} / {n_rays:,} traverse the fiducial "
          f"({nv / n_rays * 100:.4f}% of throw box).")
    if nv == 0:
        return None

    mu_ori = origins[hit]
    mu_dir = dirs[hit]
    d_in = d_in[hit]
    d_out = d_out[hit]
    cos_z = cos_z[hit]
    chord = d_out - d_in

    # Momentum: realistic depth spectrum (per-muon, angle-dependent) or a toy.
    if spectrum == 'depth':
        p_mu = sample_depth_momentum(cos_z, theta_max, rng)
    elif spectrum == 'empirical':
        p_mu = sample_empirical_momentum(nv, rng)
    else:
        p_mu = sample_momentum(nv, rng, mean=p_mean, kind=spectrum)

    # ---- Flux normalisation weight per traversing muon (Hz) ----
    # Rate through plane = area * I_v * 2pi * Z_g * <cos theta * traverse>
    # with Z_g = (1 - cos^3 theta_max) / 3 and directions sampled ~ g.
    cos_min = np.cos(theta_max)
    Z_g = (1 - cos_min**3) / 3.0
    w_flux = area * i_vertical * 2 * np.pi * Z_g * cos_z / n_rays  # Hz/muon
    # Anchor the absolute scale to the user-supplied total tunnel rate (the
    # cos(theta) weighting between muons is preserved by the rescaling).
    if target_muon_rate_hz and w_flux.sum() > 0:
        w_flux *= target_muon_rate_hz / w_flux.sum()

    # ---- In-flight decay inside the chord ----
    bg = p_mu / M_MUON
    lam = bg * CTAU_MUON                       # lab decay length (m)
    w_decay = 1.0 - np.exp(-chord / lam)       # P(decay within chord)
    u = rng.uniform(0, 1, nv)
    s_rel = -lam * np.log1p(-u * w_decay)      # decay offset from entry
    d_decay = d_in + s_rel

    # ---- Michel electron ----
    p_e, dir_e = michel_electron(p_mu, mu_dir, rng)

    # ---- Electron exit through the fiducial wall ----
    V = mu_ori + d_decay[:, None] * mu_dir
    exit_e = _first_forward_hit(mesh_fiducial, V, dir_e)
    good = ~np.isnan(exit_e[:, 0])
    d_out_e = np.full(nv, np.nan)
    d_out_e[good] = np.linalg.norm(exit_e[good] - V[good], axis=1)

    # keep electrons that exit the fiducial
    sel = good & np.isfinite(d_out_e) & (d_out_e > 1e-3)

    reco = reconstruct(mu_ori[sel], mu_dir[sel], d_in[sel], d_decay[sel],
                       dir_e[sel], d_out_e[sel], p_mu[sel], p_e[sel],
                       sigma_hit, sigma_t, rng)
    reco['w'] = w_flux[sel] * w_decay[sel]     # Hz per reconstructed event
    reco['w_flux'] = w_flux[sel]
    reco['cos_z'] = cos_z[sel]
    reco['chord'] = chord[sel]
    reco['p_mu'] = p_mu[sel]
    reco['nv'] = nv
    reco['n_reco'] = int(sel.sum())
    reco['sigma_t_ns'] = sigma_t * 1e9
    reco['muon_rate_hz'] = float(w_flux.sum())   # all traversing muons
    return reco


# ------------------------------------------------------------------
# Cutflow (full signal selection)
# ------------------------------------------------------------------

# Cut sequence shared by the in-memory printer (weighted_cutflow) and the
# streaming accumulator (run_volume_streaming), so they are guaranteed identical.
# Index 7 = collinearity stage, 10 = pre-timing (everything through the global
# pointing cut), 11 = final (timing).
COLLIN_STAGE, PRE_TIMING_STAGE = 7, 10


def cutflow_stages(r):
    """Ordered list of (name, cumulative_mask) for the full selection."""
    sep, sepo = r['sep_in'], r['sep_out']
    L = DETECTOR_THICKNESS
    out = []
    m = np.ones(len(sep), dtype=bool)
    out.append(('decay in fiducial (reco)', m.copy()))
    m = m & r['on_tracker']
    out.append(('both tracks on tracker wall', m.copy()))
    m = m & (r['p_soft'] >= P_SOFT_CUT)
    out.append((f'p_soft > {P_SOFT_CUT*1000:.0f} MeV/c', m.copy()))
    m = m & (sep >= SEP_MIN) & (sepo >= SEP_MIN)
    out.append((f'sep_in & sep_out > {SEP_MIN*1000:.0f} mm', m.copy()))
    m = m & (sep <= SEP_MAX)
    out.append((f'sep_in < {SEP_MAX:.0f} m', m.copy()))
    m = m & (r['dca'] <= DCA_CUT)
    out.append((f'DCA < {DCA_CUT*100:.0f} cm', m.copy()))
    par = r['open_a'] < THETA_PARALLEL
    m = m & ((~par) | (sepo < SEP_OUT_MAX_PARALLEL))
    out.append(('conditional max sep_outer (parallel)', m.copy()))
    gated = sepo > SEP_OUT_GATE
    m = m & ((~gated) | (r['collin'] > COLLIN_FRAC * L))
    out.append((f'collin>{COLLIN_FRAC:.2f}L if sep_out>{SEP_OUT_GATE*100:.0f}cm '
                f'({COLLIN_FRAC*L*1000:.0f}mm)', m.copy()))
    m = m & r['vtx_in']
    out.append(('vertex in fiducial (PCA)', m.copy()))
    gpt = sep < SEP_IN_POINT_GATE
    m = m & ((~gpt) | (r['pointing'] < POINT_TIGHT_SEP_IN))
    out.append(('conditional tight pointing', m.copy()))
    m = m & (r['pointing'] < POINT_GLOBAL)
    out.append((f'global pointing < {POINT_GLOBAL*1000:.0f} mrad', m.copy()))
    m = m & (r['timing_chi2'] < CHI2_TIMING_MAX)
    out.append((f'timing chi2 < {CHI2_TIMING_MAX:.0f}', m.copy()))
    return out


def print_cutflow(names, rates, nmcs, livetime_s, velo_rate, velo_nmc):
    """Print the cutflow table from per-stage rates [Hz] and MC counts; returns
    (final_rate, ul_or_None). Used by both the in-memory and streaming paths."""
    print(f"\nWeighted cutflow (rate in Hz; livetime = "
          f"{livetime_s:.3e} s -> events):")
    print(f"  {'cut':<40} {'rate [Hz]':>11} {'events':>11} {'N_mc':>7} {'rel':>7}")
    print("  " + "-" * 80)
    prev = rates[0] if rates else 0.0
    last_w, last_n = None, None
    for name, rate, nmc in zip(names, rates, nmcs):
        rel = rate / prev if prev > 0 else 0.0
        print(f"  {name:<40} {rate:>11.3e} {rate*livetime_s:>11.3e} "
              f"{nmc:>7d} {rel:>7.4f}")
        if nmc > 0:
            last_w, last_n = rate, nmc
        prev = rate
    final, n_final = rates[-1], nmcs[-1]
    print("  " + "-" * 80)
    print(f"  {'FINAL':<40} {final:>11.3e} {final*livetime_s:>11.3e} {n_final:>7d}")
    print(f"  [compare] per-track velocity (pull>-{VELO_NCUT:.2f}) on same "
          f"sample : {velo_rate:.3e} Hz  ({velo_nmc} MC)")
    ul = None
    if n_final == 0 and last_n:
        ul = 2.3 * last_w / last_n
        print(f"\n  No MC event survives -> 90% CL upper limit:")
        print(f"    < {ul:.2e} Hz  =  < {ul*livetime_s:.2e} events "
              f"(2.3 x mean weight of the {last_n} pre-veto events)")
    return final, ul


def weighted_cutflow(r, livetime_s):
    w = r['w']
    stages = cutflow_stages(r)
    names = [nm for nm, _ in stages]
    rates = [float(w[mask].sum()) for _, mask in stages]
    nmcs = [int(mask.sum()) for _, mask in stages]
    pre_timing = stages[PRE_TIMING_STAGE][1]
    vmask = pre_timing & (r['velo_pull_min'] > -VELO_NCUT)
    print_cutflow(names, rates, nmcs, livetime_s,
                  float(w[vmask].sum()), int(vmask.sum()))
    return rates[-1], stages[-1][1], pre_timing


# ------------------------------------------------------------------
# Plots
# ------------------------------------------------------------------

def make_plot(r, out_path, interactive):
    import matplotlib.pyplot as plt
    w = r['w']
    fig, ax = plt.subplots(2, 4, figsize=(24, 10))

    a = ax[0, 0]
    a.hist(r['p_e'], bins=np.linspace(0, 2, 60), weights=w,
           color='steelblue', edgecolor='k', linewidth=0.3)
    a.axvline(P_SOFT_CUT, color='red', ls='--',
              label=f'p_cut = {P_SOFT_CUT*1e3:.0f} MeV')
    a.set_xlabel('electron lab momentum (GeV/c)')
    a.set_ylabel('rate [Hz]'); a.legend(); a.set_title('Michel electron momentum')

    a = ax[0, 1]
    a.hist(np.degrees(r['open_a']), bins=60, weights=w,
           color='darkgreen', edgecolor='k', linewidth=0.3)
    a.set_xlabel('reco opening (kink) angle (deg)')
    a.set_ylabel('rate [Hz]'); a.set_title('muon-electron kink')

    a = ax[0, 2]
    a.hist(r['collin'] * 1000, bins=np.linspace(0, 50, 60), weights=w,
           color='purple', edgecolor='k', linewidth=0.3)
    a.axvline(COLLIN_FRAC * DETECTOR_THICKNESS * 1000, color='red', ls='--',
              label=f'cut = {COLLIN_FRAC:.2f}L = '
                    f'{COLLIN_FRAC*DETECTOR_THICKNESS*1e3:.0f} mm')
    a.set_xlabel('collinearity (mm)'); a.set_ylabel('rate [Hz]')
    a.legend(); a.set_title('4-hit collinearity')

    a = ax[1, 0]
    a.hist(r['dca'] * 100, bins=np.linspace(0, 30, 60), weights=w,
           color='orange', edgecolor='k', linewidth=0.3)
    a.axvline(DCA_CUT * 100, color='red', ls='--',
              label=f'cut = {DCA_CUT*100:.0f} cm')
    a.set_xlabel('DCA (cm)'); a.set_ylabel('rate [Hz]')
    a.legend(); a.set_title('track DCA')

    a = ax[1, 1]
    a.hist(r['sep_in'], bins=60, weights=w, color='teal',
           edgecolor='k', linewidth=0.3)
    a.axvline(SEP_MAX, color='red', ls='--', label=f'max = {SEP_MAX:.0f} m')
    a.set_xlabel('sep_inner (m)'); a.set_ylabel('rate [Hz]')
    a.legend(); a.set_title('inner-layer separation')

    a = ax[1, 2]
    a.hist2d(np.degrees(np.arccos(np.clip(r['cos_z'], -1, 1))), r['chord'],
             bins=[np.linspace(0, 90, 50), np.linspace(0, 8, 50)],
             weights=w, cmap='magma')
    a.set_xlabel('zenith (deg)'); a.set_ylabel('fiducial chord (m)')
    a.set_title('geometry (rate-weighted)')

    a = ax[0, 3]
    a.hist(np.clip(r['timing_chi2'], 1e-2, 1e6), bins=np.logspace(-2, 6, 60),
           weights=w, color='crimson', edgecolor='k', linewidth=0.3)
    a.axvline(CHI2_TIMING_MAX, color='black', ls='--',
              label=f'cut = {CHI2_TIMING_MAX:.0f}')
    a.set_xscale('log'); a.set_yscale('log')
    a.set_xlabel('timing chi2 (ndof=3)'); a.set_ylabel('rate [Hz]')
    a.legend(); a.set_title(f'timing consistency (sigma_t={r["sigma_t_ns"]:.1f} ns)')

    a = ax[1, 3]
    a.hist(np.clip(r['timing_rms'] * 1e9, 0, 60), bins=60, weights=w,
           color='slateblue', edgecolor='k', linewidth=0.3)
    a.axvline(r['sigma_t_ns'], color='red', ls='--',
              label=f'sigma_t = {r["sigma_t_ns"]:.1f} ns')
    a.set_xlabel('timing residual RMS (ns)'); a.set_ylabel('rate [Hz]')
    a.legend(); a.set_title('per-event timing residual')

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    print(f"\nSaved {out_path}")
    if interactive:
        plt.show()
    else:
        plt.close('all')


def make_variable_plots(r, out_dir, interactive=False):
    """One rate-weighted distribution per file (signal-style) into out_dir, with
    the ACTIVE cut thresholds drawn (single-sourced). Mirrors the signal
    per-mass plot set: Michel momentum, kink, collinearity, DCA, sep_in,
    sep_out, pointing, timing chi2, timing residual, geometry."""
    import matplotlib.pyplot as plt
    os.makedirs(out_dir, exist_ok=True)
    w = r['w']
    L = DETECTOR_THICKNESS
    collin_cut_mm = COLLIN_FRAC * L * 1000

    def save(name):
        p = os.path.join(out_dir, name)
        plt.tight_layout(); plt.savefig(p, dpi=150)
        if interactive: plt.show()
        else: plt.close('all')

    # 1. Michel electron momentum
    plt.figure(figsize=(7, 5))
    plt.hist(r['p_e'], bins=np.linspace(0, 2, 60), weights=w,
             color='steelblue', edgecolor='k', linewidth=0.3)
    plt.axvline(P_SOFT_CUT, color='red', ls='--', label=f'p_cut = {P_SOFT_CUT*1e3:.0f} MeV')
    plt.xlabel('electron lab momentum (GeV/c)'); plt.ylabel('rate [Hz]')
    plt.legend(); plt.title('Michel electron momentum'); save('michel_momentum.png')

    # 2. kink angle
    plt.figure(figsize=(7, 5))
    plt.hist(np.degrees(r['open_a']), bins=60, weights=w,
             color='darkgreen', edgecolor='k', linewidth=0.3)
    plt.xlabel('reco opening (kink) angle (deg)'); plt.ylabel('rate [Hz]')
    plt.title('muon-electron kink'); save('kink_angle.png')

    # 3. collinearity (ACTIVE cut + L/2 ceiling)
    plt.figure(figsize=(7, 5))
    plt.hist(r['collin'] * 1000, bins=np.linspace(0, max(collin_cut_mm*1.3, 150), 80),
             weights=w, color='purple', edgecolor='k', linewidth=0.3)
    plt.axvline(collin_cut_mm, color='red', ls='--',
                label=f'cut = {COLLIN_FRAC:.2f}L = {collin_cut_mm:.0f} mm')
    plt.axvline(L/2*1000, color='green', ls=':',
                label=f'signal ceiling L/2 = {L/2*1000:.0f} mm')
    plt.yscale('log')
    plt.xlabel('collinearity (mm)'); plt.ylabel('rate [Hz]')
    plt.legend(); plt.title(f'4-hit collinearity (gate sep_out>{SEP_OUT_GATE*100:.0f}cm)')
    save('collinearity.png')

    # 4. DCA
    plt.figure(figsize=(7, 5))
    plt.hist(r['dca'] * 100, bins=np.linspace(0, 30, 60), weights=w,
             color='orange', edgecolor='k', linewidth=0.3)
    plt.axvline(DCA_CUT * 100, color='red', ls='--', label=f'cut = {DCA_CUT*100:.0f} cm')
    plt.xlabel('DCA (cm)'); plt.ylabel('rate [Hz]')
    plt.legend(); plt.title('track DCA'); save('dca.png')

    # 5. sep_inner
    plt.figure(figsize=(7, 5))
    plt.hist(r['sep_in'] * 100, bins=np.linspace(0, 100, 60), weights=w,
             color='teal', edgecolor='k', linewidth=0.3)
    plt.xlabel('inner-layer separation (cm)'); plt.ylabel('rate [Hz]')
    plt.title('sep_inner'); save('sep_inner.png')

    # 6. sep_outer (ACTIVE gate)
    plt.figure(figsize=(7, 5))
    plt.hist(r['sep_out'] * 100, bins=np.linspace(0, 100, 60), weights=w,
             color='darkcyan', edgecolor='k', linewidth=0.3)
    plt.axvline(SEP_OUT_GATE * 100, color='red', ls='--',
                label=f'collinearity gate = {SEP_OUT_GATE*100:.0f} cm')
    plt.xlabel('outer-layer separation (cm)'); plt.ylabel('rate [Hz]')
    plt.legend(); plt.title('sep_outer'); save('sep_outer.png')

    # 7. pointing (tight + global cut lines, log-x)
    plt.figure(figsize=(7, 5))
    pt = np.clip(r['pointing'] * 1000, 1e-1, 3000)
    plt.hist(pt, bins=np.logspace(-1, np.log10(3000), 60), weights=w,
             color='darkorange', edgecolor='k', linewidth=0.3)
    plt.axvline(POINT_GLOBAL * 1000, color='purple', ls='-',
                label=f'global cut = {POINT_GLOBAL*1000:.0f} mrad')
    plt.axvline(POINT_TIGHT_SEP_IN * 1000, color='blue', ls=':',
                label=f'tight cut (sep_in<{SEP_IN_POINT_GATE*100:.0f}cm) = {POINT_TIGHT_SEP_IN*1000:.0f} mrad')
    plt.xscale('log'); plt.yscale('log')
    plt.xlabel('pointing angle (mrad)'); plt.ylabel('rate [Hz]')
    plt.legend(); plt.title('pointing (bisector vs IP flight)'); save('pointing.png')

    # 8. timing chi2
    plt.figure(figsize=(7, 5))
    plt.hist(np.clip(r['timing_chi2'], 1e-2, 1e6), bins=np.logspace(-2, 6, 60),
             weights=w, color='crimson', edgecolor='k', linewidth=0.3)
    plt.axvline(CHI2_TIMING_MAX, color='black', ls='--', label=f'cut = {CHI2_TIMING_MAX:.0f}')
    plt.xscale('log'); plt.yscale('log')
    plt.xlabel('timing chi2 (ndof=3)'); plt.ylabel('rate [Hz]')
    plt.legend(); plt.title(f'timing consistency (sigma_t={r["sigma_t_ns"]:.1f} ns)')
    save('timing_chi2.png')

    # 9. timing residual
    plt.figure(figsize=(7, 5))
    plt.hist(np.clip(r['timing_rms'] * 1e9, 0, 60), bins=60, weights=w,
             color='slateblue', edgecolor='k', linewidth=0.3)
    plt.axvline(r['sigma_t_ns'], color='red', ls='--', label=f'sigma_t = {r["sigma_t_ns"]:.1f} ns')
    plt.xlabel('timing residual RMS (ns)'); plt.ylabel('rate [Hz]')
    plt.legend(); plt.title('per-event timing residual'); save('timing_residual.png')

    # 10. geometry (zenith vs chord)
    plt.figure(figsize=(7, 5))
    plt.hist2d(np.degrees(np.arccos(np.clip(r['cos_z'], -1, 1))), r['chord'],
               bins=[np.linspace(0, 90, 50), np.linspace(0, 8, 50)],
               weights=w, cmap='magma')
    plt.colorbar(label='rate [Hz]')
    plt.xlabel('zenith (deg)'); plt.ylabel('fiducial chord (m)')
    plt.title('geometry (rate-weighted)'); save('geometry.png')

    print(f"Saved per-variable plots -> {out_dir}/")


# ------------------------------------------------------------------

def local_xy(points):
    """Project 3D points into the tunnel cross-section local frame (x_local,
    y_local) used by the horseshoe diagram: offset from the nearest centreline
    point onto (right, up). Vectorised over the centreline segments."""
    from grendel_geometry import path_3d_fiducial as P
    pts = np.asarray(points, dtype=float)
    m = len(pts)
    best = np.full(m, np.inf)
    bx, by = np.zeros(m), np.zeros(m)
    for i in range(len(P) - 1):
        seg = P[i + 1] - P[i]
        L = np.linalg.norm(seg)
        if L == 0:
            continue
        sh = seg / L
        wu = np.array([0., 1., 0.]) if abs(sh[1]) < 0.9 else np.array([0., 0., 1.])
        right = np.cross(sh, wu); right /= np.linalg.norm(right)
        up = np.cross(right, sh); up /= np.linalg.norm(up)
        rel = pts - P[i]
        t = np.clip(rel @ sh, 0, L)
        diff = pts - (P[i] + np.outer(t, sh))
        d2 = np.einsum('ij,ij->i', diff, diff)
        upd = d2 < best
        best[upd] = d2[upd]
        bx[upd] = diff[upd] @ right
        by[upd] = diff[upd] @ up
    return bx, by


def _draw_horseshoe(ax):
    """Outline the tunnel cross-section: outer wall coloured by surface role
    (tracker green / scintillator red) + dashed inner fiducial boundary."""
    from grendel_geometry import tunnel_profile_points, DETECTOR_THICKNESS
    outer = tunnel_profile_points(inset=0.0)
    inner = tunnel_profile_points(inset=DETECTOR_THICKNESS)
    inner = np.vstack([inner, inner[:1]])
    N_WALL, N_ARCH = 4, 32
    i_fl, i_fr = 0, 1
    i_rwt = 1 + N_WALL
    i_lwt = i_rwt + (N_ARCH - 1)
    i_lwe = i_lwt + N_WALL
    TRK, SCI = '#2ca02c', '#d62728'
    ax.plot(outer[[i_fl, i_fr], 0], outer[[i_fl, i_fr], 1], color=SCI, lw=3)
    ax.plot(outer[i_fr:i_rwt + 1, 0], outer[i_fr:i_rwt + 1, 1], color=SCI, lw=3)
    ax.plot(outer[i_rwt:i_lwt + 1, 0], outer[i_rwt:i_lwt + 1, 1], color=TRK, lw=3)
    lw_seg = np.vstack([outer[i_lwt:i_lwe + 1], outer[[i_fl]]])
    ax.plot(lw_seg[:, 0], lw_seg[:, 1], color=TRK, lw=3)
    ax.plot(inner[:, 0], inner[:, 1], color='k', lw=1.0, ls='--')


def pre_timing_mask(r, apply_collin=True):
    """Full selection up to (not including) the timing cut. With
    apply_collin=False the collinearity veto is skipped (to show its effect)."""
    m = r['on_tracker'].copy()
    m &= (r['p_soft'] >= P_SOFT_CUT)
    m &= (r['sep_in'] >= SEP_MIN) & (r['sep_out'] >= SEP_MIN)
    m &= (r['sep_in'] <= SEP_MAX)
    m &= (r['dca'] <= DCA_CUT)
    par = r['open_a'] < THETA_PARALLEL
    m &= (~par) | (r['sep_out'] < SEP_OUT_MAX_PARALLEL)
    if apply_collin:
        gated = r['sep_out'] > SEP_OUT_GATE
        m &= (~gated) | (r['collin'] > COLLIN_FRAC * DETECTOR_THICKNESS)
    m &= r['vtx_in']
    gpt = r['sep_in'] < SEP_IN_POINT_GATE
    m &= (~gpt) | (r['pointing'] < POINT_TIGHT_SEP_IN)
    m &= r['pointing'] < POINT_GLOBAL
    return m


def make_collin_sepout_plot(r, out_path, label='', interactive=False):
    """Collinearity vs sep_outer for the background, with all pre-timing cuts
    EXCEPT collinearity applied, so the conditional collinearity-veto region is
    visible. Overlays the sep_out gate and the 30 mm / L-scaled cut lines."""
    import matplotlib.pyplot as plt
    mask = pre_timing_mask(r, apply_collin=False)
    if not mask.any():
        print("[collin-sepout] no events to draw")
        return
    so_cm = r['sep_out'][mask] * 100.0
    co_mm = r['collin'][mask] * 1000.0
    w = r['w'][mask]

    fig, ax = plt.subplots(figsize=(9, 7))
    h = ax.hist2d(so_cm, co_mm,
                  bins=[np.logspace(0, np.log10(max(so_cm.max(), 1e3)), 55),
                        np.logspace(-1, 3, 55)],
                  weights=w, cmap='viridis', cmin=1e-12)
    plt.colorbar(h[3], ax=ax, label='rate [Hz] per bin')
    ax.set_ylim(1e-1, 1e3)

    # Active cut: when sep_out > SEP_OUT_GATE require collin > COLLIN_FRAC*L.
    gate_cm = SEP_OUT_GATE * 100.0
    tgt_mm = COLLIN_FRAC * DETECTOR_THICKNESS * 1000.0
    ax.axvline(gate_cm, color='red', ls='--', lw=1.8,
               label=f'gate = {gate_cm:.0f} cm')
    ax.axhline(tgt_mm, color='red', ls='-', lw=1.8,
               label=f'cut = {COLLIN_FRAC:.2f}L = {tgt_mm:.0f} mm')
    ax.axhline(DETECTOR_THICKNESS / 2 * 1000, color='green', ls=':', lw=1.5,
               label=f'signal ceiling L/2 = {DETECTOR_THICKNESS/2*1000:.0f} mm')
    ax.fill_between([gate_cm, 1e4], 1e-1, tgt_mm,
                    color='red', alpha=0.15, label='reject (muon-like)')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel('sep_outer (cm)'); ax.set_ylabel('collinearity (mm)')
    ax.set_title(f'Cosmic muon-decay: collinearity vs sep_outer\n'
                 f'{label}  (pre-timing, collinearity cut NOT applied; '
                 f'{int(mask.sum())} MC)')
    ax.legend(fontsize=8, loc='lower right')
    ax.grid(True, which='both', alpha=0.25)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    print(f"Saved {out_path}")
    plt.show() if interactive else plt.close('all')


def make_collin_ceiling_plot(r, out_path, interactive=False):
    """Cosmic-decay collinearity distribution (rate-weighted, log y), for events
    passing all cuts except the collinearity veto. Two panels: linear-x and
    log-x in collinearity. Cut (collin > COLLIN_FRAC*L) and the L/2 signal
    ceiling are marked for reference."""
    import matplotlib.pyplot as plt
    L = DETECTOR_THICKNESS
    cut_mm = COLLIN_FRAC * L * 1000.0
    ceil_mm = L / 2 * 1000.0
    mask = pre_timing_mask(r, apply_collin=False)
    if not mask.any():
        print("[collinearity] no events"); return
    co = r['collin'][mask] * 1000.0          # mm
    w = r['w'][mask]

    fig, ax = plt.subplots(1, 2, figsize=(15, 5.5))
    for col, logx in ((0, False), (1, True)):
        a = ax[col]
        if logx:
            pos = co > 0
            bins = np.logspace(np.log10(max(co[pos].min(), 1e-2)),
                               np.log10(co.max() * 1.05), 60)
            a.hist(co[pos], bins=bins, weights=w[pos], color='steelblue',
                   edgecolor='k', linewidth=0.3)
            a.set_xscale('log')
        else:
            a.hist(co, bins=60, weights=w, color='steelblue',
                   edgecolor='k', linewidth=0.3)
        a.axvline(cut_mm, color='red', ls='--', lw=1.8,
                  label=f'cut collin>{COLLIN_FRAC:.2f}L = {cut_mm:.0f} mm')
        a.axvline(ceil_mm, color='green', ls=':', lw=1.5,
                  label=f'signal ceiling L/2 = {ceil_mm:.0f} mm')
        a.set_yscale('log')
        a.set_xlabel('collinearity (mm)')
        a.set_ylabel('rate [Hz]')
        a.set_title(('log-x' if logx else 'linear-x'))
        a.legend(fontsize=8)
        a.grid(True, alpha=0.3, which='both')
    fig.suptitle(f'Cosmic muon-decay collinearity '
                 f'(L = {L*100:.0f} cm; {int(mask.sum())} MC, pre-collinearity)',
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    plt.savefig(out_path, dpi=150)
    print(f"Saved {out_path}")
    plt.show() if interactive else plt.close('all')


def make_decay_horseshoe(r, mask, out_path, label='', interactive=False):
    """Rate-weighted heatmap of decay-vertex positions in the tunnel
    cross-section (the 'horseshoe'), for the events selected by ``mask``."""
    import matplotlib.pyplot as plt
    from grendel_geometry import tunnel_profile_points
    if not mask.any():
        print("[horseshoe] no events to draw")
        return
    V = r['V_true'][mask]
    w = r['w'][mask]
    bx, by = local_xy(V)

    outer = tunnel_profile_points(inset=0.0)
    xlo, xhi = outer[:, 0].min() - 0.15, outer[:, 0].max() + 0.15
    ylo, yhi = outer[:, 1].min() - 0.15, outer[:, 1].max() + 0.15

    fig, ax = plt.subplots(figsize=(9, 8))
    h = ax.hist2d(bx, by, bins=[np.linspace(xlo, xhi, 40),
                                np.linspace(ylo, yhi, 40)],
                  weights=w, cmap='viridis', cmin=1e-12)
    plt.colorbar(h[3], ax=ax, label='decay rate [Hz] per bin')
    _draw_horseshoe(ax)
    ax.plot(0, 0, '+', color='white', ms=10, mew=1.5)
    ax.set_aspect('equal')
    ax.set_xlim(xlo, xhi); ax.set_ylim(ylo, yhi)
    ax.set_xlabel('x_local (m)'); ax.set_ylabel('y_local (m)')
    ax.set_title(f'Cosmic muon-decay vertices in the cavern cross-section\n'
                 f'{label}  (total {w.sum():.2e} Hz, {int(mask.sum())} MC events)\n'
                 f'green = tracker, red = scintillator/veto')
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    print(f"Saved {out_path}")
    plt.show() if interactive else plt.close('all')


def make_cosmic_displays(r, mask, n_show, out_dir, seed=1):
    """
    Event displays (3D + top-down + side + timing) for cosmic muon-decay
    events passing the selection. Mirrors the signal event_display layout but
    shows the incoming muon stub, the decay vertex, the outgoing electron, the
    4 hits, and a timing panel (measured hit time vs distance from the reco
    vertex, against the 'leaves vertex at c' line).
    """
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
    from event_display import _draw_cavern, _project_3d, _scatter_3d

    idxs = np.where(mask)[0]
    if len(idxs) == 0:
        print("[cosmic_display] no passing events to draw")
        return
    rng = np.random.default_rng(seed)
    pick = rng.choice(idxs, size=min(n_show, len(idxs)), replace=False)
    pick.sort()
    os.makedirs(out_dir, exist_ok=True)
    c = SPEED_OF_LIGHT

    for k, i in enumerate(pick):
        V, Vr = r['V_true'][i], r['V_reco'][i]
        P = [r['P_eo'][i], r['P_ei'][i], r['P_xi'][i], r['P_xo'][i]]
        pts = np.array([V, Vr] + P)
        # Adaptive padding: passing events are compact, so scale the view to
        # the event size (with a floor) rather than a fixed large pad.
        spread = float(max(np.ptp(pts[:, 0]), np.ptp(pts[:, 1]),
                           np.ptp(pts[:, 2])))
        pad = max(0.4, 0.5 * spread)
        xlim = (pts[:, 0].min() - pad, pts[:, 0].max() + pad)
        ylim = (pts[:, 1].min() - pad, pts[:, 1].max() + pad)
        zlim = (pts[:, 2].min() - pad, pts[:, 2].max() + pad)

        fig = plt.figure(figsize=(18, 11))
        gs = fig.add_gridspec(2, 3, height_ratios=[1.1, 1.0])
        ax3d = fig.add_subplot(gs[0, 0], projection='3d')
        ax_xz = fig.add_subplot(gs[0, 1])
        ax_yz = fig.add_subplot(gs[0, 2])
        ax_t = fig.add_subplot(gs[1, :])

        _draw_cavern(ax3d, ax_xz, ax_yz, (xlim, ylim, zlim))

        muon = np.vstack([P[0], P[1], V])   # entry-outer -> entry-inner -> vertex
        elec = np.vstack([V, P[2], P[3]])   # vertex -> exit-inner -> exit-outer
        for seg, col, lbl in ((muon, 'tab:orange', 'muon (incoming)'),
                              (elec, 'tab:green', 'electron (outgoing)')):
            _project_3d(ax3d, seg, color=col, lw=1.5)
            ax_xz.plot(seg[:, 2], seg[:, 0], color=col, lw=1.5, label=lbl)
            ax_yz.plot(seg[:, 2], seg[:, 1], color=col, lw=1.5)

        for p in P:
            _scatter_3d(ax3d, p, color='black', marker='s', s=28, zorder=9)
            ax_xz.scatter(p[2], p[0], color='black', marker='s', s=28, zorder=9)
            ax_yz.scatter(p[2], p[1], color='black', marker='s', s=28, zorder=9)
        _scatter_3d(ax3d, V, color='red', marker='*', s=140, zorder=10)
        ax_xz.scatter(V[2], V[0], color='red', marker='*', s=140, zorder=10,
                      label='true vertex')
        ax_yz.scatter(V[2], V[1], color='red', marker='*', s=140, zorder=10)
        ax_xz.scatter(Vr[2], Vr[0], color='blue', marker='x', s=70, zorder=10,
                      label='reco vertex')
        ax_yz.scatter(Vr[2], Vr[1], color='blue', marker='x', s=70, zorder=10)

        ax3d.set_xlim(*xlim); ax3d.set_ylim(*zlim); ax3d.set_zlim(*ylim)
        ax3d.set_xlabel('x (m)'); ax3d.set_ylabel('z (m)'); ax3d.set_zlabel('y (m)')
        ax3d.set_title('3D (zoom)')
        ax_xz.set_xlim(*zlim); ax_xz.set_ylim(*xlim)
        ax_xz.set_xlabel('z (m)'); ax_xz.set_ylabel('x (m)')
        ax_xz.grid(True, alpha=0.3); ax_xz.set_title('Top-down (x vs z)')
        ax_xz.legend(fontsize=8, loc='best')
        ax_yz.set_xlim(*zlim); ax_yz.set_ylim(*ylim)
        ax_yz.set_xlabel('z (m)'); ax_yz.set_ylabel('y (m)')
        ax_yz.grid(True, alpha=0.3); ax_yz.set_title('Side (y vs z)')

        # ---- timing panel ----
        Rd = np.array([np.linalg.norm(p - Vr) for p in P])     # m from reco vtx
        t_ns = r['t_meas'][i] * 1e9
        t0 = np.mean(r['t_meas'][i] - Rd / c)
        lbls = ['muon outer', 'muon inner', 'e inner', 'e outer']
        cols = ['tab:orange', 'tab:orange', 'tab:green', 'tab:green']
        mk = ['o', 's', 's', 'o']
        for j in range(4):
            ax_t.scatter(Rd[j], t_ns[j], color=cols[j], marker=mk[j], s=80,
                         zorder=6, label=lbls[j])
        Rline = np.linspace(0, Rd.max() * 1.1 + 0.1, 50)
        ax_t.plot(Rline, (t0 + Rline / c) * 1e9, 'k--',
                  label='leaves vertex at c (best t0)')
        ax_t.fill_between(Rline, (t0 + Rline / c) * 1e9 - r['sigma_t_ns'],
                          (t0 + Rline / c) * 1e9 + r['sigma_t_ns'],
                          color='gray', alpha=0.2, label=f'±σ_t')
        ax_t.set_xlabel('distance from reco vertex R (m)')
        ax_t.set_ylabel('measured hit time (ns)')
        ax_t.set_title(f'Timing: χ²={r["timing_chi2"][i]:.1f} (cut<{CHI2_TIMING_MAX:.0f}), '
                       f'RMS={r["timing_rms"][i]*1e9:.2f} ns')
        ax_t.legend(fontsize=8, ncol=2); ax_t.grid(True, alpha=0.3)

        zen = np.degrees(np.arccos(np.clip(r['cos_z'][i], -1, 1)))
        fig.suptitle(
            f"Cosmic μ-decay event ({k+1}/{len(pick)})   "
            f"p_μ={r['p_mu'][i]:.2f}, p_e={r['p_e'][i]:.2f} GeV/c   "
            f"zenith={zen:.0f}°   kink={np.degrees(r['open_a'][i]):.0f}°   "
            f"sep_in={r['sep_in'][i]:.2f} m   DCA={r['dca'][i]*100:.1f} cm   "
            f"rate={r['w'][i]:.2e} Hz",
            fontsize=10, x=0.02, ha='left')
        fig.tight_layout(rect=(0, 0, 1, 0.95))
        fname = os.path.join(out_dir, f"cosmic_decay_event_{k:02d}.png")
        fig.savefig(fname, dpi=120)
        plt.close(fig)
        print(f"  saved {fname}")


def dump_collin_passers(r, out_csv):
    """Write a CSV with the full details of every event passing the cumulative
    selection through the collinearity cut (the candidates that then face the
    vertex / pointing / timing cuts)."""
    import pandas as pd
    from grendel_geometry import local_transverse_xy, classify_points_with_basis
    L = DETECTOR_THICKNESS
    sep, sepo = r['sep_in'], r['sep_out']
    m = r['on_tracker'].copy()
    m &= r['p_soft'] >= P_SOFT_CUT
    m &= (sep >= SEP_MIN) & (sepo >= SEP_MIN)
    m &= sep <= SEP_MAX
    m &= r['dca'] <= DCA_CUT
    par = r['open_a'] < THETA_PARALLEL
    m &= (~par) | (sepo < SEP_OUT_MAX_PARALLEL)
    gated = sepo > SEP_OUT_GATE
    m &= (~gated) | (r['collin'] > COLLIN_FRAC * L)   # through collinearity
    idx = np.where(m)[0]
    print(f"\n{len(idx)} events pass the collinearity cut -> {out_csv}")
    if len(idx) == 0:
        pd.DataFrame().to_csv(out_csv, index=False)
        return
    V = r['V_true'][idx]
    bx, by = local_transverse_xy(V)
    th, s, *_ = classify_points_with_basis(V)
    gpt = sep[idx] < SEP_IN_POINT_GATE
    pass_point = ((~gpt) | (r['pointing'][idx] < POINT_TIGHT_SEP_IN)) \
        & (r['pointing'][idx] < POINT_GLOBAL)
    pass_time = r['timing_chi2'][idx] < CHI2_TIMING_MAX
    df = pd.DataFrame({
        'idx': idx,
        'x_m': V[:, 0], 'y_m': V[:, 1], 'z_m': V[:, 2],
        'x_local_m': bx, 'y_local_m': by, 's_m': s, 'theta_deg': np.degrees(th),
        'p_mu_GeV': r['p_mu'][idx], 'p_e_GeV': r['p_e'][idx],
        'p_soft_GeV': r['p_soft'][idx],
        'sep_in_cm': sep[idx] * 100, 'sep_out_cm': sepo[idx] * 100,
        'collin_mm': r['collin'][idx] * 1000,
        'kink_deg': np.degrees(r['open_a'][idx]),
        'dca_cm': r['dca'][idx] * 100,
        'pointing_mrad': r['pointing'][idx] * 1000,
        'vtx_in': r['vtx_in'][idx],
        'timing_chi2': r['timing_chi2'][idx],
        'velo_pull_min': r['velo_pull_min'][idx],
        'rate_hz': r['w'][idx],
        'pass_vtx': r['vtx_in'][idx],
        'pass_pointing': pass_point,
        'pass_timing': pass_time,
        'pass_all': r['vtx_in'][idx] & pass_point & pass_time,
    })
    df.to_csv(out_csv, index=False)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--sampler', default='volume', choices=['volume', 'plane'],
                   help="'volume' = uniform-in-fiducial decay-vertex sampler "
                        "(efficient, default); 'plane' = throw-plane (cross-check).")
    p.add_argument('--n-rays', type=int, default=2_000_000,
                   help='Number of MC samples (volume decays, or thrown rays).')
    p.add_argument('--chunk', type=int, default=200_000,
                   help='Volume-sampler chunk size (lower = less peak memory).')
    p.add_argument('--streaming', action='store_true',
                   help='Memory-bounded volume run: accumulate the cutflow per '
                        'chunk + keep only collinearity passers (no per-event '
                        'plots/displays). Use for very large N (e.g. 2e7).')
    p.add_argument('--y-throw', type=float, default=50.0,
                   help='Y altitude (m) of the throw plane (plane sampler only).')
    p.add_argument('--theta-max', type=float, default=np.deg2rad(80))
    p.add_argument('--p-mean', type=float, default=1.0,
                   help='Mean of the input momentum spectrum (GeV/c).')
    p.add_argument('--spectrum', default='depth',
                   choices=['depth', 'empirical', 'exp', 'powerlaw', 'mono'],
                   help="'depth' = realistic post-overburden spectrum "
                        "(CosmicMuonFlux); 'empirical' = measured cavern-entry "
                        "spectrum from --momentum-root; others are toy spectra.")
    p.add_argument('--momentum-root', default=EMPIRICAL_ROOT_PATH,
                   help="ROOT file with the measured cavern-entry |p| spectrum "
                        "(TH1D h_muon_momentum_at_cavern_entry, MeV); used by "
                        "--spectrum empirical.")
    p.add_argument('--gen-pmin', type=float, default=None,
                   help="Restrict empirical generation to |p| > this (GeV) and "
                        "scale the rate by the window's spectral fraction "
                        "(importance sampling; skips useless high-p muons).")
    p.add_argument('--gen-pmax', type=float, default=None,
                   help="Restrict empirical generation to |p| < this (GeV). "
                        "Must be <= the histogram's upper edge (10 GeV).")
    p.add_argument('--timing-res-ns', type=float, default=SIGMA_T_DEFAULT * 1e9,
                   help='Per-hit timing resolution (ns).')
    p.add_argument('--velo-ncut', type=float, default=VELO_NCUT,
                   help='Per-track velocity-test threshold (reject pull < -n_cut).')
    p.add_argument('--collin-frac', type=float, default=COLLIN_FRAC,
                   help='Targeted collinearity cut as a fraction of L '
                        '(sep_out>L -> collin>frac*L). Signal ceiling is 0.5.')
    p.add_argument('--i-vertical', type=float, default=I_VERTICAL_DEFAULT,
                   help='Cosmic vertical intensity I_v [m^-2 s^-1 sr^-1] '
                        '(sets only the cos-theta weighting shape).')
    p.add_argument('--muon-rate-hz', type=float, default=MUON_RATE_HZ_DEFAULT,
                   help='Total muon-through-tunnel rate to anchor to (Hz). '
                        'Set 0 to use the I_v absolute scale instead.')
    p.add_argument('--livetime-years', type=float, default=10.0)
    p.add_argument('--seed', type=int, default=42)
    p.add_argument('--interactive', action='store_true')
    p.add_argument('--out', default='cosmic_decay_observables.png')
    p.add_argument('--plot-dir', default='cosmic_g16',
                   help='Dedicated folder for one-per-file variable plots '
                        '(signal-style), with the active cut thresholds drawn.')
    p.add_argument('--n-displays', type=int, default=4,
                   help='Number of passing-event displays to draw (0 = none).')
    p.add_argument('--display-dir', default='cosmic_decay_event_displays')
    p.add_argument('--horseshoe-out', default='cosmic_decay_horseshoe.png')
    p.add_argument('--collin-plot-out', default='cosmic_decay_collin_sepout.png')
    p.add_argument('--collin-ceiling-out',
                   default='cosmic_decay_collinearity.png')
    p.add_argument('--collin-csv', default='cosmic_decay_collin_passers.csv',
                   help='CSV dump of all events passing the collinearity cut.')
    return p.parse_args()


def main():
    args = parse_args()
    if not args.interactive and not os.environ.get('MPLBACKEND'):
        matplotlib.use('Agg')

    global VELO_NCUT, COLLIN_FRAC, EMPIRICAL_ROOT_PATH, GEN_PMIN, GEN_PMAX
    VELO_NCUT = args.velo_ncut
    COLLIN_FRAC = args.collin_frac
    EMPIRICAL_ROOT_PATH = args.momentum_root
    GEN_PMIN = args.gen_pmin
    GEN_PMAX = args.gen_pmax

    # Importance-sampling weight: when generating only a restricted momentum
    # window, scale the rate by that window's fraction of the full spectrum.
    spectrum_weight = 1.0
    if args.spectrum == 'empirical' and (GEN_PMIN is not None
                                         or GEN_PMAX is not None):
        edges = _load_empirical(EMPIRICAL_ROOT_PATH)['edges']
        wlo = 0.0 if GEN_PMIN is None else GEN_PMIN
        whi = edges[-1] if GEN_PMAX is None else GEN_PMAX
        spectrum_weight = empirical_window_fraction(wlo, whi)
        print(f"  GENERATION WINDOW: {wlo:.3g} < p < {whi:.3g} GeV  "
              f"(spectral fraction {spectrum_weight:.4e}; rate scaled by it)")

    print("=" * 70)
    print("COSMIC-MUON DECAY-IN-FLIGHT BACKGROUND")
    print("=" * 70)
    spec_desc = ('realistic depth spectrum' if args.spectrum == 'depth'
                 else f'measured cavern-entry spectrum ({args.momentum_root})'
                 if args.spectrum == 'empirical'
                 else f'{args.spectrum} (mean p = {args.p_mean} GeV/c)')
    print(f"  spectrum: {spec_desc}, "
          f"cos^2 theta to {np.degrees(args.theta_max):.0f} deg")
    print(f"  timing resolution: {args.timing_res_ns:.1f} ns/hit, "
          f"chi2 cut = {CHI2_TIMING_MAX:.0f}")
    print(f"  I_vertical = {args.i_vertical} m^-2 s^-1 sr^-1, "
          f"ctau_muon = {CTAU_MUON:.1f} m")
    print(f"  detector thickness = layer spacing = "
          f"{DETECTOR_THICKNESS*100:.0f} cm (from grendel_geometry)")
    print(f"  sampler: {args.sampler}{' (streaming)' if args.streaming else ''}")

    livetime_s = args.livetime_years * 365.25 * 24 * 3600
    if args.sampler == 'volume' and args.streaming:
        run_volume_streaming(
            args.n_rays, args.theta_max, args.p_mean, args.spectrum,
            args.i_vertical, HIT_RESOLUTION, args.seed, livetime_s,
            args.collin_csv, sigma_t=args.timing_res_ns * 1e-9,
            target_muon_rate_hz=args.muon_rate_hz, chunk=args.chunk,
            depth_e_min=0.3, spectrum_weight=spectrum_weight)
        return

    if args.sampler == 'volume':
        r = run_volume(args.n_rays, args.theta_max, args.p_mean,
                       args.spectrum, args.i_vertical, HIT_RESOLUTION, args.seed,
                       sigma_t=args.timing_res_ns * 1e-9,
                       target_muon_rate_hz=args.muon_rate_hz, chunk=args.chunk,
                       spectrum_weight=spectrum_weight)
    else:
        r = run(args.n_rays, args.y_throw, args.theta_max, args.p_mean,
                args.spectrum, args.i_vertical, HIT_RESOLUTION, args.seed,
                sigma_t=args.timing_res_ns * 1e-9,
                target_muon_rate_hz=args.muon_rate_hz)
    if r is None:
        print("No cosmics traversed the fiducial.")
        return

    livetime_s = args.livetime_years * 365.25 * 24 * 3600
    print(f"\n  Traversing muons:        {r['nv']:,}")
    print(f"  With reconstructed e:    {r['n_reco']:,}")
    print(f"  Total muon rate through fiducial: {r['muon_rate_hz']:.3e} Hz")
    print(f"  Decay-in-fiducial rate (reco):    {r['w'].sum():.3e} Hz")

    final, passmask, pre_timing = weighted_cutflow(r, livetime_s)
    dump_collin_passers(r, args.collin_csv)
    print(f"\n  ==> Cosmic muon-decay background after full selection:")
    print(f"        {final:.3e} Hz  =  {final*livetime_s:.3e} events "
          f"in {args.livetime_years:.0f} yr")

    # --- Where do the surviving decays happen relative to the wall? ---
    def wpct(x, wt, ps):
        o = np.argsort(x)
        cw = np.cumsum(wt[o]) / wt[o].sum()
        return np.interp(np.array(ps) / 100.0, cw, x[o])

    print("\nDecay-vertex proximity to the wall (rate-weighted med / 90% / 99%):")
    for label, msk in (('all reconstructed decays', np.ones(len(r['w']), bool)),
                       ('passing full selection', passmask)):
        if not msk.any():
            print(f"  {label}: (none)")
            continue
        Vp = r['V_true'][msk]
        wp = r['w'][msk]
        lever = np.linalg.norm(Vp - r['P_ei'][msk], axis=1)   # entry wall->vtx
        _, dwall, _ = mesh_fiducial.nearest.on_surface(Vp)     # perp to nearest wall
        lm = wpct(lever, wp, [50, 90, 99])
        dm = wpct(dwall, wp, [50, 90, 99])
        print(f"  {label} (N={int(msk.sum())}):")
        print(f"    muon path entry-wall -> vertex : "
              f"{lm[0]*100:6.1f} / {lm[1]*100:6.1f} / {lm[2]*100:6.1f} cm")
        print(f"    perp. distance to nearest wall : "
              f"{dm[0]*100:6.1f} / {dm[1]*100:6.1f} / {dm[2]*100:6.1f} cm")

    # Per-variable plots (one file each, signal-style) in a dedicated folder.
    pdir = args.plot_dir
    os.makedirs(pdir, exist_ok=True)
    make_variable_plots(r, pdir, interactive=args.interactive)
    make_plot(r, os.path.join(pdir, 'observables_overview.png'), args.interactive)

    make_decay_horseshoe(
        r, pre_timing, os.path.join(pdir, 'decay_horseshoe.png'),
        label=f'pre-timing selection (spacing {DETECTOR_THICKNESS*100:.0f} cm)',
        interactive=args.interactive)

    make_collin_sepout_plot(
        r, os.path.join(pdir, 'collin_vs_sepout.png'),
        label=f'spacing {DETECTOR_THICKNESS*100:.0f} cm',
        interactive=args.interactive)

    make_collin_ceiling_plot(r, os.path.join(pdir, 'collinearity_log.png'),
                             interactive=args.interactive)

    if args.n_displays > 0:
        print(f"\nDrawing up to {args.n_displays} passing-event displays "
              f"(of {int(passmask.sum())} that pass) -> {args.display_dir}/")
        make_cosmic_displays(r, passmask, args.n_displays, args.display_dir,
                             seed=args.seed)


if __name__ == '__main__':
    main()
