"""
production/decay_engine/kinematics.py

Vectorized phase-space decay functions for meson → HNL production.

Supports:
  - 2-body decay: M → ℓ N  (isotropic in parent rest frame)
  - 2-body decay with longitudinal polarisation: M → ℓ N (1 + α·cosθ)
  - 3-body decay: M → H' ℓ N  (flat Lorentz-invariant phase-space)
  - 3-body decay weighted by HNLCalc dBR/dq²/dE  (meson semileptonic, with
    form-factor matrix elements)
  - 3-body decay weighted by HNLCalc dBR/dE  (tau leptonic, V-A spectrum;
    angular content of the lepton/neutrino subsystem is filled isotropically
    since dBR/dE has already been integrated over it)

All functions are vectorized over arrays of parent 4-vectors.
"""

import numpy as np


def _boost_to_lab(p4_rest, parent_E, parent_px, parent_py, parent_pz):
    """
    Boost 4-vectors from parent rest frame to lab frame.

    Parameters
    ----------
    p4_rest : ndarray, shape (N, 4)
        4-vectors (E, px, py, pz) in parent rest frame.
    parent_E, parent_px, parent_py, parent_pz : ndarray, shape (N,)
        Parent 4-momentum in lab frame.

    Returns
    -------
    p4_lab : ndarray, shape (N, 4)
        Boosted 4-vectors in lab frame.
    """
    # Boost vector: β = p_parent / E_parent
    E_p = parent_E
    beta_x = parent_px / E_p
    beta_y = parent_py / E_p
    beta_z = parent_pz / E_p
    beta2 = beta_x**2 + beta_y**2 + beta_z**2
    beta2 = np.clip(beta2, 0, 1 - 1e-12)
    gamma = 1.0 / np.sqrt(1.0 - beta2)

    # Components in rest frame
    e_r = p4_rest[:, 0]
    px_r = p4_rest[:, 1]
    py_r = p4_rest[:, 2]
    pz_r = p4_rest[:, 3]

    # p · β
    bp = beta_x * px_r + beta_y * py_r + beta_z * pz_r

    # Handle beta2 ~ 0 (parent at rest)
    at_rest = beta2 < 1e-30

    e_lab = gamma * (e_r + bp)
    px_lab = px_r.copy()
    py_lab = py_r.copy()
    pz_lab = pz_r.copy()

    moving = ~at_rest
    if np.any(moving):
        k = (gamma[moving] - 1.0) * bp[moving] / beta2[moving]
        px_lab[moving] += beta_x[moving] * k + gamma[moving] * beta_x[moving] * e_r[moving]
        py_lab[moving] += beta_y[moving] * k + gamma[moving] * beta_y[moving] * e_r[moving]
        pz_lab[moving] += beta_z[moving] * k + gamma[moving] * beta_z[moving] * e_r[moving]

    return np.column_stack([e_lab, px_lab, py_lab, pz_lab])


def decay_2body(parent_E, parent_px, parent_py, parent_pz, m_parent, m1, m2, rng=None):
    """
    Isotropic 2-body decay: parent → daughter1 (mass m1) + daughter2 (mass m2).

    Parameters
    ----------
    parent_E, parent_px, parent_py, parent_pz : ndarray, shape (N,)
        Parent 4-momenta in lab frame.
    m_parent : float
        Parent mass in GeV.
    m1, m2 : float
        Daughter masses in GeV.
    rng : numpy.random.Generator, optional

    Returns
    -------
    d1 : ndarray, shape (N, 4) — (E, px, py, pz) of daughter 1 in lab frame
    d2 : ndarray, shape (N, 4) — (E, px, py, pz) of daughter 2 in lab frame
    """
    if rng is None:
        rng = np.random.default_rng()

    N = len(parent_E)

    # Momentum magnitude in parent rest frame (2-body kinematics)
    M = m_parent
    p_star = np.sqrt(
        np.maximum((M**2 - (m1 + m2)**2) * (M**2 - (m1 - m2)**2), 0.0)
    ) / (2.0 * M)

    # Isotropic direction in rest frame
    cos_theta = rng.uniform(-1.0, 1.0, N)
    sin_theta = np.sqrt(1.0 - cos_theta**2)
    phi = rng.uniform(0.0, 2 * np.pi, N)

    px_star = p_star * sin_theta * np.cos(phi)
    py_star = p_star * sin_theta * np.sin(phi)
    pz_star = p_star * cos_theta

    E1_star_val = np.sqrt(p_star**2 + m1**2)
    E2_star_val = np.sqrt(p_star**2 + m2**2)
    E1_star = np.full(N, E1_star_val)
    E2_star = np.full(N, E2_star_val)

    d1_rest = np.column_stack([E1_star, px_star, py_star, pz_star])
    d2_rest = np.column_stack([E2_star, -px_star, -py_star, -pz_star])

    # Boost to lab
    d1_lab = _boost_to_lab(d1_rest, parent_E, parent_px, parent_py, parent_pz)
    d2_lab = _boost_to_lab(d2_rest, parent_E, parent_px, parent_py, parent_pz)

    return d1_lab, d2_lab


def _sample_polar_cos(asym, n, rng):
    """Sample cosθ from the pdf ∝ (1 + asym·cosθ) on [-1, 1].

    Inverse-CDF: with F(c) = ((c+1) + (asym/2)(c²-1)) / 2 = u, solving the
    quadratic (asym/4)c² + (1/2)c + (1/2 - asym/4 - u) = 0 for c in [-1, 1].
    asym=0 reduces to the isotropic c = 2u-1.
    """
    u = rng.random(n)
    asym = float(np.clip(asym, -1.0, 1.0))
    if abs(asym) < 1e-12:
        return 2.0 * u - 1.0
    a = asym / 4.0
    b = 0.5
    c = 0.5 - asym / 4.0 - u
    disc = np.sqrt(np.maximum(b * b - 4.0 * a * c, 0.0))
    # Root that lands in [-1, 1] (the physical branch for |asym|<=1).
    cos_theta = (-b + disc) / (2.0 * a)
    return np.clip(cos_theta, -1.0, 1.0)


def decay_2body_polarized(parent_E, parent_px, parent_py, parent_pz,
                          m_parent, m1, m2, asymmetry=0.0, rng=None):
    """2-body decay of a longitudinally polarized parent: parent → d1(m1) + d2(m2).

    The *second* daughter (mass ``m2``) is the spin-analyzed particle: its polar
    angle in the parent rest frame is drawn from ``1 + asymmetry·cosθ`` about the
    parent's lab-momentum direction (the longitudinal-polarization axis), rather
    than isotropically. ``asymmetry = analyzing_power × P_parent ∈ [-1, 1]``;
    ``asymmetry = 0`` is exactly the isotropic :func:`decay_2body`.

    Returns
    -------
    d1, d2 : ndarray, shape (N, 4) — lab-frame (E, px, py, pz). d2 is the
        analyzed daughter (mass m2).
    """
    if rng is None:
        rng = np.random.default_rng()

    N = len(parent_E)
    M = m_parent
    p_star = np.sqrt(
        np.maximum((M**2 - (m1 + m2)**2) * (M**2 - (m1 - m2)**2), 0.0)
    ) / (2.0 * M)

    # Polarization axis = parent lab-momentum direction (longitudinal pol).
    p_mag = np.sqrt(parent_px**2 + parent_py**2 + parent_pz**2)
    safe = np.where(p_mag > 0, p_mag, 1.0)
    nx, ny, nz = parent_px / safe, parent_py / safe, parent_pz / safe
    # Parents at rest get an arbitrary fixed axis (asymmetry then irrelevant).
    at_rest = p_mag == 0
    nz = np.where(at_rest, 1.0, nz)
    nx = np.where(at_rest, 0.0, nx)
    ny = np.where(at_rest, 0.0, ny)

    cos_theta = _sample_polar_cos(asymmetry, N, rng)
    sin_theta = np.sqrt(np.maximum(1.0 - cos_theta**2, 0.0))
    phi = rng.uniform(0.0, 2 * np.pi, N)

    # Build an orthonormal basis (n, e1, e2) per event to place d2 at angle θ
    # from the polarization axis n.
    # Choose a reference not parallel to n.
    ref_z = np.abs(nz) < 0.9
    rx = np.where(ref_z, 0.0, 1.0)
    ry = np.zeros(N)
    rz = np.where(ref_z, 1.0, 0.0)
    # e1 = n × ref, normalized
    e1x = ny * rz - nz * ry
    e1y = nz * rx - nx * rz
    e1z = nx * ry - ny * rx
    e1n = np.sqrt(e1x**2 + e1y**2 + e1z**2)
    e1n = np.where(e1n > 0, e1n, 1.0)
    e1x, e1y, e1z = e1x / e1n, e1y / e1n, e1z / e1n
    # e2 = n × e1
    e2x = ny * e1z - nz * e1y
    e2y = nz * e1x - nx * e1z
    e2z = nx * e1y - ny * e1x

    dirx = cos_theta * nx + sin_theta * (np.cos(phi) * e1x + np.sin(phi) * e2x)
    diry = cos_theta * ny + sin_theta * (np.cos(phi) * e1y + np.sin(phi) * e2y)
    dirz = cos_theta * nz + sin_theta * (np.cos(phi) * e1z + np.sin(phi) * e2z)

    px2 = p_star * dirx
    py2 = p_star * diry
    pz2 = p_star * dirz
    E2_star = np.sqrt(p_star**2 + m2**2)
    E1_star = np.sqrt(p_star**2 + m1**2)

    d2_rest = np.column_stack([np.full(N, 0.0) + E2_star, px2, py2, pz2])
    d1_rest = np.column_stack([np.full(N, 0.0) + E1_star, -px2, -py2, -pz2])

    d1_lab = _boost_to_lab(d1_rest, parent_E, parent_px, parent_py, parent_pz)
    d2_lab = _boost_to_lab(d2_rest, parent_E, parent_px, parent_py, parent_pz)
    return d1_lab, d2_lab


def decay_3body_flat(parent_E, parent_px, parent_py, parent_pz,
                     m_parent, m1, m2, m3, rng=None):
    """
    Flat 3-body phase-space decay (Dalitz): parent → d1(m1) + d2(m2) + d3(m3).

    No matrix-element weighting — uniform sampling in phase space using the
    sequential 2-body method:
      1. Sample m23 uniformly in [(m2+m3), (M-m1)]
      2. Decay parent → d1 + (23-system) as 2-body
      3. Decay (23-system) → d2 + d3 as 2-body

    This produces flat phase-space (up to the Jacobian |p*_1| × |p*_23|
    which is already the phase-space measure).

    Parameters
    ----------
    parent_E, parent_px, parent_py, parent_pz : ndarray, shape (N,)
    m_parent : float
    m1, m2, m3 : float
        Daughter masses in GeV.
    rng : numpy.random.Generator, optional

    Returns
    -------
    d1, d2, d3 : ndarray, shape (N, 4) — lab frame 4-vectors
    """
    if rng is None:
        rng = np.random.default_rng()

    N = len(parent_E)
    M = m_parent

    # Sample invariant mass of (23) system with correct phase-space weight.
    # The phase-space weight for flat sampling in m23 is |p*_1| × |p*_23|.
    # We use accept-reject on m23 to get the correct distribution.
    m23_min = m2 + m3
    m23_max = M - m1
    if m23_max <= m23_min:
        raise RuntimeError("No available 3-body phase space (at or below threshold).")

    # Use accept-reject: weight = p_star_1(m23) × p_star_23(m23)
    # Accept-reject envelope: max weight (found numerically per-event, but
    # since all events have same masses, one max suffices).

    # Precompute max weight by sampling densely
    m23_test = np.linspace(m23_min + 1e-6, m23_max - 1e-6, 1000)
    lam1 = (M**2 - m1**2 - m23_test**2)**2 - 4 * m1**2 * m23_test**2
    lam23 = (m23_test**2 - m2**2 - m3**2)**2 - 4 * m2**2 * m3**2
    p1_test = np.sqrt(np.maximum(lam1, 0.0)) / (2.0 * M)
    p23_test = np.sqrt(np.maximum(lam23, 0.0)) / (2.0 * np.maximum(m23_test, 1e-15))
    w_test = p1_test * p23_test
    w_max = w_test.max() * 1.01  # small safety margin
    if w_max <= 0:
        raise RuntimeError("No available 3-body phase space weight (at threshold).")

    # Accept-reject loop
    m23 = np.empty(N)
    accepted = np.zeros(N, dtype=bool)
    n_iter = 0
    max_iter = 10000
    while not accepted.all():
        n_iter += 1
        if n_iter > max_iter:
            raise RuntimeError("Failed to sample 3-body decay after max iterations.")
        n_need = (~accepted).sum()
        m23_cand = rng.uniform(m23_min, m23_max, n_need)
        lam1_c = (M**2 - m1**2 - m23_cand**2)**2 - 4 * m1**2 * m23_cand**2
        lam23_c = (m23_cand**2 - m2**2 - m3**2)**2 - 4 * m2**2 * m3**2
        p1_cand = np.sqrt(np.maximum(lam1_c, 0.0)) / (2.0 * M)
        p23_cand = np.sqrt(np.maximum(lam23_c, 0.0)) / (2.0 * np.maximum(m23_cand, 1e-15))
        w_cand = p1_cand * p23_cand
        u = rng.random(n_need) * w_max
        accept_mask = u < w_cand
        idx = np.where(~accepted)[0]
        m23[idx[accept_mask]] = m23_cand[accept_mask]
        accepted[idx[accept_mask]] = True

    # Step 1: decay parent → d1 + system23 in parent rest frame
    # Momentum of d1 in parent rest frame
    p1_star = np.sqrt(np.maximum(
        (M**2 - (m1 + m23)**2) * (M**2 - (m1 - m23)**2), 0.0
    )) / (2.0 * M)

    cos_th1 = rng.uniform(-1.0, 1.0, N)
    sin_th1 = np.sqrt(1.0 - cos_th1**2)
    phi1 = rng.uniform(0.0, 2 * np.pi, N)

    d1_px = p1_star * sin_th1 * np.cos(phi1)
    d1_py = p1_star * sin_th1 * np.sin(phi1)
    d1_pz = p1_star * cos_th1
    d1_E = np.sqrt(p1_star**2 + m1**2)

    sys23_E = np.sqrt(p1_star**2 + m23**2)
    sys23_px = -d1_px
    sys23_py = -d1_py
    sys23_pz = -d1_pz

    # Boost d1 and sys23 from parent rest frame to lab
    d1_rest = np.column_stack([d1_E, d1_px, d1_py, d1_pz])
    d1_lab = _boost_to_lab(d1_rest, parent_E, parent_px, parent_py, parent_pz)

    sys23_rest = np.column_stack([sys23_E, sys23_px, sys23_py, sys23_pz])
    sys23_lab = _boost_to_lab(sys23_rest, parent_E, parent_px, parent_py, parent_pz)

    # Step 2: decay sys23 → d2 + d3 in sys23 rest frame
    p23_star = np.sqrt(np.maximum(
        (m23**2 - (m2 + m3)**2) * (m23**2 - (m2 - m3)**2), 0.0
    )) / (2.0 * m23)

    cos_th2 = rng.uniform(-1.0, 1.0, N)
    sin_th2 = np.sqrt(1.0 - cos_th2**2)
    phi2 = rng.uniform(0.0, 2 * np.pi, N)

    d2_px_r = p23_star * sin_th2 * np.cos(phi2)
    d2_py_r = p23_star * sin_th2 * np.sin(phi2)
    d2_pz_r = p23_star * cos_th2
    d2_E_r = np.sqrt(p23_star**2 + m2**2)
    d3_E_r = np.sqrt(p23_star**2 + m3**2)

    d2_in_sys23 = np.column_stack([d2_E_r, d2_px_r, d2_py_r, d2_pz_r])
    d3_in_sys23 = np.column_stack([d3_E_r, -d2_px_r, -d2_py_r, -d2_pz_r])

    # Boost d2, d3 from sys23 rest frame to lab
    d2_lab = _boost_to_lab(
        d2_in_sys23,
        sys23_lab[:, 0], sys23_lab[:, 1], sys23_lab[:, 2], sys23_lab[:, 3],
    )
    d3_lab = _boost_to_lab(
        d3_in_sys23,
        sys23_lab[:, 0], sys23_lab[:, 1], sys23_lab[:, 2], sys23_lab[:, 3],
    )

    return d1_lab, d2_lab, d3_lab


# ---------------------------------------------------------------------------
# Matrix-element-weighted 3-body samplers
#
# These replace decay_3body_flat for the channels where HNLCalc provides a
# differential branching expression. The flat sampler integrates correctly
# (uniform LIPS) but distributes events flat in the Dalitz plot, which biases
# the HNL boost spectrum once form factors / V-A are folded in. The weighted
# samplers below use accept-reject on the HNLCalc dBR string so the
# kinematics carry the same matrix-element shape as the rate.
# ---------------------------------------------------------------------------


def _eval_dbr(dbr_expr, locals_dict):
    """Vectorised eval of an HNLCalc differential-BR string.

    HNLCalc constructs the dBR as a Python expression string whose free names
    are q, q2, energy, mass / m3 (both refer to the HNL mass — the
    pseudoscalar template uses ``mass``, the vector template uses ``m3``),
    coupling, and np. We bind both ``mass`` and ``m3`` to the same HNL-mass
    scalar so either template form evaluates, and pass np in globals so the
    expression works element-wise on numpy arrays.
    """
    if "mass" in locals_dict and "m3" not in locals_dict:
        locals_dict = dict(locals_dict, m3=locals_dict["mass"])
    return eval(dbr_expr, {"np": np, "__builtins__": {}}, locals_dict)


def _orthonormal_basis(nx, ny, nz):
    """For each unit vector (nx, ny, nz), pick (e1, e2) perpendicular to it.

    Used to place a daughter at polar angle θ wrt the boost axis and uniform
    azimuth φ around it. Same logic as decay_2body_polarized.
    """
    n = len(nx)
    ref_z = np.abs(nz) < 0.9
    rx = np.where(ref_z, 0.0, 1.0)
    ry = np.zeros(n)
    rz = np.where(ref_z, 1.0, 0.0)
    e1x = ny * rz - nz * ry
    e1y = nz * rx - nx * rz
    e1z = nx * ry - ny * rx
    e1n = np.sqrt(e1x ** 2 + e1y ** 2 + e1z ** 2)
    e1n = np.where(e1n > 0, e1n, 1.0)
    e1x, e1y, e1z = e1x / e1n, e1y / e1n, e1z / e1n
    e2x = ny * e1z - nz * e1y
    e2y = nz * e1x - nx * e1z
    e2z = nx * e1y - ny * e1x
    return (e1x, e1y, e1z), (e2x, e2y, e2z)


def decay_3body_weighted_dq2dE(parent_E, parent_px, parent_py, parent_pz,
                               m_parent, m1, m2, m3, dbr_expr,
                               coupling=1.0, rng=None,
                               ceiling_grid=80, safety=1.2,
                               max_iter=10000):
    """Matrix-element-weighted 3-body decay parent → d1(m1) + d2(m2) + d3(m3).

    Samples (q², E_N) ∝ dBR/dq²/dE with q² = (p_d2 + p_d3)² and E_N the d3
    energy in the parent rest frame, using accept-reject on the HNLCalc
    expression ``dbr_expr``. Conventions match HNLCalc.integrate_3body_br
    with integration="dq2dE": d3 is the HNL (mass = m3 = HNL mass), d2 is
    the charged lepton, d1 is the daughter meson.

    Other degrees of freedom (d1 direction in parent rest frame and the HNL
    azimuth in the (d2,d3) rest frame) are uniform — they are not constrained
    by the q²/E_N differential.

    Parameters
    ----------
    dbr_expr : str
        HNLCalc dBR/dq²/dE string from get_3body_dbr_pseudoscalar /
        get_3body_dbr_vector. Eval'd with locals q, q2, energy, mass, coupling.
    coupling : float
        Coupling value used when HNLCalc built the expression (the same
        constant the rate integrator passes through; 1.0 for our U²=1
        bookkeeping).
    ceiling_grid : int
        Grid resolution per axis for the dBR ceiling search before
        accept-reject. The ceiling is multiplied by ``safety`` before use.
    """
    if rng is None:
        rng = np.random.default_rng()

    N = len(parent_E)
    M = m_parent

    q2min = (m2 + m3) ** 2
    q2max = (M - m1) ** 2
    if q2max <= q2min:
        raise RuntimeError("No 3-body phase space (at or below threshold).")

    # Ceiling search on a (q², u) grid where u in [0,1] maps to E_N within
    # the q²-dependent physical range. The acceptance weight at each (q², E)
    # is dbr × (ENmax(q²) - ENmin(q²)); the (q2max - q2min) factor is constant
    # and absorbs into the ceiling.
    q2_grid = np.linspace(q2min + 1e-9, q2max - 1e-9, ceiling_grid)
    u_grid = np.linspace(1e-6, 1.0 - 1e-6, ceiling_grid)
    Q2, U = np.meshgrid(q2_grid, u_grid, indexing="ij")
    Q = np.sqrt(Q2)
    E2st = (Q2 - m2 ** 2 + m3 ** 2) / (2.0 * Q)
    E3st = (M ** 2 - Q2 - m1 ** 2) / (2.0 * Q)
    # Boundary on m23² from the second-stage decay (see HNLCalc integrator).
    A = E2st + E3st
    B1 = np.sqrt(np.maximum(E2st ** 2 - m3 ** 2, 0.0))
    B2 = np.sqrt(np.maximum(E3st ** 2 - m1 ** 2, 0.0))
    m232min = A ** 2 - (B1 + B2) ** 2
    m232max = A ** 2 - (B1 - B2) ** 2
    ENmin_grid = (m232min + Q2 - m2 ** 2 - m1 ** 2) / (2.0 * M)
    ENmax_grid = (m232max + Q2 - m2 ** 2 - m1 ** 2) / (2.0 * M)
    E_grid = ENmin_grid + U * (ENmax_grid - ENmin_grid)
    w_grid = _eval_dbr(
        dbr_expr,
        {"q": Q, "q2": Q2, "energy": E_grid, "mass": m3, "coupling": coupling},
    )
    w_grid = np.where(np.isfinite(w_grid) & (w_grid > 0), w_grid, 0.0)
    span_grid = np.maximum(ENmax_grid - ENmin_grid, 0.0)
    ceiling = float((w_grid * span_grid).max()) * safety
    if not np.isfinite(ceiling) or ceiling <= 0:
        raise RuntimeError(
            f"dbr ceiling not positive (max={ceiling}). dbr_expr may be "
            f"degenerate at this (M, m1, m2, m3) = "
            f"({M}, {m1}, {m2}, {m3})."
        )

    # Vectorised accept-reject.
    q2_acc = np.empty(N)
    E_acc = np.empty(N)
    accepted = np.zeros(N, dtype=bool)
    for _ in range(max_iter):
        if accepted.all():
            break
        n_need = int((~accepted).sum())
        q2c = rng.uniform(q2min, q2max, n_need)
        qc = np.sqrt(q2c)
        E2 = (q2c - m2 ** 2 + m3 ** 2) / (2.0 * qc)
        E3 = (M ** 2 - q2c - m1 ** 2) / (2.0 * qc)
        A = E2 + E3
        B1 = np.sqrt(np.maximum(E2 ** 2 - m3 ** 2, 0.0))
        B2 = np.sqrt(np.maximum(E3 ** 2 - m1 ** 2, 0.0))
        mmin = A ** 2 - (B1 + B2) ** 2
        mmax = A ** 2 - (B1 - B2) ** 2
        ENmin = (mmin + q2c - m2 ** 2 - m1 ** 2) / (2.0 * M)
        ENmax = (mmax + q2c - m2 ** 2 - m1 ** 2) / (2.0 * M)
        uc = rng.random(n_need)
        Ec = ENmin + uc * (ENmax - ENmin)
        span = ENmax - ENmin
        wc = _eval_dbr(
            dbr_expr,
            {"q": qc, "q2": q2c, "energy": Ec, "mass": m3, "coupling": coupling},
        )
        wc = np.where(np.isfinite(wc) & (wc > 0), wc, 0.0)
        u_acc = rng.uniform(0.0, ceiling, n_need)
        ok = u_acc < (wc * span)
        idx = np.where(~accepted)[0]
        q2_acc[idx[ok]] = q2c[ok]
        E_acc[idx[ok]] = Ec[ok]
        accepted[idx[ok]] = True
    if not accepted.all():
        raise RuntimeError(
            f"dq2dE accept-reject failed to converge ({(~accepted).sum()} "
            f"of {N} unaccepted after {max_iter} iterations; ceiling may be "
            f"too loose)."
        )

    # Reconstruct full kinematics in parent rest frame from (q², E_N).
    q_acc = np.sqrt(q2_acc)
    # d1 in parent rest frame: isotropic.
    lam_d1 = (M ** 2 - (m1 + q_acc) ** 2) * (M ** 2 - (m1 - q_acc) ** 2)
    p_d1 = np.sqrt(np.maximum(lam_d1, 0.0)) / (2.0 * M)
    E_d1 = np.sqrt(p_d1 ** 2 + m1 ** 2)
    cos_t1 = rng.uniform(-1.0, 1.0, N)
    sin_t1 = np.sqrt(np.maximum(1.0 - cos_t1 ** 2, 0.0))
    phi1 = rng.uniform(0.0, 2 * np.pi, N)
    p_d1x = p_d1 * sin_t1 * np.cos(phi1)
    p_d1y = p_d1 * sin_t1 * np.sin(phi1)
    p_d1z = p_d1 * cos_t1

    # sys23 (lN system) 4-vec in parent rest frame.
    E_sys23 = M - E_d1
    p_sys23x = -p_d1x
    p_sys23y = -p_d1y
    p_sys23z = -p_d1z
    # Boost params from sys23 rest frame to parent rest frame:
    gamma = E_sys23 / q_acc
    p_sys23_mag = p_d1  # magnitude, opposite direction
    beta = np.where(E_sys23 > 0, p_sys23_mag / E_sys23, 0.0)

    # In sys23 rest frame, HNL energy and momentum magnitude.
    E_d3_star = (q2_acc + m3 ** 2 - m2 ** 2) / (2.0 * q_acc)
    p_d3_star = np.sqrt(np.maximum(E_d3_star ** 2 - m3 ** 2, 0.0))

    # Recover cos(θ*) of HNL wrt sys23 boost direction from E_N:
    #   E_N = γ E*_d3 + β γ |p*_d3| cos(θ*)
    bg_p3 = beta * gamma * p_d3_star
    safe = bg_p3 > 1e-15
    cos_ts = np.where(safe, (E_acc - gamma * E_d3_star) / np.where(safe, bg_p3, 1.0), 0.0)
    cos_ts = np.clip(cos_ts, -1.0, 1.0)
    sin_ts = np.sqrt(np.maximum(1.0 - cos_ts ** 2, 0.0))
    phi_s = rng.uniform(0.0, 2 * np.pi, N)

    # sys23 boost axis is along -d1 (sys23 momentum direction in parent frame).
    safe_d1 = p_d1 > 0
    nx = np.where(safe_d1, p_sys23x / np.where(safe_d1, p_sys23_mag, 1.0), 0.0)
    ny = np.where(safe_d1, p_sys23y / np.where(safe_d1, p_sys23_mag, 1.0), 0.0)
    nz = np.where(safe_d1, p_sys23z / np.where(safe_d1, p_sys23_mag, 1.0), 1.0)
    (e1x, e1y, e1z), (e2x, e2y, e2z) = _orthonormal_basis(nx, ny, nz)
    cphi = np.cos(phi_s)
    sphi = np.sin(phi_s)
    dx = cos_ts * nx + sin_ts * (cphi * e1x + sphi * e2x)
    dy = cos_ts * ny + sin_ts * (cphi * e1y + sphi * e2y)
    dz = cos_ts * nz + sin_ts * (cphi * e1z + sphi * e2z)

    p_d3x_star = p_d3_star * dx
    p_d3y_star = p_d3_star * dy
    p_d3z_star = p_d3_star * dz
    d3_in_sys23 = np.column_stack([
        np.full(N, 0.0) + E_d3_star,
        p_d3x_star, p_d3y_star, p_d3z_star,
    ])

    # d2 (lepton) in sys23 rest frame: back-to-back with d3.
    E_d2_star = np.sqrt(p_d3_star ** 2 + m2 ** 2)
    d2_in_sys23 = np.column_stack([
        np.full(N, 0.0) + E_d2_star,
        -p_d3x_star, -p_d3y_star, -p_d3z_star,
    ])

    # Boost d2, d3 from sys23 rest frame to parent rest frame.
    d3_parent = _boost_to_lab(d3_in_sys23, E_sys23, p_sys23x, p_sys23y, p_sys23z)
    d2_parent = _boost_to_lab(d2_in_sys23, E_sys23, p_sys23x, p_sys23y, p_sys23z)
    d1_parent = np.column_stack([E_d1, p_d1x, p_d1y, p_d1z])

    # Boost all three from parent rest frame to lab.
    d1_lab = _boost_to_lab(d1_parent, parent_E, parent_px, parent_py, parent_pz)
    d2_lab = _boost_to_lab(d2_parent, parent_E, parent_px, parent_py, parent_pz)
    d3_lab = _boost_to_lab(d3_parent, parent_E, parent_px, parent_py, parent_pz)

    return d1_lab, d2_lab, d3_lab


def decay_3body_weighted_dE(parent_E, parent_px, parent_py, parent_pz,
                            m_parent, m1, m2, m3, dbr_expr,
                            coupling=1.0, rng=None,
                            ceiling_grid=400, safety=1.2,
                            max_iter=10000):
    """3-body decay weighted by dBR/dE_3 (tau leptonic mode).

    For tau → lepton + neutrino + N, HNLCalc supplies dBR/dE already
    integrated over the lepton/neutrino angular content. We accept-reject on
    E_N (HNL energy in tau rest frame), then place the HNL direction
    isotropically in the tau rest frame and fill the lepton/neutrino
    subsystem isotropically in its own rest frame.

    The HNL lab kinematics are then exact under the assumed matrix element;
    the lepton/neutrino angular distribution relative to the HNL is *not*
    faithful (it has been integrated out by the dE differential). Downstream
    consumers in this codebase only read the HNL 4-vector, so this is OK.
    If you ever need the lepton kinematics with V-A correlations, switch to
    a dq²/dE expression and use decay_3body_weighted_dq2dE instead.
    """
    if rng is None:
        rng = np.random.default_rng()

    N = len(parent_E)
    M = m_parent

    Emin = m3
    Emax = (M ** 2 + m3 ** 2 - (m1 + m2) ** 2) / (2.0 * M)
    if Emax <= Emin:
        raise RuntimeError("No 3-body phase space for dE sampler.")

    # 1D ceiling.
    E_grid = np.linspace(Emin + 1e-9, Emax - 1e-9, ceiling_grid)
    w_grid = _eval_dbr(
        dbr_expr,
        {"energy": E_grid, "mass": m3, "coupling": coupling},
    )
    w_grid = np.where(np.isfinite(w_grid) & (w_grid > 0), w_grid, 0.0)
    ceiling = float(w_grid.max()) * safety
    if not np.isfinite(ceiling) or ceiling <= 0:
        raise RuntimeError(
            f"dE ceiling not positive (max={ceiling}); "
            f"(M, m1, m2, m3) = ({M}, {m1}, {m2}, {m3})."
        )

    # Accept-reject on E_N.
    E_acc = np.empty(N)
    accepted = np.zeros(N, dtype=bool)
    for _ in range(max_iter):
        if accepted.all():
            break
        n_need = int((~accepted).sum())
        Ec = rng.uniform(Emin, Emax, n_need)
        wc = _eval_dbr(
            dbr_expr,
            {"energy": Ec, "mass": m3, "coupling": coupling},
        )
        wc = np.where(np.isfinite(wc) & (wc > 0), wc, 0.0)
        u = rng.uniform(0.0, ceiling, n_need)
        ok = u < wc
        idx = np.where(~accepted)[0]
        E_acc[idx[ok]] = Ec[ok]
        accepted[idx[ok]] = True
    if not accepted.all():
        raise RuntimeError("dE accept-reject failed to converge.")

    # HNL direction isotropic in parent rest frame; |p_N| set by E_N.
    p_d3 = np.sqrt(np.maximum(E_acc ** 2 - m3 ** 2, 0.0))
    cos_t3 = rng.uniform(-1.0, 1.0, N)
    sin_t3 = np.sqrt(np.maximum(1.0 - cos_t3 ** 2, 0.0))
    phi3 = rng.uniform(0.0, 2 * np.pi, N)
    p_d3x = p_d3 * sin_t3 * np.cos(phi3)
    p_d3y = p_d3 * sin_t3 * np.sin(phi3)
    p_d3z = p_d3 * cos_t3
    d3_parent = np.column_stack([E_acc, p_d3x, p_d3y, p_d3z])

    # (d1, d2) subsystem in parent rest frame: invariant mass squared
    # m12² = M² + m3² - 2 M E_N, momentum equal and opposite to d3.
    m12sq = M ** 2 + m3 ** 2 - 2.0 * M * E_acc
    m12 = np.sqrt(np.maximum(m12sq, 0.0))
    E_sys12 = M - E_acc
    # Decay (d1, d2) isotropically in their own rest frame.
    lam12 = (m12sq - (m1 + m2) ** 2) * (m12sq - (m1 - m2) ** 2)
    p_star = np.sqrt(np.maximum(lam12, 0.0)) / (2.0 * np.maximum(m12, 1e-15))
    cos_a = rng.uniform(-1.0, 1.0, N)
    sin_a = np.sqrt(np.maximum(1.0 - cos_a ** 2, 0.0))
    phi_a = rng.uniform(0.0, 2 * np.pi, N)
    p1x = p_star * sin_a * np.cos(phi_a)
    p1y = p_star * sin_a * np.sin(phi_a)
    p1z = p_star * cos_a
    E1_star = np.sqrt(p_star ** 2 + m1 ** 2)
    E2_star = np.sqrt(p_star ** 2 + m2 ** 2)
    d1_in_sys12 = np.column_stack([E1_star, p1x, p1y, p1z])
    d2_in_sys12 = np.column_stack([E2_star, -p1x, -p1y, -p1z])
    # Boost (d1, d2) from sys12 rest frame (3-momentum -p_d3) to parent rest.
    d1_parent = _boost_to_lab(d1_in_sys12, E_sys12, -p_d3x, -p_d3y, -p_d3z)
    d2_parent = _boost_to_lab(d2_in_sys12, E_sys12, -p_d3x, -p_d3y, -p_d3z)

    # Boost everything from parent rest frame to lab.
    d1_lab = _boost_to_lab(d1_parent, parent_E, parent_px, parent_py, parent_pz)
    d2_lab = _boost_to_lab(d2_parent, parent_E, parent_px, parent_py, parent_pz)
    d3_lab = _boost_to_lab(d3_parent, parent_E, parent_px, parent_py, parent_pz)

    return d1_lab, d2_lab, d3_lab
