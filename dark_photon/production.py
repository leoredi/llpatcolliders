"""BC1 dark-photon production: A' four-vectors at the CMS IP (eps^2 = 1).

Two GRENDEL-relevant channels, mirroring the single-mechanism production
layers of BC4/BC10 but summed over the open channels at each mass:

  * Meson decay (low mass): light mesons pi0/eta/omega from pp SoftQCD
    (committed Pythia 8.315 spectrum, ``data/spectra/meson_softqcd.npz``) decay
    P -> A' gamma or omega -> A' pi0. The A' four-vector is built by the exact
    two-body decay in the meson rest frame boosted to the lab.
  * Drell-Yan (m_A >= 1.65 GeV): q qbar -> A' from a committed per-mass
    Pythia NewGaugeBoson spectrum (``data/spectra/dy/``); see
    ``generator/dy_pythia.py``. The lower bound is the NNPDF4.0 Q-validity
    floor; the perturbative calculation is not extrapolated below it.

Proton bremsstrahlung is deliberately omitted: like the BC10 kaon-tower note
in the paper, its forward approximation does not supply the transverse
acceptance GRENDEL requires.

Output: one headerless CSV per mass, columns ``weight, E, px, py, pz`` with
``weight`` the production cross section [pb] at ``eps^2 = 1`` carried by that
A' (so ``sum(weight)`` is the total sigma(pp -> A' + X) at eps^2 = 1, and the
sensitivity scan multiplies by the physical eps^2). This is exactly the
contract consumed by ``analysis.format_bridge.load_combined_csv``.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

_ROOT = Path(__file__).resolve().parent

from dark_photon import model  # noqa: E402

MESON_SPECTRUM = _ROOT / "data" / "spectra" / "meson_softqcd.npz"
DY_SPECTRUM_DIR = _ROOT / "data" / "spectra" / "dy"
MB_TO_PB = 1.0e9

# Parent -> (meson mass, recoil-partner mass, production-BR fn at eps^2 = 1,
# spectrum histogram key).
_MESON_CHANNELS = {
    "pi0":   (model.M_PI0,   0.0,        model.PRODUCTION_PARENTS["pi0"][1],   "hist_pi0"),
    "eta":   (model.M_ETA,   0.0,        model.PRODUCTION_PARENTS["eta"][1],   "hist_eta"),
    "omega": (model.M_OMEGA, model.M_PI0, model.PRODUCTION_PARENTS["omega"][1], "hist_omega"),
}


def _mass_label(m_a: float) -> str:
    return f"{m_a:.3f}".replace(".", "p")


def _sample_meson_4vectors(hist, pt_edges, y_edges, m_meson, n, rng):
    """Sample n meson four-vectors (E, px, py, pz) from a committed (pT,y) hist."""
    flat = hist.ravel().astype(float)
    cdf = np.cumsum(flat)
    if cdf[-1] <= 0:
        return None
    cdf /= cdf[-1]
    idx = np.clip(np.searchsorted(cdf, rng.random(n)), 0, len(flat) - 1)
    i_pt, i_y = np.unravel_index(idx, hist.shape)
    pt = rng.uniform(pt_edges[i_pt], pt_edges[i_pt + 1])
    y = rng.uniform(y_edges[i_y], y_edges[i_y + 1])
    phi = rng.uniform(0.0, 2.0 * np.pi, n)
    mt = np.sqrt(pt * pt + m_meson * m_meson)
    E = mt * np.cosh(y)
    pz = mt * np.sinh(y)
    px = pt * np.cos(phi)
    py = pt * np.sin(phi)
    return np.column_stack([E, px, py, pz])


def _two_body_decay_to_dp(parent_p4, m_parent, m_a, m_recoil, rng):
    """Boost an isotropic two-body decay parent -> A'(m_a) + X(m_recoil) to lab.

    Returns the A' lab four-vectors (n,4) as [E, px, py, pz].
    """
    n = len(parent_p4)
    # Rest-frame A' energy and momentum.
    E_star = (m_parent ** 2 + m_a ** 2 - m_recoil ** 2) / (2.0 * m_parent)
    lam = model._kallen(m_parent ** 2, m_a ** 2, m_recoil ** 2)
    p_star = np.sqrt(max(lam, 0.0)) / (2.0 * m_parent)
    cos_t = rng.uniform(-1.0, 1.0, n)
    sin_t = np.sqrt(np.clip(1.0 - cos_t ** 2, 0.0, None))
    phi = rng.uniform(0.0, 2.0 * np.pi, n)
    px = p_star * sin_t * np.cos(phi)
    py = p_star * sin_t * np.sin(phi)
    pz = p_star * cos_t
    E = np.full(n, E_star)
    # Boost from parent rest frame to lab using the parent four-vector.
    Ep, ppx, ppy, ppz = parent_p4[:, 0], parent_p4[:, 1], parent_p4[:, 2], parent_p4[:, 3]
    beta = np.column_stack([ppx, ppy, ppz]) / Ep[:, None]
    b2 = np.sum(beta * beta, axis=1)
    gamma = 1.0 / np.sqrt(np.clip(1.0 - b2, 1e-15, None))
    pvec = np.column_stack([px, py, pz])
    bp = np.sum(beta * pvec, axis=1)
    factor = np.where(b2 > 0, (gamma - 1.0) * bp / np.where(b2 > 0, b2, 1.0) + gamma * E, 0.0)
    p_lab = pvec + beta * factor[:, None]
    E_lab = gamma * (E + bp)
    return np.column_stack([E_lab, p_lab])


def generate_meson_channel(m_a, n_per_channel, rng, spectrum_path=MESON_SPECTRUM):
    """A' four-vectors + weights [pb, eps^2=1] from open meson channels."""
    data = np.load(spectrum_path)
    pt_edges, y_edges = data["pt_edges"], data["y_edges"]
    sigma_inel_pb = float(data["sigma_inel_mb"]) * MB_TO_PB
    out_p4, out_w = [], []
    for name, (m_par, m_rec, br_fn, hist_key) in _MESON_CHANNELS.items():
        if m_a >= m_par - m_rec:
            continue                      # channel kinematically closed
        br = br_fn(m_a)                    # production BR at eps^2 = 1
        if br <= 0:
            continue
        n_par_per_evt = float(data[f"n_{name}_per_evt"])
        # sigma(pp -> A' + X) via this meson [pb, eps^2 = 1]
        sigma_channel = sigma_inel_pb * n_par_per_evt * br
        parent_p4 = _sample_meson_4vectors(
            data[hist_key], pt_edges, y_edges, m_par, n_per_channel, rng)
        if parent_p4 is None:
            continue
        a_p4 = _two_body_decay_to_dp(parent_p4, m_par, m_a, m_rec, rng)
        out_p4.append(a_p4)
        out_w.append(np.full(len(a_p4), sigma_channel / len(a_p4)))
    if not out_p4:
        return np.empty((0, 4)), np.empty(0)
    return np.concatenate(out_p4), np.concatenate(out_w)


def generate_dy_channel(m_a, rng, spectrum_dir=DY_SPECTRUM_DIR):
    """A' four-vectors + weights [pb, eps^2=1] from the committed DY spectrum.

    Returns empty arrays when no committed DY spectrum exists for this mass
    (the meson channel then stands alone). See ``generator/dy_pythia.py``.
    """
    path = spectrum_dir / f"dy_{_mass_label(m_a)}.npz"
    if not path.exists():
        return np.empty((0, 4)), np.empty(0)
    data = np.load(path)
    pt_edges, y_edges = data["pt_edges"], data["y_edges"]
    sigma_dy_pb = float(data["sigma_dy_pb"])      # eps^2 = 1
    hist = data["hist"]
    flat = hist.ravel().astype(float)
    cdf = np.cumsum(flat)
    if cdf[-1] <= 0:
        return np.empty((0, 4)), np.empty(0)
    cdf /= cdf[-1]
    n = int(data["n_sample"]) if "n_sample" in data.files else int(flat.sum())
    n = max(n, 1)
    idx = np.clip(np.searchsorted(cdf, rng.random(n)), 0, len(flat) - 1)
    i_pt, i_y = np.unravel_index(idx, hist.shape)
    pt = rng.uniform(pt_edges[i_pt], pt_edges[i_pt + 1])
    y = rng.uniform(y_edges[i_y], y_edges[i_y + 1])
    phi = rng.uniform(0.0, 2.0 * np.pi, n)
    mt = np.sqrt(pt * pt + m_a * m_a)
    p4 = np.column_stack([mt * np.cosh(y), pt * np.cos(phi),
                          pt * np.sin(phi), mt * np.sinh(y)])
    return p4, np.full(n, sigma_dy_pb / n)


def generate(m_a, n_per_channel, rng):
    """Combined A' four-vectors + weights [pb, eps^2=1] over all open channels."""
    p4_m, w_m = generate_meson_channel(m_a, n_per_channel, rng)
    p4_d, w_d = generate_dy_channel(m_a, rng)
    p4 = np.concatenate([p4_m, p4_d]) if len(p4_d) else p4_m
    w = np.concatenate([w_m, w_d]) if len(w_d) else w_m
    return p4, w


def write_csv(m_a, out_dir, n_per_channel, rng):
    """Write the headerless ``weight,E,px,py,pz`` CSV for one mass."""
    p4, w = generate(m_a, n_per_channel, rng)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"mA_{_mass_label(m_a)}.csv"
    rows = np.column_stack([w, p4])       # weight, E, px, py, pz
    np.savetxt(path, rows, delimiter=",", fmt="%.8e")
    return path, len(w), float(w.sum())
