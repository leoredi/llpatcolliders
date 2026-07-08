r"""BC10 fermiophilic-ALP model layer (couplings -> production, lifetime, BRs).

PBC benchmark **BC10**: a pseudoscalar axion-like particle ``a`` whose dominant
coupling is to SM fermions (the genuine charged-track "ALP" benchmark).  The
single parameter that this analysis sets a limit on is the inverse decay
constant ``1/f`` (universal Wilson coefficients ``c_f = 1``); it controls
production, lifetime, and the visible branching ratios *simultaneously*.

Coupling convention (fixes the meaning of the 1/f axis): the effective
pseudoscalar Yukawa g_aff = c_f m_f / f of Bauer-Neubert-Thamm
(arXiv:1708.00443), i.e.

    Gamma(a -> f fbar) = N_c (c_f m_f)^2 m_a (1/f)^2 beta / (8 pi).

The state-of-the-art BC10 phenomenology paper, Garcia-Kahlhoefer-Ovchynnikov-
Zaporozhchenko (GKOZ) arXiv:2310.03524, normalises its Lagrangian with
d_mu a/(2 f_GKOZ) instead, so our axis maps onto theirs as
1/f (here) = 2/f_GKOZ.  All external inputs below have been converted into the
BNT convention and cross-checked on physical (convention-free) anchors; see
data/alpinist/PROVENANCE.md.

External inputs (replacing the analytic placeholders of the first version):

1. **Hadronic + gamma gamma decay widths (data-driven).**  Digitized GKOZ
   tables shipped by ALPINIST (github.com/jjerhot/ALPINIST, BSD-3, pinned in
   data/alpinist/COMMIT_SHA.txt): per-channel Gamma/(1/f)^2 vs m_a over
   0.01--3.01 GeV, with the physical eta/eta' mixing poles and the 2 m_c
   perturbative onset.  Above the table ceiling the perturbative quark-level
   sum continues the width, normalised to the table at the seam (the raw sum
   is ~1.5x the table there -- charm mass-scheme sensitivity).  Anchors
   validating the normalisation: BR(a->mumu) ~ 9% at m_a = 1 GeV (GKOZ state
   "<10% for m_a >~ 1 GeV") and the quark-level sum matching the table's own
   2 m_c onset step to ~5%.

2. **b -> s a production coupling (one-loop RG, finite terms).**  The
   flavour-violating coefficient is evaluated with the ALPINIST implementation
   of the GKOZ RG equations (tools/compute_cbs_alpinist.py, ported under
   BSD-3):

       g_bs = CBS_EFF * (1/f),   CBS_EFF = 3.518383e-4  (Lambda_UV = 1 TeV)

   replacing the leading-log-with-C_TOP=1 placeholder (which was 2.5x larger
   in amplitude -- the finite one-loop and RG terms partially cancel the small
   log at 1 TeV).  Cross-check: GKOZ Table 1 quotes |C_bs| = 1.8e-3 at
   f_GKOZ = 1 GeV, i.e. c_bs = 2 |C_bs| / m_b ~ 7.5e-4 -- 8.6e-4 (m_b pole vs
   MSbar), whose BNT-convention half is 3.7e-4 -- 4.3e-4, within 6--20% of
   CBS_EFF.  The island's production-limited (lower) edge inherits this
   residual as <~20% in 1/f.

3. **Production channels: the full B -> K^(i) a tower** (K, K0*(700),
   K0*(1430), K*(892), K*(1410), K*(1680), K1(1270), K1(1400), K2*(1430)),
   for both B+ and B0, with the Boiarska et al. (arXiv:1904.10447) form
   factors as implemented in ALPINIST.  GKOZ find the tower gives ~4x the
   K + K*(892) rate alone.  (Bs -> phi a is not included, matching ALPINIST.)

Leptonic couplings carry the mild RG enhancement C_ll(mu_w)/C(Lambda) = CLL_RG
from the same evaluation.  Decays to gamma gamma and to all-neutral hadronic
final states (3pi0, pi0 pi0 eta(') with eta(') -> neutrals, K0 K0 pi0) produce
no prompt charged tracks at the vertex and are counted in the total width but
not in the visible channels (see visible_channel_weights / visible_fraction).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from particle import Particle

# --------------------------------------------------------------------------
# Physical constants (PDG 2024 central values)
# --------------------------------------------------------------------------
HBAR_C_GEV_M = 1.973269804e-16   # GeV * m  (hbar c)
HBAR_GEV_S = 6.582119569e-25     # GeV * s

ALPHA_S_HAD = 0.30               # fixed alpha_s for the (1+a_s/pi) QCD factor

# Lepton + quark masses (GeV).  Current quark masses for the perturbative width.
M_E = Particle.from_pdgid(11).mass * 1e-3
M_MU = Particle.from_pdgid(13).mass * 1e-3
M_TAU = Particle.from_pdgid(15).mass * 1e-3
M_PION = Particle.from_pdgid(211).mass * 1e-3

QUARK_MASSES = {                 # MS-bar current masses (GeV)
    "u": 0.00216, "d": 0.00467, "s": 0.0934, "c": 1.27, "b": 4.18,
}

# B masses and lifetimes for the B -> K^(i) a production rates.
M_BPLUS = Particle.from_pdgid(521).mass * 1e-3
M_B0 = Particle.from_pdgid(511).mass * 1e-3
TAU_BPLUS_S = 1.638e-12
TAU_B0_S = 1.519e-12

# --------------------------------------------------------------------------
# External-input anchors (provenance: module docstring + data/alpinist/)
# --------------------------------------------------------------------------
# b -> s a coefficient per unit (1/f): ALPINIST above_EW.C_qq_ij_mu at
# Lambda = 1 TeV, universal c_f = 1 (tools/compute_cbs_alpinist.py).
CBS_EFF = 3.518383e-4
# Lepton coupling RG enhancement C_ll(mu_w)/C(Lambda), same evaluation.
CLL_RG = 1.053552

# Reference inverse decay constant at which production weights and ctau are
# tabulated; the (m_a, 1/f) scan rescales off this point (production ~ (1/f)^2,
# ctau ~ (1/f)^-2), so the choice is just a numerical anchor.
INV_F_REF = 1.0e-3               # GeV^-1   (f_ref = 1 TeV)

# Digitized GKOZ (arXiv:2310.03524) decay-width tables via ALPINIST.
ALPINIST_DATA_DIR = Path(__file__).resolve().parent / "data" / "alpinist"

# Charged-decay fractions of the promptly decaying eta/eta' (PDG 2024): the
# pi0 pi0 eta(') hadronic channels are vertex-visible only through these.
ETA_CHARGED_FRAC = 0.272         # eta -> pi+pi-pi0 (0.230) + pi+pi-gamma (0.042)
ETAP_CHARGED_FRAC = 0.78         # eta' -> pi+pi-eta (0.426) + pi+pi-gamma (0.295) + ...

_TABLE_CACHE: dict[str, tuple[np.ndarray, np.ndarray]] = {}


def _width_table(name):
    """(m_a, Gamma/(1/f)^2) columns of a digitized GKOZ table (BNT convention;
    normalisation anchors in the module docstring)."""
    if name not in _TABLE_CACHE:
        arr = np.loadtxt(ALPINIST_DATA_DIR / f"digitized_width_{name}.txt")
        _TABLE_CACHE[name] = (arr[:, 0], arr[:, 1])
    return _TABLE_CACHE[name]


def _table_width(name, m_a):
    """Interpolated Gamma/(1/f)^2 [GeV^3]; zero outside the tabulated range."""
    m_grid, g_grid = _width_table(name)
    m_a = float(m_a)
    if m_a < m_grid[0] or m_a > m_grid[-1]:
        return 0.0
    return float(np.interp(m_a, m_grid, g_grid))


def table_mass_max(name="TotalHad"):
    return float(_width_table(name)[0][-1])


# --------------------------------------------------------------------------
# a -> f fbar partial widths, total width, lifetime
# --------------------------------------------------------------------------

def _two_body_fermion_width(m_a, m_f, n_c, inv_f, c_f=1.0):
    """Gamma(a -> f fbar) [GeV] for the pseudoscalar coupling g = c_f m_f / f."""
    m_a = float(m_a)
    if m_a <= 2.0 * m_f:
        return 0.0
    beta = np.sqrt(1.0 - (2.0 * m_f / m_a) ** 2)
    return n_c * (c_f * m_f) ** 2 * m_a * inv_f ** 2 * beta / (8.0 * np.pi)


def _perturbative_hadronic_width(m_a, inv_f, c_f=1.0):
    """Quark-level sum (N_c = 3, 1 + alpha_s/pi)."""
    qcd = 1.0 + ALPHA_S_HAD / np.pi
    w = 0.0
    for q in ("u", "d", "s", "c"):
        w += _two_body_fermion_width(m_a, QUARK_MASSES[q], 3, inv_f, c_f)
    return w * qcd


_SEAM_NORM = None


def _matched_perturbative_hadronic_width(m_a, inv_f, c_f=1.0):
    """Perturbative continuation above the GKOZ table ceiling, normalised to
    the data-driven table at the matching point (the raw quark-level sum is
    ~1.5x the table there -- charm mass-scheme sensitivity; the shape is kept,
    the normalisation is anchored to the data-driven side of the seam)."""
    global _SEAM_NORM
    if _SEAM_NORM is None:
        m_seam = table_mass_max()
        w_tab = _table_width("TotalHad", m_seam) * INV_F_REF ** 2
        w_pqcd = _perturbative_hadronic_width(m_seam, INV_F_REF)
        _SEAM_NORM = w_tab / w_pqcd if w_pqcd > 0.0 else 1.0
    return _SEAM_NORM * _perturbative_hadronic_width(m_a, inv_f, c_f)


def alp_partial_widths(m_a, inv_f, c_f=1.0):
    """Per-channel partial widths [GeV] as a dict.

    Leptonic channels (ee, mumu, tautau) are analytic with the CLL_RG running
    factor.  ``hadronic`` (total, incl. all-neutral modes) and ``gammagamma``
    are the data-driven GKOZ tables; above the table ceiling ``hadronic``
    falls back to the perturbative quark-level sum and ``gammagamma`` to zero
    (negligible there).
    """
    g_lep = CLL_RG * c_f
    w = {
        "ee": _two_body_fermion_width(m_a, M_E, 1, inv_f, g_lep),
        "mumu": _two_body_fermion_width(m_a, M_MU, 1, inv_f, g_lep),
        "tautau": _two_body_fermion_width(m_a, M_TAU, 1, inv_f, g_lep),
    }
    if m_a <= table_mass_max():
        w["hadronic"] = _table_width("TotalHad", m_a) * inv_f ** 2 * c_f ** 2
        w["gammagamma"] = _table_width("2Gamma", m_a) * inv_f ** 2 * c_f ** 2
    else:
        w["hadronic"] = _matched_perturbative_hadronic_width(m_a, inv_f, c_f)
        w["gammagamma"] = 0.0
    return w


def hadronic_invisible_width(m_a, inv_f, c_f=1.0):
    """Width into hadronic final states with no prompt charged track at the
    decay vertex: 3pi0 and K0 K0bar pi0 entirely, pi0 pi0 eta(') for the
    neutral eta(') decay fractions.  (K_S/K_L decay centimetres-to-metres
    downstream, not at the ALP vertex, so K0 modes give no vertex tracks.)
    Zero above the table ceiling (perturbative region: quark final states)."""
    if m_a > table_mass_max():
        return 0.0
    w = (_table_width("3Pi0", m_a)
         + _table_width("2K0Pi0", m_a)
         + (1.0 - ETA_CHARGED_FRAC) * _table_width("2Pi0Eta", m_a)
         + (1.0 - ETAP_CHARGED_FRAC) * _table_width("2Pi0EtaPrim", m_a))
    return w * inv_f ** 2 * c_f ** 2


def alp_total_width(m_a, inv_f, c_f=1.0):
    return float(sum(alp_partial_widths(m_a, inv_f, c_f).values()))


def alp_ctau(m_a, inv_f, c_f=1.0):
    """Lab-frame c*tau [m] of the ALP at rest-frame width Gamma_tot."""
    gamma = alp_total_width(m_a, inv_f, c_f)
    if gamma <= 0.0:
        return np.inf
    return HBAR_C_GEV_M / gamma


def alp_branchings(m_a, inv_f=INV_F_REF, c_f=1.0):
    """Branching ratios per channel (coupling-independent ratios)."""
    w = alp_partial_widths(m_a, inv_f, c_f)
    tot = sum(w.values())
    if tot <= 0.0:
        return {k: 0.0 for k in w}
    return {k: v / tot for k, v in w.items()}


# --------------------------------------------------------------------------
# Visible decay channels for the charged-track reconstruction
# --------------------------------------------------------------------------
# Each channel is reconstructed as a two-charged-track final state (the best-two
# reduction in the acceptance MC selects exactly two tracks).  For multi-prong
# modes (tautau, hadronic) the two charged daughters are the proxy for the two
# leading charged tracks.  The tau decay length (c*tau ~ 87 um * beta*gamma) is
# far below the 3 mm hit resolution, so each tau is treated as one charged track
# along its flight direction.  gammagamma and the all-neutral hadronic modes
# leave no vertex tracks: they enter the total width (lifetime) but not the
# visible weights.
#
# (pdg-pair used for the two charged daughters, daughter mass).
VISIBLE_CHANNELS = {
    "ee": (11, M_E),
    "mumu": (13, M_MU),
    "tautau": (15, M_TAU),
    "hadronic": (211, M_PION),   # proxy: two leading charged hadrons ~ pions
}


def visible_channel_weights(m_a, inv_f=INV_F_REF, c_f=1.0):
    """{channel: BR} over the track-producing channels (ratios coupling
    independent).  ``hadronic`` here is the *charged-visible* hadronic BR:
    total hadronic minus the all-neutral modes."""
    w = alp_partial_widths(m_a, inv_f, c_f)
    tot = sum(w.values())
    if tot <= 0.0:
        return {}
    w_vis_had = max(w["hadronic"] - hadronic_invisible_width(m_a, inv_f, c_f), 0.0)
    out = {}
    for ch in ("ee", "mumu", "tautau"):
        if w[ch] > 0.0:
            out[ch] = w[ch] / tot
    if w_vis_had > 0.0:
        out["hadronic"] = w_vis_had / tot
    return out


def visible_fraction(m_a, inv_f=INV_F_REF, c_f=1.0):
    """Total BR into channels with >= 2 prompt charged tracks (multiplies the
    production yield in the sensitivity scan; coupling independent)."""
    return float(sum(visible_channel_weights(m_a, inv_f, c_f).values()))


# --------------------------------------------------------------------------
# Production: BR(B -> K^(i) a) over the kaon tower (b -> s a penguin)
# --------------------------------------------------------------------------
# Amplitudes as in ALPINIST alp_production_BRs/ds_production_BRs (BSD-3;
# form factors from Boiarska et al. arXiv:1904.10447):
#
#   Gamma(B -> K_i a) = |g_bs|^2 m_bs2(K_i) |M_i(m_a)|^2
#                        * sqrt(lambda(m_B^2, m_K^2, m_a^2)) / (16 pi m_B^3)
#
# with g_bs = CBS_EFF * (1/f), m_bs2 the (m_b -+ m_s)^2 quark-density factor
# (pseudoscalar/axial daughters: -, scalar/vector/tensor: +).

M_B_QUARK = QUARK_MASSES["b"]
M_S_QUARK = QUARK_MASSES["s"]

# Kaon-tower masses (GeV): charged partner for B+, neutral partner for B0
# where they differ (ALPINIST exotic_constants m_K dict / PDG).
KAON_TOWER = {
    #  name          m(charged)  m(neutral)
    "K":            (0.493677, 0.497611),
    "K0star_700":   (0.845, 0.845),
    "K0star_1430":  (1.425, 1.425),
    "Kstar_892":    (0.89176, 0.89555),
    "Kstar_1410":   (1.414, 1.414),
    "Kstar_1680":   (1.718, 1.718),
    "K1_1270":      (1.253, 1.253),
    "K1_1400":      (1.403, 1.403),
    "K2star_1430":  (1.4273, 1.4324),
}


def _kallen(a, b, c):
    return a * a + b * b + c * c - 2.0 * (a * b + b * c + c * a)


def _lambda_kallen_m(m1, m2, m3):
    """Kallen lambda of the squared masses, lambda(m1^2, m2^2, m3^2)."""
    return _kallen(m1 * m1, m2 * m2, m3 * m3)


def g_bs(inv_f):
    """Flavour-violating b-s-a coupling [GeV^-1] at inverse decay constant
    1/f (linear rescale off the Lambda = 1 TeV evaluation of CBS_EFF)."""
    return CBS_EFF * inv_f


def _M_BP(m_B, m_K, q, name):
    """B -> K (pseudoscalar) via the scalar density, f_0 form factor."""
    f_0, m_fit = 0.33, 6.12
    return 0.5 * (m_B ** 2 - m_K ** 2) / (M_B_QUARK - M_S_QUARK) \
        * f_0 / (1.0 - q ** 2 / m_fit ** 2)


def _M_BS(m_B, m_K, q, name):
    """B -> K0* (scalar) via the pseudoscalar density."""
    if name == "K0star_700":
        f_p, a, b = 0.46, 1.6, 1.35
    else:  # K0star_1430
        f_p, a, b = 0.17, 4.4, 6.4
    return 0.5 * (m_B ** 2 - m_K ** 2 - q ** 2) / (M_B_QUARK + M_S_QUARK) \
        * f_p / (1.0 - a * q ** 2 / m_B ** 2 + b * (q ** 2 / m_B ** 2) ** 2)


def _M_BV(m_B, m_K, q, name):
    """B -> K* (vector), A_0 form factor."""
    if name == "Kstar_892":
        a_0 = 1.364 / (1.0 - q ** 2 / m_B ** 2) - 0.990 / (1.0 - q ** 2 / 36.78)
    elif name == "Kstar_1410":
        a_0 = ((1.0 - 2.0 * m_K ** 2 / (m_B ** 2 + m_K ** 2 - q ** 2)) * 0.22
               + m_K / m_B * 0.28) / (1.0 - q ** 2 / m_B ** 2)
    else:  # Kstar_1680
        a_0 = ((1.0 - 2.0 * m_K ** 2 / (m_B ** 2 + m_K ** 2 - q ** 2)) * 0.18
               + m_K / m_B * 0.24) / (1.0 - q ** 2 / m_B ** 2)
    return -0.5 * np.sqrt(max(_lambda_kallen_m(m_B, m_K, q), 0.0)) \
        / (M_B_QUARK + M_S_QUARK) * a_0


def _M_BA(m_B, m_K, q, name):
    """B -> K1 (axial vector)."""
    th = -0.593

    def v_ab(f0, a, b):
        return f0 / (1.0 - a * (q / m_B) ** 2 + b * (q / m_B) ** 4)

    if name == "K1_1270":
        v_0 = np.sin(th) * 1.31 * v_ab(0.22, 2.4, 1.78) \
            + np.cos(th) * 1.34 * v_ab(-0.45, 1.34, 0.64)
    else:  # K1_1400
        v_0 = np.cos(th) * 1.31 * v_ab(0.22, 2.4, 1.78) \
            - np.sin(th) * 1.34 * v_ab(-0.45, 1.34, 0.64)
    return 0.5 * np.sqrt(max(_lambda_kallen_m(m_B, m_K, q), 0.0)) \
        / (M_B_QUARK - M_S_QUARK) * v_0 / m_K


def _M_BT(m_B, m_K, q, name):
    """B -> K2* (tensor)."""
    a_t = 0.23 / ((1.0 - (q / m_B) ** 2)
                  * (1.0 - 1.23 * (q / m_B) ** 2 + 0.76 * (q / m_B) ** 4))
    return -0.5 * np.sqrt(1.0 / 6.0) * (_lambda_kallen_m(m_B, m_K, q) / (m_B * m_K)) \
        / (M_B_QUARK + M_S_QUARK) * a_t


# name -> (matrix element, quark-density sign: -1 => (m_b - m_s), +1 => (m_b + m_s))
_TOWER_AMPLITUDES = {
    "K":            (_M_BP, -1),
    "K0star_700":   (_M_BS, +1),
    "K0star_1430":  (_M_BS, +1),
    "Kstar_892":    (_M_BV, +1),
    "Kstar_1410":   (_M_BV, +1),
    "Kstar_1680":   (_M_BV, +1),
    "K1_1270":      (_M_BA, -1),
    "K1_1400":      (_M_BA, -1),
    "K2star_1430":  (_M_BT, +1),
}


def _parent_props(parent):
    if parent == "B+":
        return M_BPLUS, TAU_BPLUS_S, 0
    if parent == "B0":
        return M_B0, TAU_B0_S, 1
    raise ValueError(f"unknown parent {parent!r}")


def kaon_mass(kaon, parent="B+"):
    _, _, col = _parent_props(parent)
    return KAON_TOWER[kaon][col]


def width_B_to_Ka(m_a, inv_f, parent="B+", kaon="K"):
    """Gamma(B -> K_i a) [GeV] for one kaon-tower channel."""
    m_a = float(m_a)
    m_B, _, _ = _parent_props(parent)
    m_K = kaon_mass(kaon, parent)
    if m_a >= m_B - m_K:
        return 0.0
    lam = _lambda_kallen_m(m_B, m_K, m_a)
    if lam <= 0.0:
        return 0.0
    amp_fn, sign = _TOWER_AMPLITUDES[kaon]
    m_q = M_B_QUARK + sign * M_S_QUARK
    amp2 = (g_bs(inv_f) * m_q * amp_fn(m_B, m_K, m_a, kaon)) ** 2
    return amp2 * np.sqrt(lam) / (16.0 * np.pi * m_B ** 3)


def br_B_to_Ka(m_a, inv_f, parent="B+", kaon="K"):
    """Branching ratio of one B -> K_i a channel."""
    _, tau_s, _ = _parent_props(parent)
    return width_B_to_Ka(m_a, inv_f, parent, kaon) / (HBAR_GEV_S / tau_s)


def br_B_to_K_a(m_a, inv_f, parent="B+"):
    """BR of the ground-state B -> K a channel (kept for tests/back-compat)."""
    return br_B_to_Ka(m_a, inv_f, parent, "K")


def br_B_tower_total(m_a, inv_f, parent="B+"):
    """Sum of BR(B -> K_i a) over the whole kaon tower."""
    return sum(br_B_to_Ka(m_a, inv_f, parent, k) for k in KAON_TOWER)


# Production channels: (parent label, FONLL species pdg) x kaon tower.
# Fragmentation fractions reuse the HNL convention (FRAG_B).
PRODUCTION_PARENTS = (
    ("B+", 521),
    ("B0", 511),
)
