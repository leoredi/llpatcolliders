r"""BC10 fermiophilic-ALP model layer (couplings -> production, lifetime, BRs).

PBC benchmark **BC10**: a pseudoscalar axion-like particle ``a`` whose dominant
coupling is to SM fermions (the genuine charged-track "ALP" benchmark). The
single parameter that this analysis sets a limit on is the inverse decay
constant ``1/f`` (with universal Wilson coefficients ``c_f = 1``); it controls
production, lifetime, and the visible branching ratios *simultaneously*.

Conventions and formulae (all referenced, none invented):

Lagrangian (Bauer, Neubert, Thamm, "Collider probes of axion-like particles",
JHEP 12 (2017) 044, arXiv:1708.00443; and the EFT review Bauer et al.,
arXiv:2012.12272):

    L_int = (d_mu a / f) sum_f c_f  fbar gamma^mu gamma_5 f .

Integrating by parts and using the axial-current divergence
d_mu(fbar gamma^mu gamma_5 f) = 2 i m_f fbar gamma_5 f gives an effective
pseudoscalar Yukawa with coupling g_aff = c_f m_f / f, and the partial width

    Gamma(a -> f fbar) = N_c c_f^2 m_f^2 m_a / (8 pi f^2) * sqrt(1 - 4 m_f^2/m_a^2)

(BNT 2017 Eq. 2.13; N_c = 3 for quarks, 1 for leptons).  This is the source the
2025 PBC report (arXiv:2505.00947) and the unified FIP calculation
(arXiv:2311.00507) use for the BC10 decay table.

Flavour-violating production b -> s a (top-W penguin, leading log; BNT 2017 and
"Flavour probes of ALPs", arXiv:2110.10698):

    L_FCNC = (d_mu a / f) h_sb  sbar_L gamma^mu b_L + h.c. ,
    h_sb   = C_TOP/(16 pi^2) * (m_t^2 / v^2) * V_tb V_ts^* * log(Lambda_UV^2/m_t^2)

with the ALP-photon/loop O(1) coefficient C_TOP = 1 for the universal c_f = 1
case (the top, c_t = 1, dominates the penguin).  The two-body width of the
pseudoscalar via the vector current (only the q_mu f_0 piece survives between
two pseudoscalars):

    Gamma(B -> K a) = |p_a|/(8 pi m_B^2) * (h_sb/f)^2
                      * (m_B^2 - m_K^2)^2 * f_0(m_a^2)^2 ,
    |p_a| = lambda^{1/2}(m_B^2, m_K^2, m_a^2) / (2 m_B) .

CAVEATS (see EXTERNAL_INPUTS_NEEDED.md):
  * The absolute b->s a normalization (C_TOP, Lambda_UV, and any finite
    matching piece) and the B->K scalar form factor f_0(q^2) should be matched
    to the BC10 tabulation in arXiv:2311.00507 / ALPINIST (arXiv:2105.10806).
    The reach scales as 1/f ~ C_TOP^{-1} f_0^{-1}; the island below is reported
    with the leading-log value and that scaling is stated.
  * The hadronic width in the 2 m_pi -- ~2 GeV window is genuinely data-driven
    (ALP-pi/eta/eta' mixing). Here it is approximated perturbatively (quark
    level, current masses, 1+alpha_s/pi); this is documented as a placeholder.
"""

from __future__ import annotations

import numpy as np
from particle import Particle

# --------------------------------------------------------------------------
# Physical constants (PDG 2024 central values)
# --------------------------------------------------------------------------
HBAR_C_GEV_M = 1.973269804e-16   # GeV * m  (hbar c)
HBAR_GEV_S = 6.582119569e-25     # GeV * s

V_HIGGS = 246.21965              # GeV   (electroweak vev)
M_TOP = 172.5                    # GeV   (top pole mass)
V_TB_VTS = 0.0405                # |V_tb V_ts^*|  (PDG CKM fit)

ALPHA_S_HAD = 0.30               # fixed alpha_s for the (1+a_s/pi) QCD factor

# Lepton + quark masses (GeV).  Current quark masses for the perturbative width.
M_E = Particle.from_pdgid(11).mass * 1e-3
M_MU = Particle.from_pdgid(13).mass * 1e-3
M_TAU = Particle.from_pdgid(15).mass * 1e-3
M_PION = Particle.from_pdgid(211).mass * 1e-3

QUARK_MASSES = {                 # MS-bar current masses (GeV)
    "u": 0.00216, "d": 0.00467, "s": 0.0934, "c": 1.27, "b": 4.18,
}

# B / K masses and lifetimes for the B -> K a production rate.
M_BPLUS = Particle.from_pdgid(521).mass * 1e-3
M_B0 = Particle.from_pdgid(511).mass * 1e-3
M_KPLUS = Particle.from_pdgid(321).mass * 1e-3
M_K0 = Particle.from_pdgid(311).mass * 1e-3
TAU_BPLUS_S = 1.638e-12
TAU_B0_S = 1.519e-12

# B -> K scalar form factor f_0(q^2): f_0(0) with a single-pole rise tuned to
# the lattice value f_0(q^2_max) ~ 0.8 (HPQCD/HFLAV B->K).  See EXTERNAL_INPUTS.
F0_B_K_AT_0 = 0.330
F0_B_K_POLE2 = 37.5              # GeV^2

# Penguin O(1) coefficient and UV scale entering the leading log.
C_TOP = 1.0
LAMBDA_UV_GEV = 1000.0           # 1 TeV reference UV scale (fixes the log)

# Reference inverse decay constant at which production weights and ctau are
# tabulated; the (m_a, 1/f) scan rescales off this point (production ~ (1/f)^2,
# ctau ~ (1/f)^-2), so the choice is just a numerical anchor.
INV_F_REF = 1.0e-3               # GeV^-1   (f_ref = 1 TeV)


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


def alp_partial_widths(m_a, inv_f, c_f=1.0):
    """Per-channel partial widths [GeV] as a dict.

    Leptonic channels (ee, mumu, tautau) are exact.  ``hadronic`` is the
    perturbative quark-level sum over open q = u,d,s,c with N_c = 3 and a
    (1 + alpha_s/pi) factor, gated at m_a > 2 m_pi (no open hadronic phase space
    below).  This is a documented placeholder for the data-driven spectral
    function in the resonance region (see module docstring / EXTERNAL_INPUTS).
    """
    w = {
        "ee": _two_body_fermion_width(m_a, M_E, 1, inv_f, c_f),
        "mumu": _two_body_fermion_width(m_a, M_MU, 1, inv_f, c_f),
        "tautau": _two_body_fermion_width(m_a, M_TAU, 1, inv_f, c_f),
    }
    w_had = 0.0
    if m_a > 2.0 * M_PION:
        qcd = 1.0 + ALPHA_S_HAD / np.pi
        for q in ("u", "d", "s", "c"):
            w_had += _two_body_fermion_width(m_a, QUARK_MASSES[q], 3, inv_f, c_f)
        w_had *= qcd
    w["hadronic"] = w_had
    return w


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
# along its flight direction.  a -> gamma gamma (loop-induced, O(alpha/4pi)
# relative to the open leptonic mode) is omitted; it is negligible for
# m_a > 2 m_mu and only mildly lengthens ctau.
#
# (pdg-pair used for the two charged daughters, daughter mass).
VISIBLE_CHANNELS = {
    "ee": (11, M_E),
    "mumu": (13, M_MU),
    "tautau": (15, M_TAU),
    "hadronic": (211, M_PION),   # proxy: two leading charged hadrons ~ pions
}


def visible_channel_weights(m_a, inv_f=INV_F_REF, c_f=1.0):
    """Return {channel: BR} over the track-producing channels (ratios are
    coupling independent).  All BC10 channels in the B->K a window are charged,
    so the tracker-visible fraction is ~1 (gamma gamma omitted)."""
    br = alp_branchings(m_a, inv_f, c_f)
    return {k: br[k] for k in VISIBLE_CHANNELS if br.get(k, 0.0) > 0.0}


# --------------------------------------------------------------------------
# Production: BR(B -> K a) via the b -> s a top penguin
# --------------------------------------------------------------------------

def _kallen(a, b, c):
    return a * a + b * b + c * c - 2.0 * (a * b + b * c + c * a)


def h_sb(inv_f, c_top=C_TOP, lambda_uv=LAMBDA_UV_GEV):
    """Dimensionless flavour-violating coefficient of (d_mu a/f) sbar_L g^mu b_L.

    Leading-log top penguin; the (1/f) is carried separately (amplitude uses
    h_sb * inv_f).  See module docstring for the reference and caveats.
    """
    log = np.log(lambda_uv ** 2 / M_TOP ** 2)
    return c_top / (16.0 * np.pi ** 2) * (M_TOP ** 2 / V_HIGGS ** 2) * V_TB_VTS * log


def f0_B_K(q2):
    """B -> K scalar form factor f_0(q^2) (single-pole, lattice-anchored)."""
    return F0_B_K_AT_0 / (1.0 - q2 / F0_B_K_POLE2)


def width_B_to_K_a(m_a, inv_f, m_B, m_K, c_top=C_TOP, lambda_uv=LAMBDA_UV_GEV):
    """Gamma(B -> K a) [GeV]."""
    m_a = float(m_a)
    if m_a >= m_B - m_K:
        return 0.0
    lam = _kallen(m_B ** 2, m_K ** 2, m_a ** 2)
    if lam <= 0.0:
        return 0.0
    p_a = np.sqrt(lam) / (2.0 * m_B)
    coupling = h_sb(inv_f, c_top, lambda_uv) * inv_f          # GeV^-1
    amp2 = (coupling * (m_B ** 2 - m_K ** 2) * f0_B_K(m_a ** 2)) ** 2  # GeV^2
    return p_a / (8.0 * np.pi * m_B ** 2) * amp2


def br_B_to_K_a(m_a, inv_f, parent="B+", c_top=C_TOP, lambda_uv=LAMBDA_UV_GEV):
    """Branching ratio BR(B -> K a) for parent 'B+' (K+) or 'B0' (K0)."""
    if parent == "B+":
        m_B, m_K, tau_s = M_BPLUS, M_KPLUS, TAU_BPLUS_S
    elif parent == "B0":
        m_B, m_K, tau_s = M_B0, M_K0, TAU_B0_S
    else:
        raise ValueError(f"unknown parent {parent!r}")
    gamma = width_B_to_K_a(m_a, inv_f, m_B, m_K, c_top, lambda_uv)
    gamma_B = HBAR_GEV_S / tau_s
    return gamma / gamma_B


# Parents used for production (pseudoscalar B -> K a channels).  Fragmentation
# fractions reuse the HNL convention (FRAG_B); B_s -> phi a (and the B -> K* a
# vector channels) would add an O(1) uplift and are flagged as an extension.
PRODUCTION_PARENTS = (
    # (parent label, FONLL species pdg, K daughter mass, B mass)
    ("B+", 521, M_KPLUS, M_BPLUS),
    ("B0", 511, M_K0, M_B0),
)
