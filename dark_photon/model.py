"""BC1 dark-photon model layer: widths, branching ratios, and lifetime.

Minimal vector portal (PBC benchmark BC1, Beacham:2019nyx): a massive dark
photon ``A'`` with kinetic mixing ``epsilon`` to the SM photon. The single
parameter ``epsilon^2`` fixes both the production rate (each production BR
scales as ``epsilon^2``) and the lifetime (total width scales as
``epsilon^2``), so at fixed mass the excluded region is a closed island in
``(m_A', epsilon^2)`` -- the same single-coupling structure as BC4/BC10.

Physics:
  * Leptonic widths are the exact tree-level vector-current result
    ``Gamma(A'->l+l-) = (alpha/3) eps^2 m_A' sqrt(1-4r)(1+2r)``, ``r=m_l^2/m_A'^2``.
  * Hadronic width uses the measured R-ratio,
    ``Gamma(A'->had) = Gamma(A'->mu+mu-) x R(m_A')``. Exclusive branching
    fractions for ``0.2 <= m_A' <= 1.7`` GeV are taken from the DeLiVeR
    vector-meson-dominance tabulation (arXiv:2201.01788), committed at
    ``data/br_tables/dp_brs_deliver.csv``; above 1.7 GeV a perturbative-QCD
    R-ratio with one-loop alpha_s running is used.

Coupling convention: all widths and production BRs below are quoted at
``epsilon^2 = 1``; the sensitivity scan multiplies by the physical
``epsilon^2``. Ported from the reviewed vector-portal package on the
``DP`` development branch (``dark_photon/generator/dp_meson_brs.py``).
"""
from __future__ import annotations

import math
from functools import lru_cache
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
ALPHA_EM = 1.0 / 137.036            # fine-structure constant
HBAR_C_GEV_M = 0.197326980e-15      # hbar c in GeV m

# Particle masses (GeV), PDG 2024
M_ELECTRON = 0.000511
M_MUON = 0.10566
M_TAU = 1.77686
M_PIPLUS = 0.13957
M_PI0 = 0.13498
M_KPLUS = 0.49368
M_K0 = 0.49761
M_PROTON = 0.93827
M_ETA = 0.54785
M_OMEGA = 0.78266

# PDG parent branching ratios used for meson production
BR_PI0_GAMGAM = 0.98823
BR_ETA_GAMGAM = 0.3941
BR_OMEGA_PI0GAM = 0.0828

# Kinematic ceiling of the tabulated (VMD) decay description.
M_A_MAX_VMD = 1.70                  # GeV

_BR_TABLE = Path(__file__).resolve().parent / "data" / "br_tables" / "dp_brs_deliver.csv"

# DeLiVeR exclusive-channel columns, mapped to a stable internal order.
_DELIVER_CHANNELS = (
    "ee", "mumu", "tau", "pipi", "pi3", "pi4c", "pi4n",
    "KKc", "KKn", "PiGamma", "EtaGamma", "ppbar", "had_other",
)


# ---------------------------------------------------------------------------
# DeLiVeR VMD branching-ratio table
# ---------------------------------------------------------------------------
@lru_cache(maxsize=1)
def _load_deliver():
    """Return ``(masses, {channel: br_array})`` from the committed VMD table."""
    cols: list[str] = []
    rows: list[list[float]] = []
    with open(_BR_TABLE) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if not cols:
                cols = [c.strip() for c in line.split(",")]
                continue
            rows.append([float(x) for x in line.split(",")])
    data = np.asarray(rows)
    idx = {c: i for i, c in enumerate(cols)}
    masses = data[:, idx["mass_GeV"]]
    table = {ch: data[:, idx[ch]] for ch in _DELIVER_CHANNELS if ch in idx}
    return masses, table


# ---------------------------------------------------------------------------
# Partial widths (at epsilon^2 = 1)
# ---------------------------------------------------------------------------
def _beta_factor(m_a: float, m_l: float) -> float:
    if m_a <= 2.0 * m_l:
        return 0.0
    r = (m_l / m_a) ** 2
    return math.sqrt(1.0 - 4.0 * r) * (1.0 + 2.0 * r)


def width_leptonic(m_a: float, m_l: float) -> float:
    """Gamma(A' -> l+ l-) at eps^2 = 1 [GeV]."""
    return (ALPHA_EM / 3.0) * m_a * _beta_factor(m_a, m_l)


def r_ratio(m_a: float) -> float:
    """R(s) = sigma(e+e- -> had)/sigma(e+e- -> mu+mu-) at sqrt(s) = m_A'."""
    if m_a < 2.0 * M_PIPLUS:
        return 0.0
    if m_a <= M_A_MAX_VMD:
        masses, table = _load_deliver()
        br_mumu = float(np.interp(m_a, masses, table["mumu"]))
        if br_mumu <= 0:
            return 0.0
        br_had = float(np.interp(m_a, masses, table.get(
            "had_other", np.zeros_like(masses))))
        # BRqcd is the summed hadronic fraction; reconstruct from 1 - leptons.
        br_ee = float(np.interp(m_a, masses, table["ee"]))
        br_tau = float(np.interp(m_a, masses, table["tau"]))
        br_had_total = max(0.0, 1.0 - br_ee - br_mumu - br_tau)
        return br_had_total / br_mumu
    return _r_ratio_pqcd(m_a)


def _r_ratio_pqcd(m_a: float) -> float:
    """Perturbative R-ratio with one-loop alpha_s (m_A' > 1.7 GeV)."""
    lam = 0.2
    if m_a <= lam:
        return 0.0
    nf = 3 + (m_a > 2 * 1.27) + (m_a > 2 * 4.18)
    als = min(0.4, 12 * math.pi / ((33 - 2 * nf) * math.log((m_a / lam) ** 2)))
    eq2 = (2.0 / 3.0) ** 2 + 2 * (1.0 / 3.0) ** 2      # u, d, s
    if m_a > 2 * 1.27:
        eq2 += (2.0 / 3.0) ** 2                         # c
    if m_a > 2 * 4.18:
        eq2 += (1.0 / 3.0) ** 2                         # b
    return 3.0 * eq2 * (1.0 + als / math.pi)


def width_hadronic(m_a: float) -> float:
    """Gamma(A' -> hadrons) at eps^2 = 1 [GeV] via the R-ratio."""
    return width_leptonic(m_a, M_MUON) * r_ratio(m_a)


def partial_widths(m_a: float) -> dict:
    """All partial widths at eps^2 = 1 [GeV]."""
    return {
        "ee": width_leptonic(m_a, M_ELECTRON),
        "mumu": width_leptonic(m_a, M_MUON),
        "tautau": width_leptonic(m_a, M_TAU),
        "hadrons": width_hadronic(m_a),
    }


def total_width(m_a: float) -> float:
    """Total A' width at eps^2 = 1 [GeV]."""
    return sum(partial_widths(m_a).values())


# ---------------------------------------------------------------------------
# Decay branching ratios (visible-mode composition for the acceptance MC)
# ---------------------------------------------------------------------------
def branching_ratios(m_a: float) -> dict:
    """Decay branching ratios of the A', normalized to sum to 1.

    For ``0.2 <= m_A' <= 1.7`` GeV the exclusive VMD channels of the DeLiVeR
    table are returned (ee, mumu, tau, and the resolved hadronic exclusive
    modes). Outside that window the composition is leptonic widths plus an
    inclusive ``hadrons`` mode from the R-ratio.
    """
    if 0.2 <= m_a <= M_A_MAX_VMD:
        masses, table = _load_deliver()
        br = {ch: float(np.interp(m_a, masses, table[ch]))
              for ch in table}
        br = {k: v for k, v in br.items() if v > 0.0}
        s = sum(br.values())
        return {k: v / s for k, v in br.items()} if s > 0 else br
    w = partial_widths(m_a)
    tot = sum(w.values())
    if tot <= 0:
        return {}
    return {k: v / tot for k, v in w.items() if v > 0.0}


# ---------------------------------------------------------------------------
# Lifetime
# ---------------------------------------------------------------------------
def ctau(m_a: float, eps2: float) -> float:
    """A' proper decay length c*tau [m] at kinetic mixing eps^2."""
    g0 = total_width(m_a)
    if g0 <= 0 or eps2 <= 0:
        return float("inf")
    return HBAR_C_GEV_M / (eps2 * g0)


def ctau_eps2_1(m_a: float) -> float:
    """c*tau [m] at eps^2 = 1; the scan rescales as c*tau ~ 1/eps^2."""
    return ctau(m_a, 1.0)


# ---------------------------------------------------------------------------
# Meson-production branching ratios (at epsilon^2 = 1)
# ---------------------------------------------------------------------------
def _kallen(a: float, b: float, c: float) -> float:
    return a * a + b * b + c * c - 2 * (a * b + a * c + b * c)


def br_pseudoscalar_to_dp_gamma(m_a: float, m_P: float, br_P_gamgam: float) -> float:
    """BR(P -> A' gamma) at eps^2 = 1 for a pseudoscalar P (pi0, eta).

    ``BR(P->A'gamma) = 2 eps^2 BR(P->gamma gamma) (1 - m_A'^2/m_P^2)^3``
    (arXiv:2005.01515). Transition form factor set to unity (< few % for
    m_A' well below m_P). Returns the eps^2 = 1 value.
    """
    if m_a >= m_P:
        return 0.0
    return 2.0 * br_P_gamgam * (1.0 - (m_a / m_P) ** 2) ** 3


def br_omega_to_dp_pi0(m_a: float) -> float:
    """BR(omega -> A' pi0) at eps^2 = 1 (arXiv:2005.01515)."""
    if m_a >= M_OMEGA - M_PI0:
        return 0.0
    p_gamma = (M_OMEGA ** 2 - M_PI0 ** 2) / (2.0 * M_OMEGA)
    lam = _kallen(M_OMEGA ** 2, M_PI0 ** 2, m_a ** 2)
    if lam <= 0:
        return 0.0
    p_dp = math.sqrt(lam) / (2.0 * M_OMEGA)
    return BR_OMEGA_PI0GAM * (p_dp / p_gamma) ** 3


# Production parent registry: label -> (parent mass, production-BR fn at eps^2=1,
# recoil-partner mass in the two-body meson decay P/V -> A' + X).
PRODUCTION_PARENTS = {
    "pi0": (M_PI0,
            lambda m: br_pseudoscalar_to_dp_gamma(m, M_PI0, BR_PI0_GAMGAM),
            0.0),                     # pi0 -> A' gamma
    "eta": (M_ETA,
            lambda m: br_pseudoscalar_to_dp_gamma(m, M_ETA, BR_ETA_GAMGAM),
            0.0),                     # eta -> A' gamma
    "omega": (M_OMEGA, br_omega_to_dp_pi0, M_PI0),   # omega -> A' pi0
}
